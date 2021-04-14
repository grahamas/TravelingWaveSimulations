
const DEFAULT_SLOPE_MIN = 1e-4
const DEFAULT_VEL_JITTER = 1.5

### Single wave measurements across time ###

struct SpatiotemporalWaveMeasurements{T}
    velocities::Vector{T}
    slopes::Vector{T}
    maxes::Vector{T}
    distance::T
    duration::T
end
function SpatiotemporalWaveMeasurements(pf::Persistent{T}) where {T}
    SpatiotemporalWaveMeasurements{T}( # FIXME had {T} before
        get_velocities(pf),
        [slope_val(wave) for wave in pf.waveforms],
        [max(wave) for wave in pf.waveforms],
        slope_loc(pf.waveforms[end]) - slope_loc(pf.waveforms[begin]),
        pf.t[end] - pf.t[begin]
       )
end

#######################################################
### Single persistent wave classification functions ###
#######################################################

struct WaveClassifications
    traveling::Bool
    positive_velocity::Bool
    decaying::Bool
    growing::Bool
    oscillating::Bool
    positive_slope::Bool
    sane::Bool
end
function WaveClassifications(pf::Persistent; kwargs...)
    measurements = SpatiotemporalWaveMeasurements(pf)
    WaveClassifications(measurements; kwargs...)
end


######################################
### Whole execution classification ###
######################################
abstract type AbstractExecutionClassifications end
struct ExecutionClassifications <: AbstractExecutionClassifications
    has_propagation::Bool
    has_oscillation::Bool
    farthest_propagation_is_decaying::Bool
    farthest_propagation_is_oscillating::Bool
    persistently_active_near_origin::Bool
    reaches_steady_state::Bool
end 
function Base.show(io::IO, ec::ExecutionClassifications)
    nv_pairs = [(name, getfield(ec, name)) for name in fieldnames(ExecutionClassifications)]
    names, values = zip(nv_pairs...)
    nt = NamedTuple{Tuple(names)}(values)
    Base.show(io, nt)
end


##################################
### Propagation Classification ###
##################################

struct RunningFront{T_VAL, VEL<:Union{T_VAL,Missing}, T_LOC}
    starting_location::T_LOC
    previous_location::T_LOC
    previous_time::T_LOC
    previous_velocity::VEL
    slope::T_VAL
end

RF_vel_missing(::Type{RunningFront{T1,IGN,T2}}) where {T1,T2,IGN} = RunningFront{T1,Missing,T2}
RF_possibilities(::Type{RF}) where {T1,V,T2,RF<:RunningFront{T1, V, T2}} = Union{RF,RF_vel_missing(RF), Nothing}
RF_possibilities(::Type{T}) where {T<:Number} = Union{RunningFront{T,T,T},RunningFront{T,Missing,T},Nothing}


_front_loc(front::AxisArray) = axes_keys(front) |> only |> only
_front_slope(front::AxisArray) = only(front)
predict_location(rf::RunningFront{T,T}, current_time::T) where T = rf.previous_location + rf.previous_velocity * (current_time - rf.previous_time)
predict_location(rf::RunningFront{T,Missing}, current_time) where T = rf.previous_location
calc_jitter(rf::RunningFront{T,T}, current_time::T, const_jitter, vel_jitter) where T = const_jitter + abs(rf.previous_velocity * (current_time - rf.previous_time)) * vel_jitter
calc_jitter(rf::RunningFront{T,Missing}, current_time, const_jitter, vel_jitter) where T = const_jitter
function start_running_front(front::AxisArray, current_time)
    # Must be axisarray with only the slope at its location
    loc = _front_loc(front)
    val = _front_slope(front)
    #@info "Start at $loc"
    RunningFront(loc, loc, current_time, missing, val)
end
function extend_running_front(rf::RunningFront, front::AxisArray, current_time) 
    loc = _front_loc(front)
    val = _front_slope(front)
    #@info "Continue to $(loc): dist of $(rf.starting_location - loc)"
    RunningFront(
        rf.starting_location,
        loc,
        current_time,
        (loc - rf.previous_location) / (current_time - rf.previous_time),
        val
    )
end

function front_continues_rf(front, running_front, current_time, const_jitter, vel_jitter)
    """Return true if front can continue running_front.
        Needs:
            1. slope has same sign
            2. projected location is within const_jitter, vel_jitter
    """
    slope = _front_slope(front)
    rf_slope = running_front.slope
    if approx_sign(rf_slope) != 0 && approx_sign(slope) != approx_sign(rf_slope)
        return false
    end

    loc = _front_loc(front)
    predicted_loc = predict_location(running_front, current_time)
    jitter = calc_jitter(running_front, current_time, const_jitter, vel_jitter)
    if !(predicted_loc - jitter <= loc <= predicted_loc + jitter)
        return false
    end

    return true
end

# Resolve when either:
#    RF1 ... F ... RF2
# or
#    F1 .... RF ... F2
# where F is the location of a front, and RF is a running front.
# so when you have RF1 F RF2, you decide to which RF F belongs.
# If F could only belong to RF1, then update RF to have F and place it in F's location, and unset unresolved_front_idx
# If F could only belong to RF2, then terminate RF1, mark F as unresolved_front_idx, and RF2 as preceding_running_front
# If F could belong to either.... warn and give it to RF1
# If F could belong to neither, then set it as a new RF, unset unresolved_front_idx, and set RF2 as preceding running front
function reconcile_running_front!(running_fronts, new_running_front::Nothing,
                                    frame, unresolved_front_idx,
                                    preceding_running_front,
                                    current_time,
                                    const_jitter, vel_jitter)
    # haven't come across RF2
    return (preceding_running_front, unresolved_front_idx, nothing)    
end
function reconcile_running_front!(running_fronts, 
                                  new_running_front::RunningFront,
                                  frame, unresolved_front_idx,
                                  preceding_running_front::Nothing,
                                  current_time,
                                  const_jitter, vel_jitter)
    # there's no RF1
    return (new_running_front, unresolved_front_idx, nothing)
end
function reconcile_running_front!(running_fronts, 
                                    new_running_front::RunningFront,
                                    frame, unresolved_front_idx::Nothing,
                                    preceding_running_front::RunningFront,
                                    current_time,
                                    const_jitter, vel_jitter)
    # There's no F in RF1 F RF2, so just set RF2 as the new preceding front.
    # (also note there can't be F RF1 RF2... RF1 is just orphaned)
    return (new_running_front, unresolved_front_idx, nothing)
end
function reconcile_running_front!(running_fronts, 
                                    rf2::RunningFront,
                                    frame, front_idx::Int,
                                    rf1::RunningFront,
                                    current_time,
                                    const_jitter, vel_jitter)#::Tuple{Union{RF,RFM,Nothing},Union{Int,Nothing},Union{RF,RFM,Nothing}}
    # Now we have RF1 F RF2  (actually! could be F RF1 RF2, but works anyway)
    @assert rf1 != rf2
    f = frame[[front_idx]]
    f_continues_rf1 = front_continues_rf(f, rf1, current_time, const_jitter, vel_jitter)
    f_continues_rf2 = front_continues_rf(f, rf2, current_time, const_jitter, vel_jitter)
    if f_continues_rf1 && !f_continues_rf2
        running_fronts[front_idx] = extend_running_front(rf1, f, current_time)
        return (rf2, nothing, running_fronts[front_idx])
    elseif !f_continues_rf1 && f_continues_rf2
        return (rf2, front_idx, nothing)
    elseif f_continues_rf1 && f_continues_rf2
        # @warn "Ambiguous front continuation! reconciling RF"
        # @show rf1, rf2, f
        if abs(predict_location(rf1, current_time) - _front_loc(f)) <= abs(predict_location(rf2, current_time) - _front_loc(f))
            running_fronts[front_idx] = extend_running_front(rf1, f, current_time)
            return (rf2, nothing, running_fronts[front_idx])
        else
            running_fronts[front_idx] = extend_running_front(rf2, f, current_time)
            return (nothing, nothing, running_fronts[front_idx])
        end
    else
        running_fronts[front_idx] = start_running_front(f, current_time)
        return (rf2, nothing, running_fronts[front_idx])
    end
end
# When you have F1 RF F2  (or RF F1 F2...)
# If RF could only belong to F1 then add F1 to RF and put it at F1
# If RF could only belong to F2, then set preceding_running_front to RF, and unresolved_front_idx to F2, and put a new RF on F1
# If RF could extend either... warn and give it to F1
# if RF could not belong to either, then unset preceding_running_front, set unresolved_front_idx to F2, and put a new RF on F1

function reconcile_front!(running_fronts, front_idx, 
            frame, unresolved_front_idx::Nothing, 
            preceding_running_front, 
            current_time,
            const_jitter, vel_jitter)
    # don't have F1, so just set unresolved_front_idx to be that next time
    return (preceding_running_front, front_idx, nothing)
end
function reconcile_front!(running_fronts, front_idx, 
        frame, unresolved_front_idx::Int, 
        preceding_running_front::Nothing, 
        current_time,
        const_jitter, vel_jitter)
    # have F1 but no RF, so F1 is orphaned.
    running_fronts[unresolved_front_idx] = start_running_front(frame[[unresolved_front_idx]], current_time)
    return (preceding_running_front, front_idx, running_fronts[unresolved_front_idx])
end
function reconcile_front!(running_fronts, f2_idx, 
        frame, f1_idx::Int, 
        rf::RunningFront, 
        current_time,
        const_jitter, vel_jitter)#::Tuple{RF_possibilities,Union{Int,Nothing},Union{RF,RFM}}
    @assert f1_idx != f2_idx
    f1 = frame[[f1_idx]]
    f2 = frame[[f2_idx]]
    f1_continues_rf = front_continues_rf(f1, rf, current_time, const_jitter, vel_jitter)
    f2_continues_rf = front_continues_rf(f2, rf, current_time, const_jitter, vel_jitter)
    if f1_continues_rf && !f2_continues_rf
        running_fronts[f1_idx] = extend_running_front(rf, f1, current_time)
        return (nothing, f2_idx, running_fronts[f1_idx])
    elseif !f1_continues_rf && f2_continues_rf
        running_fronts[f1_idx] = start_running_front(f1, current_time)
        return (rf, f2_idx, running_fronts[f1_idx])
    elseif f1_continues_rf && f2_continues_rf
        #@warn "Ambiguous front continuation! reconciling F"
        # @show f1, f2, rf
        if abs(predict_location(rf, current_time) - _front_loc(f1)) <= abs(predict_location(rf, current_time) - _front_loc(f2))
            running_fronts[f1_idx] = extend_running_front(rf, f1, current_time)
            return (nothing, f2_idx, running_fronts[f1_idx])
        else
            running_fronts[f1_idx] = start_running_front(f1, current_time)
            running_fronts[f2_idx] = extend_running_front(rf, f2, current_time)
            return (nothing, nothing, running_fronts[f2_idx])
        end
    else
        running_fronts[f1_idx] = start_running_front(f1, current_time)
        return (nothing, f2_idx, running_fronts[f1_idx])
    end
end



function continue_fronts!(running_fronts::AxisVector, 
        current_time,
        dframe::AxisVector,
        frame::AxisVector, 
        return_fn::Function, 
        d1_ghost_op, slope_min, 
        const_jitter, vel_jitter)::Bool
    deriv!(dframe, frame, d1_ghost_op)
    locations = axes_keys(frame) |> only
    prev_slope = dframe[begin]
    preceding_running_front = nothing
    unresolved_front_idx = nothing
    for idx in LinearIndices(locations)[begin+2:end-2]
        this_running_front = running_fronts[idx]
        running_fronts[idx] = nothing
        preceding_running_front, unresolved_front_idx, new_rf = reconcile_running_front!(running_fronts, this_running_front, 
                 dframe, unresolved_front_idx, 
                 preceding_running_front, current_time, const_jitter, vel_jitter)
        return_fn(new_rf) && return true
        p2_slope, prev_slope, slope, next_slope, n2_slope = dframe[idx-2:idx+2]
        if abs(slope) > slope_min &&
                ((p2_slope <= prev_slope < slope >= next_slope >= n2_slope) || 
                (p2_slope >= prev_slope > slope <= next_slope <= n2_slope))
            preceding_running_front, unresolved_front_idx, new_rf = reconcile_front!(running_fronts, idx, 
                    dframe, unresolved_front_idx,
                    preceding_running_front, 
                    current_time,
                    const_jitter, vel_jitter
                )
            return_fn(new_rf) && return true
        end
    end
    if !isnothing(unresolved_front_idx)
        if !isnothing(preceding_running_front) &&
            front_continues_rf(dframe[[unresolved_front_idx]], preceding_running_front, current_time, const_jitter, vel_jitter)
            running_fronts[unresolved_front_idx] = extend_running_front( preceding_running_front, dframe[[unresolved_front_idx]], current_time)
        else
            running_fronts[unresolved_front_idx] = start_running_front(dframe[[unresolved_front_idx]], current_time)
        end
        return return_fn(running_fronts[unresolved_front_idx])
    end
    return false
end


## ASSUMPTIONS:
##  1) fronts cannot "pass" eachother because a) they would deform into a single
##     front first, and 2) dt is sufficiently small
##     explicitly, a front can't move further than the distance between two fronts.
##  thus, a front in one frame must either appear de novo or come from an adjacent front
# continue_fronts!(running_fronts, current_time, frame, return_fn(::RunningFront)::Union{STHG,Nothing})::Union{STHG,Nothing}
# Not really a classification....

function is_propagated(u, t, integrator)
    # return true if is_propagated, or is completely flat
    p = integrator.p
    propagated = continue_fronts!(p.running_fronts, t, p.dframe_cache, population(u,1), p.return_fn, p.d1_ghost_op, p.slope_min, p.const_jitter, p.vel_jitter)
    if propagated
        integrator.p.has_propagation[1] = true 
    end
    flattened = (all(abs.(p.dframe_cache) .< sqrt(eps())) && t > 1.0)  
    return propagated || flattened
end

export MinimalPropagationClassification
struct MinimalPropagationClassification <: AbstractExecutionClassifications
    has_propagation::Bool
end
MinimalPropagationClassification(::AbstractExecution; kwargs...) = error("MinimalPropagationClassification should not be pre-reduced")
function MinimalPropagationClassification(exec::AbstractExecution{T,SIM}) where {T,M,S,IV,ALG,DT,SV_IDX,CB<:Tuple{typeof(is_propagated),<:NamedTuple},GR,SIM<:Simulation{T,M,S,IV,ALG,DT,SV_IDX,CB,GR}}
    # callback saved propagation
    return MinimalPropagationClassification(exec.solution.p.has_propagation[1])
end
MinimalPropagationClassification(::Missing) = missing
function MinimalPropagationClassification(
        l_frames::AbstractArray{<:AxisVector{T}}, 
        ts,
        xs, 
        periodic::Bool; 
        slope_min=DEFAULT_SLOPE_MIN,
        const_jitter=2abs(xs[2] - xs[1]),  # typically a few dx
        vel_jitter=DEFAULT_VEL_JITTER,
        min_dist_for_propagation=0.4(xs[end] - xs[begin])
    ) where T
    has_traveled_dist(rf::RunningFront) = abs(rf.previous_location - rf.starting_location) >= min_dist_for_propagation
    has_traveled_dist(::Nothing) = false
    return_fn = has_traveled_dist
    d1_ghost_op = make_ghost_op(T, xs, 1, periodic)
    dframe_cache = copy(l_frames[1])
    running_fronts = AxisVector{RF_possibilities(T)}(RF_possibilities(T)[nothing for _ in xs], collect(xs))
    continue_fronts!(running_fronts, ts[begin], dframe_cache, l_frames[1],
        return_fn, d1_ghost_op, slope_min, const_jitter, vel_jitter)
    for (t, frame) in zip(ts[begin+1:end], l_frames[begin+1:end])
        if continue_fronts!(running_fronts, t, dframe_cache, 
                frame, return_fn, d1_ghost_op, slope_min, const_jitter, vel_jitter) 
            return MinimalPropagationClassification(true)
        end
    end
    return MinimalPropagationClassification(false)
    # in-place continue_fronts onto empty running_fronts
    # loop over l_frames
    #   in-place continue_fronts on running_fronts
    #   return nothing from continue_fronts if condition not met
    #   if continue_fronts returns true, then return MinimalPropagationClassification(true)
    # return MinimalPropagationClassification(false)
end
function MinimalPropagationClassification(l_frames::AbstractArray{<:AbstractMatrix}, ts, multipop_xs, args...; kwargs...)
    MinimalPropagationClassification(
            [population(multipop_frame,1) for multipop_frame in l_frames], 
            ts, 
            multipop_xs, 
            args...; kwargs...)
end

function MinimalPropagationClassification(exec::Execution; kwargs...)
    l_frames::AbstractArray{<:AxisArray} = exec.solution.u
    ts = exec.solution.t
    multipop_xs = frame_xs(exec)
    periodic = exec.simulation.space isa AbstractPeriodicLattice
    MinimalPropagationClassification(l_frames, ts, multipop_xs, periodic; kwargs...)
end

function MinimalPropagationClassification(sol::DiffEqBase.AbstractTimeseriesSolution; kwargs...)
    l_frames = sol.u
    ts = sol.t
    xs = axes_keys(first(l_frames))[1]
    MinimalPropagationClassification(l_frames, ts, xs, true; kwargs...)
end

