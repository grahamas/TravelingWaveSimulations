
const DEFAULT_SLOPE_MIN = 1e-4
const DEFAULT_MAX_JITTER = 10.

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

last_quartile_dxs(fin) = ((3*fin)÷4):fin

function WaveClassifications(measurements::SpatiotemporalWaveMeasurements; 
                                 velocity_threshold,
                                 n_traveling_frames_threshold)
    if length(measurements.maxes) < 4
        # not long enough to be classified
        return WaveClassifications(false, false, false, false, false, false, false)
    end
    traveling = is_traveling(measurements.velocities, velocity_threshold, n_traveling_frames_threshold)
    unidirectional_travel = all(measurements.velocities .>= 0) || all(measurements.velocities .<= 0)
    decaying = is_decaying(measurements.maxes)
    growing = is_growing(measurements.maxes)
    oscillating = is_oscillating(measurements.maxes)
    positive_slope = measurements.slopes[1] > 0

    WaveClassifications(
        traveling,
        traveling && all(measurements.velocities .> 0),
        decaying,
        growing, 
        oscillating,
        positive_slope,
        unidirectional_travel #Are there other assumptions?
    )
end

function is_traveling(velocities::Vector{<:AbstractFloat}, 
                                 velocity_threshold,
                                 n_traveling_frames_threshold)
    sum(abs.(velocities) .> velocity_threshold) > n_traveling_frames_threshold
end
function is_decaying(maxes::Vector{<:AbstractFloat})
    fin = length(maxes)
    all(diff(maxes[last_quartile_dxs(fin)]) .<= -1e-8)
end
function is_growing(maxes::Vector{<:AbstractFloat})
    fin = length(maxes)
    all(diff(maxes[last_quartile_dxs(fin)]) .>= 1e-8)
end
function is_oscillating(maxes::Vector{<:AbstractFloat})
    # Checks if, within the second half, at least a third
    # are increasing, and a third decreasing
    fin = length(maxes)
    second_half = fin÷2:fin
    abs_eps = 1e-10
    just_increasing = maxes[second_half] .>= abs_eps
    just_decreasing = maxes[second_half] .<= -abs_eps
    min_num = length(second_half) ÷ 3
    sum(just_increasing) > min_num && sum(just_decreasing) > min_num
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
function ExecutionClassifications(wavefronts::WS, 
                                 ts::TS,
                                 xs::XS,
                                 final_frame::AbstractArray{T};
                                 max_resting=5e-2,
                                 origin_radius=20.0,
                                 velocity_threshold,
                                 n_traveling_frames_threshold) where {T,
                                    WS <: AbstractVector{<:AbstractVector{<:Wavefront}},
                                    TS <: AbstractVector{T},
                                    XS <: AbstractVector{T}
                                 }
    @assert origin_radius < xs[end]
    persistent_fronts = link_persistent_fronts(wavefronts, ts)
    all_fronts_velocities = get_velocities.(persistent_fronts)
    all_fronts_is_traveling = is_traveling.(all_fronts_velocities, velocity_threshold, n_traveling_frames_threshold)
    has_propagation = any(all_fronts_is_traveling)

    # pf_measurements = SpatiotemporalWaveMeasurements.(persistent_fronts)
    # pf_classifications = WaveClassifications.(pf_measurements; wave_kwargs...)

    # # TODO calculate first four bools with regard to propagation
    # has_propagation = any(map((cls) -> cls.traveling, pf_classifications))
    # has_oscillation = any(map((cls) -> cls.oscillating, pf_classifications))
    # farthest_propagation_is_decaying, farthest_propagation_is_oscillating = if has_propagation
    #     _, max_dx = findmax(map(x -> abs(x.distance), pf_measurements))
    #     (pf_classifications[max_dx].decaying, pf_classifications[max_dx].oscillating)
    # else
    #     (false, false)
    # end
    # persistently_active_near_origin = check_has_activity_near_origin(population(final_frame, 1),
    #                                                            #FIXME: E population assumption
    #                                                            max_resting,
    #                                                            xs,
    #                                                           origin_radius)
    # reaches_steady_state = if length(wavefronts[end]) == 0
    #     true 
    #     #FIXME: invalid if moving plateau... but that's fundamentally a steady state
    #     # could also check ts... if it ends before cutoff, then it had to have reached
    #     # steady state
    # else
    #     length(wavefronts[end]) == length(wavefronts[end-1]) && all(wavefronts[end] .== wavefronts[end-1]) 
    # end
    # ExecutionClassifications(
    #     has_propagation,
    #     has_oscillation,
    #     farthest_propagation_is_decaying,
    #     farthest_propagation_is_oscillating,
    #     persistently_active_near_origin,
    #     reaches_steady_state
    # )
    ExecutionClassifications(
        has_propagation,
        false,
        false,
        false,
        false,
        false
    )
end

function implies_origin_activity(wavefront, max_resting, radius)
    left, slope, right = points = left(wavefront), slope(wavefront), right(wavefront)
    active_point_near_origin = any(map(points) do point
        _val(point) > max_resting && (-radius <= _loc(point) <= radius)
    end)
    elevated_spanning_front = _loc(left) < -radius && _loc(right) > radius && _val(left) > max_resting && _val(right) > max_resting
    return active_point_near_origin || elevated_spanning_front
end


function check_has_activity_near_origin(wavefronts::AbstractVector{<:Wavefront}, max_resting, xs, radius)
    return any(implies_origin_activity.(wavefronts, Ref(max_resting), Ref(radius)))
end

function check_has_activity_near_origin(frame::AbstractVector{T}, max_resting, xs, radius) where {T <: Number}
    left_dx = findfirst(xs .>= -radius)
    right_dx = findfirst(xs .>= radius)
    return any(frame[left_dx:right_dx] .> max_resting)
end

function ExecutionClassifications(exec::Execution; slope_min=DEFAULT_SLOPE_MIN, periodic=true, kwargs...)
    @warn "should dispatch to solution"
    l_frames = exec.solution.u
    ts = timepoints(exec)
    xs = frame_xs(exec)
    l_frame_fronts = Vector{<:Wavefront{Float64}}[substantial_fronts(frame, periodic, slope_min) for frame in l_frames] #arr of arrs of fronts
    final_frame = l_frames[end]
    ExecutionClassifications(l_frame_fronts, ts, xs, final_frame; kwargs...)
end

# Handle case where already reduced to fronts
function ExecutionClassifications(exec::AugmentedExecution{T,W}; kwargs...) where {T, W <: AbstractArray{<:Wavefront}}
    @warn "should dispatch to solution"
    @assert exec.solution.t[end] > 0.0 # Needs solution to have final frame
    l_frame_fronts = exec.saved_values.saveval #arr of arrs of fronts
    xs = frame_xs(exec)
    ExecutionClassifications(l_frame_fronts, exec.saved_values.t, xs, exec.solution.u[end]; kwargs...)
end

function ExecutionClassifications(exec::ReducedExecution{T,W}; kwargs...) where {T, W <: AbstractArray{<:Wavefront}}
    @warn "should dispatch to solution"
    error("Needs final frame, but given ReducedExecution with no frames")
    @assert exec.solution.t[end] > 0.0 # Needs solution to have final frame
    l_frame_fronts = exec.saved_values.saveval #arr of arrs of fronts
    ExecutionClassifications(l_frame_fronts, exec.saved_values.t, exec.solution.u[end]; kwargs...)
end

using DiffEqBase
function ExecutionClassifications(sol::DiffEqBase.AbstractTimeseriesSolution; slope_min=DEFAULT_SLOPE_MIN, periodic=true, kwargs...)
    l_frames = sol.u
    @assert l_frames[1] isa AxisArray
    ts = sol.t
    xs = keys.(axes(l_frames[1]))
    @assert length(xs) == 2
    xs = xs[1]
    l_frame_fronts = Vector{<:Wavefront{Float64}}[substantial_fronts(frame, periodic, slope_min) for frame in l_frames]
    final_frame = l_frames[end]
    ExecutionClassifications(l_frame_fronts, ts, xs, final_frame; kwargs...)
end



##################################
### Propagation Classification ###
##################################

struct RunningFront{T_VAL, VEL<:Union{T_VAL,Missing}, T_LOC}
    starting_location::T_LOC
    previous_location::T_LOC
    previous_velocity::VEL
    slope::T_VAL
end
const RF = RunningFront{Float64, Float64, Float64}
const RFM = RunningFront{Float64, Missing, Float64}
const RF_RFM_NOTHING = Union{RF, RFM, Nothing}
_front_loc(front::AxisArray) = axes_keys(front) |> only |> only
_front_slope(front::AxisArray) = only(front)
predict_location(rf::RunningFront{T,T}) where T = rf.previous_location + rf.previous_velocity
predict_location(rf::RunningFront{T,Missing}) where T = rf.previous_location
function start_running_front(front::AxisArray)
    # Must be axisarray with only the slope at its location
    loc = _front_loc(front)
    val = _front_slope(front)
    RunningFront(loc, loc, missing, val)
end
function extend_running_front(rf::RunningFront, front::AxisArray)
    loc = _front_loc(front)
    val = _front_slope(front)
    RunningFront(
        rf.starting_location,
        loc,
        rf.previous_location - loc,
        val
    )
end

function front_continues_rf(front, running_front, max_jitter)
    """Return true if front can continue running_front.
        Needs:
            1. slope has same sign
            2. projected location is within max_jitter
    """
    slope = _front_slope(front)
    rf_slope = running_front.slope
    if approx_sign(rf_slope) != 0 && approx_sign(slope) != approx_sign(rf_slope)
        return false
    end

    loc = _front_loc(front)
    predicted_loc = predict_location(running_front)
    if !(predicted_loc - max_jitter <= loc <= predicted_loc + max_jitter)
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
                                    max_jitter)
    # haven't come across RF2
    return (preceding_running_front, unresolved_front_idx, nothing)    
end
function reconcile_running_front!(running_fronts, 
                                  new_running_front::RunningFront,
                                  frame, unresolved_front_idx,
                                  preceding_running_front::Nothing,
                                  max_jitter)
    # there's no RF1
    return (new_running_front, unresolved_front_idx, nothing)
end
function reconcile_running_front!(running_fronts, 
                                    new_running_front::RunningFront,
                                    frame, unresolved_front_idx::Nothing,
                                    preceding_running_front::RunningFront,
                                    max_jitter)
    # There's no F in RF1 F RF2, so just set RF2 as the new preceding front.
    # (also note there can't be F RF1 RF2... RF1 is just orphaned)
    return (new_running_front, unresolved_front_idx, nothing)
end
function reconcile_running_front!(running_fronts, 
                                    rf2::RunningFront,
                                    frame, front_idx::Int,
                                    rf1::RunningFront,
                                    max_jitter)::Tuple{Union{RF,RFM},Union{Int,Nothing},Union{RF,RFM,Nothing}}
    # Now we have RF1 F RF2  (actually! could be F RF1 RF2, but works anyway)
    @assert rf1 != rf2
    f = frame[[front_idx]]
    f_continues_rf1 = front_continues_rf(f, rf1, max_jitter)
    f_continues_rf2 = front_continues_rf(f, rf2, max_jitter)
    if f_continues_rf1 && !f_continues_rf2
        running_fronts[front_idx] = extend_running_front(rf1, f)
        return (rf2, nothing, running_fronts[front_idx])
    elseif !f_continues_rf1 && f_continues_rf2
        return (rf2, front_idx, nothing)
    elseif f_continues_rf1 && f_continues_rf2
        #@warn "Ambiguous front continuation!"
        running_fronts[front_idx] = extend_running_front(rf1, f)
        return (rf2, nothing, running_fronts[front_idx])
    else
        running_fronts[front_idx] = start_running_front(f)
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
            max_jitter)
    # don't have F1, so just set unresolved_front_idx to be that next time
    return (preceding_running_front, front_idx, nothing)
end
function reconcile_front!(running_fronts, front_idx, 
        frame, unresolved_front_idx::Int, 
        preceding_running_front::Nothing, 
        max_jitter)
    # have F1 but no RF, so F1 is orphaned.
    running_fronts[unresolved_front_idx] = start_running_front(frame[[unresolved_front_idx]])
    return (preceding_running_front, front_idx, running_fronts[unresolved_front_idx])
end
function reconcile_front!(running_fronts, f2_idx, 
        frame, f1_idx::Int, 
        rf::RunningFront, 
        max_jitter)::Tuple{RF_RFM_NOTHING,Int,Union{RF,RFM}}
    @assert f1_idx != f2_idx
    f1 = frame[[f1_idx]]
    f2 = frame[[f2_idx]]
    f1_continues_rf = front_continues_rf(f1, rf, max_jitter)
    f2_continues_rf = front_continues_rf(f2, rf, max_jitter)
    if f1_continues_rf && !f2_continues_rf
        running_fronts[f1_idx] = extend_running_front(rf, f1)
        return (nothing, f2_idx, running_fronts[f1_idx])
    elseif !f1_continues_rf && f2_continues_rf
        running_fronts[f1_idx] = start_running_front(f1)
        return (rf, f2_idx, running_fronts[f1_idx])
    elseif f1_continues_rf && f2_continues_rf
        #@warn "Ambiguous front continuation!"
        running_fronts[f1_idx] = extend_running_front(rf, f1)
        return (nothing, f2_idx, running_fronts[f1_idx])
    else
        running_fronts[f1_idx] = start_running_front(f1)
        return (nothing, f2_idx, running_fronts[f1_idx])
    end
end



function continue_fronts!(running_fronts::AxisVector, 
        dframe::AxisVector,
        frame::AxisVector, 
        return_fn::Function, 
        d1_ghost_op, slope_min, 
        max_jitter)::Bool
    deriv!(dframe, frame, d1_ghost_op)
    locations = axes_keys(frame) |> only
    prev_slope = dframe[begin]
    preceding_running_front = nothing
    unresolved_front_idx = nothing
    for idx in LinearIndices(locations)[begin+1:end-1]
        this_running_front = running_fronts[idx]
        running_fronts[idx] = nothing
        preceding_running_front, unresolved_front_idx, new_rf = reconcile_running_front!(running_fronts, this_running_front, 
                 dframe, unresolved_front_idx, 
                 preceding_running_front, max_jitter)
        return_fn(new_rf) && return true
        prev_slope, slope, next_slope = dframe[idx-1:idx+1]
        if abs(slope) > slope_min &&
                (prev_slope < slope >= next_slope) || 
                (prev_slope > slope <= next_slope)
            preceding_running_front, unresolved_front_idx, new_rf = reconcile_front!(running_fronts, idx, 
                    dframe, unresolved_front_idx,
                    preceding_running_front, 
                    max_jitter
                )
            return_fn(new_rf) && return true
        end
    end
    if !isnothing(unresolved_front_idx)
        if !isnothing(preceding_running_front) &&
            front_continues_rf(dframe[[unresolved_front_idx]], preceding_running_front, max_jitter)
            running_fronts[unresolved_front_idx] = extend_running_front( preceding_running_front, dframe[[unresolved_front_idx]])
        else
            running_fronts[unresolved_front_idx] = start_running_front(dframe[[unresolved_front_idx]])
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
# continue_fronts!(running_fronts, frame, return_fn(::RunningFront)::Union{STHG,Nothing})::Union{STHG,Nothing}
# Not really a classification....

export MinimalPropagationClassification
struct MinimalPropagationClassification <: AbstractExecutionClassifications
    has_propagation::Bool
end
MinimalPropagationClassification(::Union{AugmentedExecution,ReducedExecution}; kwargs...) = error("MinimalPropagationClassification should not be pre-reduced")
function MinimalPropagationClassification(l_frames::AbstractArray{<:AxisVector{T}},
        xs::AbstractArray, periodic; 
        slope_min=DEFAULT_SLOPE_MIN,
        max_jitter=DEFAULT_MAX_JITTER,
        min_dist_for_propagation
    ) where T
    has_traveled_dist(rf::RunningFront) = abs(rf.previous_location - rf.starting_location) >= min_dist_for_propagation
    has_traveled_dist(::Nothing) = false
    return_fn = has_traveled_dist
    d1_ghost_op = make_ghost_op(T, xs, 1, periodic)
    dframe_cache = copy(l_frames[1])
    running_fronts = AxisVector{RF_RFM_NOTHING}(RF_RFM_NOTHING[nothing for _ in xs], collect(xs))
    continue_fronts!(running_fronts, dframe_cache, l_frames[1],
        return_fn, d1_ghost_op, slope_min, max_jitter)
    for frame in l_frames[begin+1:end]
        if continue_fronts!(running_fronts, dframe_cache, 
                frame, return_fn, d1_ghost_op, slope_min, max_jitter) 
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
function MinimalPropagationClassification(l_frames::AbstractArray{<:AbstractMatrix}, multipop_xs, args...; kwargs...)
    MinimalPropagationClassification([population(multipop_frame,1) for multipop_frame in l_frames], multipop_xs, args...; kwargs...)
end

function MinimalPropagationClassification(exec::Execution; kwargs...)
    l_frames::AbstractArray{<:AxisArray} = exec.solution.u
    multipop_xs = frame_xs(exec)
    periodic = exec.simulation.space isa AbstractPeriodicLattice
    MinimalPropagationClassification(l_frames, multipop_xs, periodic; kwargs...)
end

function MinimalPropagationClassification(sol::DiffEqBase.AbstractTimeseriesSolution; kwargs...)
    l_frames = sol.u
    xs = axes_keys(first(l_frames))[1]
    MinimalPropagationClassification(l_frames, xs, true; kwargs...)
end