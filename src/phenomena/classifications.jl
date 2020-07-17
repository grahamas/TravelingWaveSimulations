
### Single wave measurements across time ###

struct SpatiotemporalWaveMeasurements{T}
    velocities::Vector{T}
    slopes::Vector{T}
    maxes::Vector{T}
    distance::T
    duration::T
end
function SpatiotemporalWaveMeasurements(pf::Persistent{W,T}) where {W,T}
    SpatiotemporalWaveMeasurements{T}(
        get_velocities(pf),
        [wave.slope.val for wave in pf.waveforms],
        [max(wave).val for wave in pf.waveforms],
        pf.waveforms[end].slope.loc - pf.waveforms[begin].slope.loc,
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
                              kwargs...)
    if length(measurements.maxes) < 4
        # not long enough to be classified
        return WaveClassifications(false, false, false, false, false, false, false)
    end
    traveling = is_traveling(measurements.velocities; kwargs...)
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

function is_traveling(velocities::Vector{<:AbstractFloat};
                       velocity_threshold=1e-8, n_traveling_frames_threshold=5, kwargs...)
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
struct ExecutionClassifications
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
                                 kwargs...) where {T,
                                    WS <: AbstractVector{<:AbstractVector{<:Wavefront}},
                                    TS <: AbstractVector{T},
                                    XS <: AbstractVector{T}
                                 }
    @assert origin_radius < xs[end]
    persistent_fronts = link_persistent_fronts(wavefronts, ts)
    pf_measurements = SpatiotemporalWaveMeasurements.(persistent_fronts)
    pf_classifications = pf_measurements .|> (x) -> WaveClassifications(x; kwargs...)

    # TODO calculate first four bools with regard to propagation
    has_propagation = any(map((cls) -> cls.traveling, pf_classifications))
    has_oscillation = any(map((cls) -> cls.oscillating, pf_classifications))
    farthest_propagation_is_decaying, farthest_propagation_is_oscillating = if has_propagation
        _, max_dx = findmax(map(x -> abs(x.distance), pf_measurements))
        (pf_classifications[max_dx].decaying, pf_classifications[max_dx].oscillating)
    else
        (false, false)
    end
    persistently_active_near_origin = check_has_activity_near_origin(population(final_frame, 1),
                                                               #FIXME: E population assumption
                                                               max_resting,
                                                               xs,
                                                              origin_radius)
    reaches_steady_state = if length(wavefronts[end]) == 0
        true 
        #FIXME: invalid if moving plateau... but that's fundamentally a steady state
        # could also check ts... if it ends before cutoff, then it had to have reached
        # steady state
    else
        length(wavefronts[end]) == length(wavefronts[end-1]) && all(wavefronts[end] .== wavefronts[end-1]) 
    end
    ExecutionClassifications(
        has_propagation,
        has_oscillation,
        farthest_propagation_is_decaying,
        farthest_propagation_is_oscillating,
        persistently_active_near_origin,
        reaches_steady_state
    )


end

function implies_origin_activity(wavefront, max_resting, radius)
    left, slope, right = points = wavefront.left, wavefront.slope, wavefront.right
    active_point_near_origin = any(map(points) do point
        point.val > max_resting && (-radius <= point.loc <= radius)
    end)
    elevated_spanning_front = left.loc < -radius && right.loc > radius && left.val > max_resting && right.val > max_resting
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

function ExecutionClassifications(exec::Execution; kwargs...)
    @warn "should dispatch to solution"
    l_frames = exec.solution.u
    ts = timepoints(exec)
    xs = frame_xs(exec)
    l_frame_fronts = substantial_fronts.(l_frames, Ref(xs)) #arr of arrs of fronts
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
function ExecutionClassifications(sol::DiffEqBase.AbstractTimeseriesSolution; kwargs...)
    l_frames = sol.u
    @assert l_frames[1] isa AxisArray
    ts = sol.t
    xs = keys.(axes(l_frames[1]))
    @assert length(xs) == 2
    xs = xs[1]
    l_frame_fronts = substantial_fronts.(l_frames, Ref(xs)) #arr of arrs of fronts
    final_frame = l_frames[end]
    ExecutionClassifications(l_frame_fronts, ts, xs, final_frame; kwargs...)
end
