
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
function WaveClassifications(pf::Persistent)
    measurements = SpatiotemporalWaveMeasurements(pf)
    WaveClassifications(measurements)
end
function WaveClassifications(measurements::SpatiotemporalWaveMeasurements)
    traveling = is_traveling(measurements.velocities)
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

function is_traveling(velocities::Vector{<:AbstractFloat}, min_num_traveling_frames=5)
    sum(abs.(velocities) .> 1e-8) > min_num_traveling_frames
end
function is_decaying(maxes::Vector{<:AbstractFloat})
    fin = length(maxes)
    all(diff(maxes[3*fin÷4:fin]) .<= 1e-10)
end
function is_growing(maxes::Vector{<:AbstractFloat})
    fin = length(maxes)
    all(diff(maxes[3*fin÷4:fin]) .>= -1e-10)
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
export ExecutionClassifications
struct ExecutionClassifications
    has_propagation::Bool
    has_oscillation::Bool
    farthest_propagation_is_decaying::Bool
    farthest_propagation_is_oscillating::Bool
    persistently_active_near_origin::Bool
    reaches_steady_state::Bool
end

function ExecutionClassifications(wavefronts::WS, 
                                 ts::TS;
                                 max_resting=5e-2,
                                 left_length=15.0) where {T,
                                    WS <: AbstractVector{<:AbstractVector{<:Wavefront}},
                                    TS <: AbstractVector{T}
                                 }
    persistent_fronts = link_persistent_fronts(wavefronts, ts)
    pf_measurements = SpatiotemporalWaveMeasurements.(persistent_fronts)
    pf_classifications = WaveClassifications.(pf_measurements)

    # TODO calculate first four bools with regard to propagation
    has_propagation = any(map((cls) -> cls.traveling, pf_classifications))
    has_oscillation = any(map((cls) -> cls.oscillating, pf_classifications))
    farthest_propagation_is_decaying, farthest_propagation_is_oscillating = if has_propagation
        _, max_dx = findmax(map(x -> abs(x.distance), pf_measurements))
        (pf_classifications[max_dx].decaying, pf_classifications[max_dx].oscillating)
    else
        (false, false)
    end
    @show length.(wavefronts[end])
    persistently_active_near_origin = check_has_activity_near_left(wavefronts[end],
                                                               max_resting,
                                                               left_length)
    reaches_steady_state = all(wavefronts[end] .== wavefronts[end-1]) 
    ExecutionClassifications(
        has_propagation,
        has_oscillation,
        farthest_propagation_is_decaying,
        farthest_propagation_is_oscillating,
        persistently_active_near_origin,
        reaches_steady_state
    )


end

function check_has_activity_near_left(wavefronts, max_resting, left_length)
   long_enough_activation = false
   elevated_through = 0.0
   cur_wavefront_idx = 1 # FIXME should be firstidx
   while !long_enough_activation && (elevated_through < left_length)
       this_wavefront = wavefronts[cur_wavefront_idx]
       if this_wavefront.left.val <= max_resting
           break
       else
           elevated_through = this_wavefront.slope.loc
           if elevated_through >= left_length
               long_enough_activation = true
           else
               cur_wavefront_idx += 1
           end
       end
   end
   return long_enough_activation
end

function ExecutionClassifications(exec::Execution)
    l_frames = exec.solution.u
    ts = timepoints(exec)
    xs = [x[1] for x in space(exec).arr]
    l_frame_fronts = substantial_fronts.(l_frames, Ref(xs)) #arr of arrs of fronts
    ExecutionClassifications(l_frame_fronts, ts)
end

# Handle case where already reduced to fronts
function ExecutionClassifications(exec::Union{AugmentedExecution{T,W},ReducedExecution{T,W}}) where {T, W <: AbstractArray{<:Wavefront}}
    l_frame_fronts = exec.saved_values.saveval #arr of arrs of fronts
    @show length.(l_frame_fronts)
    ExecutionClassifications(l_frame_fronts, exec.saved_values.t)
end

function ExecutionClassifications(nt::NamedTuple; params...)
    ExecutionClassifications(nt.wavefronts, nt.wavefronts_t; params...)
end

