
### Single wave measurements across time ###

struct SpatiotemporalWaveMeasurements{T}
    velocities::Vector{T}
    slopes::Vector{T}
    maxes::Vector{T}
end
function SpatiotemporalWaveMeasurements(pf::Persistent{W,T}) where {W,T}
    SpatiotemporalWaveMeasurements{T}(
        get_velocities(pf),
        [slope(wave) for wave in pf.waveforms],
        [max(wave) for wave in pf.waveforms]
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

function is_traveling(velocities::Vector{T<:AbstractFloat})
    sum(abs.(velocities) .> 1e-8)
end
function is_decaying(maxes::Vector{T<:AbstractFloat})
    fin = length(maxes)
    all(diff(maxes[3*fin÷4:fin]) .<= 1e-10)
end
function is_growing(maxes::Vector{T<:AbstractFloat})
    fin = length(maxes)
    all(diff(maxes[3*fin÷4:fin]) .>= -1e-10)
end
function is_oscillating(maxes::Vector{T<:AbstractFloat})
    fin = length(maxes)
    second_half = fin÷2:fin
    abs_eps = 1e-10
    just_increasing = maxes[second_half] .>= abs_eps
    just_decreasing = maxes[second_half] .<= -abs_eps
    min_num = length(second_half) ÷ 3
    sum(just_increasing) > min_num && sum(just_decreasing) > min_num
end

### Whole execution classification ###


function get_execution_properties(exec::Execution; params...)
    l_frames = exec.solution.u
    ts = timepoints(exec)
    xs = [x[1] for x in space(exec).arr]
    l_frame_fronts = TravelingWaveSimulations.substantial_fronts.(l_frames, Ref(xs))
    get_execution_properties(l_frame_fronts, ts; params...)
end

function get_execution_properties(exec::Union{AugmentedExecution,ReducedExecution{<:Wavefront}}; params...)
    get_execution_properties(exec.saved_values.saveval, exec.saved_values.t; params...)
end

function get_execution_properties(nt::NamedTuple; params...)
    get_execution_properties(nt.wavefronts, nt.wavefronts_t; params...)
end

struct ExecutionProperties
    epileptic::Bool
    traveling_solitary::Bool
    decaying::Bool
    velocity::Union{Float64,Missing}
    velocity_error::Union{ Float64,Missing}
end
ExecutionProperties(; epileptic, traveling_solitary, decaying, velocity, velocity_error) = ExecutionProperties(epileptic, traveling_solitary, decaying, velocity, velocity_error)

