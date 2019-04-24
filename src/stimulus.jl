
abstract type AbstractStimulus{T,D} <: AbstractParameter{T} end

@memoize Dict function make_mutator(stimulus_arr::AbstractArray{<:AbstractStimulus{T}}, space::AbstractSpace) where T
    stimulus_mutators = [make_stimulus(stim, space) for stim in stimulus_arr]
    function stimulus_mutator!(dA::AbstractArray{T,D}, A::AbstractArray{T,D}, t::T) where {T,D}
        for (i, stimulus!) in enumerate(stimulus_mutators)
            stimulus!(view_slice_last(dA, i), t)
        end
    end
end


function distance(x,y) where {T <: Real}
    sqrt(sum((x .- y).^2))
end


struct NoStimulus{T,N} <: AbstractStimulus{T,N} end
function make_stimulus(nostim::NoStimulus{T,N}, space::AbstractSpace{T,N}) where {T,N,AT<: AbstractArray{T,N}}
    (val,t) -> return
end

struct GaussianNoiseStimulus{T,N} <: AbstractStimulus{T,N}
    mean::T
    sd::T
end
function GaussianNoiseStimulus{T,N}(; SNR::T=0.0, mean::T=0.0) where {T,N}
    sd = sqrt(1/10 ^ (SNR / 10))
    GaussianNoiseStimulus{T,N}(mean, sd)
end
function gaussian_noise!(val::AT, mean::T, sd::T) where {T, AT<:AbstractArray{T}} # assumes signal power is 0db
    randn!(val)
    val .*= sd
    val .+= mean
end
function make_stimulus(wns::GaussianNoiseStimulus{T}, space::AbstractSpace{T}) where {T}
    (val,t) -> gaussian_noise!(val, wns.mean, wns.sd) # Not actually time dependent
end

abstract type TransientBumpStimulus{T,N} <: AbstractStimulus{T,N} end
function make_stimulus(bump::TBS, space::AbstractSpace{T,N}) where {T, N, TBS<:TransientBumpStimulus{T,N}}
    bump_frame = on_frame(bump, space)
    onset = bump.time_window[1]
    offset = bump.time_window[2]
    function stimulus!(val, t)
        if onset <= t < offset
            val .+= bump_frame
        end
    end
end

struct SharpBumpStimulus{T,N} <: TransientBumpStimulus{T,N}
    width::T
    strength::T
    time_window::Tuple{T,T}
    center::NTuple{N,T}
end

function SharpBumpStimulus{T,N}(; strength::T=nothing, width::T=nothing,
        duration=nothing, time_window=nothing, center=NTuple{N,T}(zero(T) for i in 1:N)) where {T,N}
    if time_window == nothing
        return SharpBumpStimulus{T,N}(width, strength, (zero(T), duration), center)
    else
        @assert duration == nothing
        return SharpBumpStimulus{T,N}(width, strength, time_window, center)
    end
end

function on_frame(sbs::SharpBumpStimulus{T,N}, space::AbstractSpace{T,N}) where {T,N}
    pop = one_pop(space)
    frame = zeros(T,size(pop)...)
    half_width = sbs.width / 2.0
    frame[distance.(pop, Ref(sbs.center)) .<= half_width] .= sbs.strength
    return frame
end


struct NoisyStimulus{T,N,STIM} <: AbstractStimulus{T,N}
    noise::GaussianNoiseStimulus{T,N}
    stimulus::STIM
end
function NoisyStimulus{T,N}(; stim_type::Type=SharpBumpStimulus{T,N}, SNR::T, mean::T=0.0, kwargs...) where {T,N}
    @show stim_type
    NoisyStimulus{T,N,stim_type}(
            GaussianNoiseStimulus{T,N}(SNR = SNR, mean=mean),
            stim_type(; kwargs...)
        )
end
@memoize Dict function make_stimulus(noisy_stimulus::NoisyStimulus{T}, space::AbstractSpace) where {T}
    noise_mutator = make_stimulus(noisy_stimulus.noise, space)
    stim_mutator = make_stimulus(noisy_stimulus.stimulus, space)
    (val, t) -> (noise_mutator(val,t); stim_mutator(val,t))
end
