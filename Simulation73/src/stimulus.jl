
abstract type AbstractStimulusParameter{T} <: AbstractParameter{T} end
abstract type AbstractStimulusAction{T,D} <: AbstractSpaceAction{T,D} end

# Naturally map stimulus arrays
function (stim_params::AbstractArray{<:AbstractStimulusParameter{T}})(space::AbstractSpace{T,N}) where {T,N}
    return map(stim_params) do param
        param(space)
    end
end
# Stimulus actions are applied IN ORDER
function (stims::AbstractArray{<:AbstractStimulusAction})(args...)
    for stim in stims
        stim(args...)
    end
end

struct NoStimulusParameter{T} <: AbstractStimulusParameter{T} end
struct NoStimulusAction{T,N} <: AbstractStimulusAction{T,N} end
function (nostim::NoStimulusParameter{T})(space::AbstractSpace{T,N}) where {T,N}
    NoStimulusAction{T,N}()
end
(ns::NoStimulusAction)(args...) = nothing
