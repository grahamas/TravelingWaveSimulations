
abstract type AbstractAction{T} end
abstract type AbstractInteraction{T} <: AbstractAction{T} end
abstract type AbstractSpaceAction{T,N} <: AbstractAction{T} end
abstract type AbstractSpaceInteraction{T,N} <: AbstractSpaceAction{T,N} end
abstract type AbstractPopulationP{NP,AP} end

struct NullAction{T} <: AbstractAction{T} end
(na::NullAction)(args...) = nothing

struct CompositeAction{T,P} <: AbstractAction{T}
    arr::Vector{P}
end
CompositeAction(arr::Vector{A}) where {T, A<:AbstractAction{T}} = CompositeAction{T,A}(arr)
function (cp::CompositeParameter)(args...)
    actions = [param(args...) for param in cp.arr]
    types = Union{typeof.(actions)...}
    actions = Array{types}(actions)
    CompositeAction(actions)
end
function (ca::CompositeAction)(args...)
    for action in ca.arr
        action(args...)
    end
end

