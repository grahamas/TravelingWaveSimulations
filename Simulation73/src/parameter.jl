
export NoOpParameter, MaybeNoOp, MaybeComposite, MaybeNoOpComposite



struct NoOpParameter{T} <: AbstractParameter{T} end
NoOpParameter(param::AbstractParameter{T}) where T = NoOpParameter{T}()
(np::NoOpParameter{T})(space::AbstractSpace{T}) where {T} = NullAction{T}()
(ap::Type{<:AbstractParameter})(np::NoOpParameter) = np

const MaybeNoOp{T,P<:AbstractParameter{T}} = Union{P,NoOpParameter{T}}

struct CompositeParameter{T,P} <: AbstractParameter{T}
    arr::Vector{P}
end
function CompositeParameter(arr::Vector{P}) where {T, P<:AbstractParameter{T}}
    CompositeParameter{T,P}(arr)
end
function (ap::Type{<:AbstractParameter})(cp::CP) where {CP <: CompositeParameter}
    CompositeParameter(ap.(cp.arr))
end
const MaybeComposite{T,P<:AbstractParameter{T}} = Union{P,CompositeParameter{T,<:P}}
const MaybeNoOpComposite{T,P<:AbstractParameter{T}} = Union{P,CompositeParameter{T,<:MaybeNoOp{<:P}}}

struct NullifyParameter{T} <: AbstractParameter{T} end
NullifyParameter(param::AbstractParameter{T}) where T = NullifyParameter{T}()
(np::NullifyParameter{T})(space::AbstractSpace{T}) where {T} = NullAction{T}()
(ap::Type{<:AbstractParameter})(np::NullifyParameter) = np
