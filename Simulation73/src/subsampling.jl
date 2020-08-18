
# Subsampling space can definitely be necessary
# Subsampling time may not be, maybe rather use DESolution methods.

# Supposing adaptive time, specify "saveat" -- simple.
# Optional, allow specifying "tstops" which forces non-interpolation.


import Base: to_indices, _maybetail, @_inline_meta, tail, getindex

abstract type AbstractSubsampler{D} end
function coordinate_indices(lattice::LATTICE, 
        subsampler_arr::AbstractArray{<:AbstractSubsampler}) where {LATTICE<:AbstractLattice}
     idxs_arr = map(subsampler_arr) do subsampler
        coordinate_indices(lattice, subsampler)
    end
    all_keep = ones(Bool, size(lattice)...)
    for idxs in idxs_arr
        keep = zeros(Bool, size(lattice)...)
        keep[idxs] .= true
        all_keep .&= keep
    end
    all_indices = CartesianIndices(lattice)
    return all_indices[all_keep]
end
function subsample(lattice::LATTICE, sub) where {LATTICE<:AbstractLattice}
    idxs = coordinate_indices(lattice, sub)
    LATTICE(lattice.arr[idxs])
end
subsample(space::AbstractLattice, ::Nothing) = space
function getindex(lattice::LATTICE, sub) where {T,D,LATTICE<:AbstractLattice{T,D}}
    idxs = coordinate_indices(lattice, sub)
    LATTICE(lattice.arr[idxs])
end

@with_kw struct IndexSubsampler{D} <: AbstractSubsampler{D}
	strides::NTuple{D,Int}
end
@with_kw struct RightCutFromValue{D,T} <: AbstractSubsampler{D}
	cut::NTuple{D,T}
end
export RightCutFromValue

@with_kw struct RightCutProportionFromValue{D,T} <: AbstractSubsampler{D}
    cut::NTuple{D,T}
    proportion::NTuple{D,T}
end
export RightCutProportionFromValue

"IndexInfo specifies the index of 0 and the stride (`Δ`) between indices."
@with_kw struct IndexInfo{D,N}
	Δ::NTuple{N,D}
	origin_idx::CartesianIndex{N}
end

"""
	StrideToEnd(stride, start=1)

StrideToEnd is a custom index type that acts like `start:stride:end`, to circumvent the fact that you can't put `start:stride:end` into a variable.
"""
struct StrideToEnd{STOP}
    stride::Int
	start::Int
    stop::STOP
    StrideToEnd{STOP}(stride::Int, start::Int=1, stop::STOP=nothing) where {STOP<:Union{Nothing,Int}} = new(stride, start,stop)
end
StrideToEnd(stride, start, stop::S=nothing) where S = StrideToEnd{S}(stride, start, stop)
to_indices(A, inds, I::Tuple{StrideToEnd{Int}, Vararg{Any}})=
(@_inline_meta; (I[1].start:I[1].stride:I[1].stop, to_indices(A, _maybetail(inds), tail(I))...))
to_indices(A, inds, I::Tuple{StrideToEnd{Nothing}, Vararg{Any}})=
	(@_inline_meta; (I[1].start:I[1].stride:inds[1][end], to_indices(A, _maybetail(inds), tail(I))...))
to_indices(A, inds, I::Tuple{NTuple{N,StrideToEnd}, Vararg{Any}}) where N = to_indices(A, inds, (I[1]..., _maybetail(I)...))
getindex(A::AbstractArray, S::StrideToEnd) = getindex(A, to_indices(A, (S,))...)

function coordinate_indices(lattice::AbstractLattice{T}, subsampler::RightCutFromValue{D,T}) where {T,D}
    axis_dxs = map(zip(subsampler.cut, coordinate_axes(lattice))) do (value, axis)
        findfirst(axis .>= value)
    end
    stride = 1
    cuts = Tuple(axis_dxs)
    return StrideToEnd.(stride, cuts)
end
function coordinate_indices(lattice::AbstractLattice{T}, subsampler::RightCutProportionFromValue{D,T}) where {T,D}
    start_stop = map(zip(subsampler.cut, subsampler.proportion, coordinate_axes(lattice))) do (value, proportion, axis)
        axis_begin = findfirst(axis .>= value)
        remaining_span = length(axis) - axis_begin
        axis_end = floor(Int, axis_begin + remaining_span * proportion)
        (axis_begin, axis_end)
    end
    axis_starts, axis_stops = zip(start_stop...)
    return StrideToEnd.(1, axis_starts, axis_stops)
end

function coordinate_indices(::Any, subsampler::IndexSubsampler{D}) where D
    return StrideToEnd.(subsampler.strides, 1)
end

function coordinate_indices(lattice, sub_dxs::AbstractArray{<:CartesianIndex}) 
    #@warn "subsampling one pop"
    #FIXME don't subsample one pop
    idxs = population(sub_dxs, 1)
end
