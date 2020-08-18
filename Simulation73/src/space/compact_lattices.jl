
abstract type AbstractCompactLattice{T,N_ARR,N_CDT} <: AbstractLattice{T,N_ARR,N_CDT} end

@doc """
A lattice of points, each an `N_ARR` coordinate, evenly distributed.

Can be constructed using default AbstractLattice constructor
"""
struct CompactLattice{T,N_ARR} <: AbstractCompactLattice{T,N_ARR,N_ARR}
    arr::Array{NTuple{N_ARR,T},N_ARR}
end
difference(lattice::AbstractCompactLattice, edge) = abs_difference(edge)

# TODO: Should be completely uniform in periodic case?
export simpson_weights
function simpson_weights(lattice::CompactLattice{T,1}) where T
    @assert all(size(lattice) .% 2 .== 1)
    w = ones(T, size(lattice)...)
    w[1] = 0.5
    w[end] = 0.5
    return w
end
function simpson_weights(lattice::CompactLattice{T,2}) where T
    # http://mathfaculty.fullerton.edu/mathews/n2003/SimpsonsRule2DMod.html
    @assert all(size(lattice) .% 2 .== 1)
    w = ones(T, size(lattice)...)
    w[2:2:end-1,:] .*= 4.0
    w[3:2:end-2,:] .*= 2.0
    w[:,2:2:end-1] .*= 4.0
    w[:,3:2:end-2] .*= 2.0
    return w
end
Base.getindex(lat::CompactLattice, dx) = getindex(lat.arr, dx)

# Use default lattice plotting
