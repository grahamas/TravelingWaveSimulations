
abstract type AbstractPeriodicLattice{T,N_ARR,N_CDT} <: AbstractLattice{T,N_ARR,N_CDT} end

@doc """
A Lattice of points with `extent` describing the length along each dimension and `n_points` describing the number of points representing each dimension.
"""
struct PeriodicLattice{T,N_ARR} <: AbstractPeriodicLattice{T,N_ARR,N_ARR}
    arr::Array{NTuple{N_ARR,T},N_ARR}
end
function (t::Type{<:PeriodicLattice})(; extent, n_points)
    t(.-extent ./ 2, extent ./ 2 .- (extent ./ n_points), n_points)
end
difference(p_lattice::PeriodicLattice, edge) = abs_difference_periodic(edge, extent(p_lattice))

# @recipe function f(lattice::PeriodicLattice{T,1}, values; val_lim=(0.0, 1.0)) where T
#     θ = coordinate_axes(lattice)[1] .* (2π / extent(lattice)[1])
#     #y = values .* sin.(θ)
#     #x = values .* cos.(θ)
#     seriestype := :scatter
#     projection := :polar
#     min_val, max_val = val_lim
#     min_r = max(abs(min_val), abs(max_val))
#     inner_r = min_r * 2
#     r_values = values .+ inner_r
#     #lims := (0, inner_r + max_val)

#     aspect_ratio := :equal

#     (θ |> collect, r_values)
# end

# @recipe function f(lattice::PeriodicLattice{T,2}, values::Array{T,2}; val_lim=nothing) where T
#     (x, y) = coordinate_axes(lattice) .|> collect
#     seriestype := :heatmap
#     if val_lim != nothing
#
#     (x,y,values)
# end

# TODO: Should be completely uniform in periodic case?
export simpson_weights
function simpson_weights(lattice::PeriodicLattice{T,1}) where T
    @assert all(size(lattice) .% 2 .== 0)
    w = ones(T, size(lattice)...)
    return w
end
function simpson_weights(lattice::PeriodicLattice{T,2}) where T
    # http://mathfaculty.fullerton.edu/mathews/n2003/SimpsonsRule2DMod.html
    @assert all(size(lattice) .% 2 .== 0)
    w = ones(T, size(lattice)...)
    w[1:2:end-1,:] .*= 4.0
    w[2:2:end,:] .*= 2.0
    w[:,1:2:end-1] .*= 4.0
    w[:,2:2:end] .*= 2.0
    return w
end
