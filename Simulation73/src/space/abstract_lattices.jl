
abstract type AbstractLattice{T,N_ARR,N_CDT} <: AbstractSpace{T,N_ARR,N_CDT} end
(t::Type{<:AbstractLattice})(start::Tuple, stop::Tuple, n_points::Tuple) = t(discrete_lattice(start, stop, n_points))
(t::Type{<:AbstractLattice{T,1}})(start::Number, stop::Number, n_points::Int) where T = t((start,),(stop,),(n_points,))
(t::Type{<:AbstractLattice})(; extent, n_points) = t(.-extent ./ 2, extent ./ 2, n_points)

Base.CartesianIndices(lattice::AbstractLattice) = CartesianIndices(lattice.arr)
Base.step(space::AbstractLattice) = extent(space) ./ (size(space) .- 1)
Base.size(lattice::AbstractLattice) = size(lattice.arr)
Base.size(lattice::AbstractLattice, d::Int) = size(lattice.arr, d)
Base.zeros(lattice::AbstractLattice{T}) where T = zeros(T,size(lattice)...)
start(space::AbstractLattice) = space.arr[1]
stop(space::AbstractLattice) = space.arr[end]
extent(space::AbstractLattice) = stop(space) .- start(space)


"""
    discrete_segment(extent, n_points)

Return an object containing `n_points` equidistant coordinates of a segment of length `extent` centered at 0. If you want 0 to be an element of the segment, make sure `n_points` is odd.

# Example
```jldoctest
julia> seg = discrete_segment(0.0, 5.0, 7);

julia> length(seg) == 7
true

julia> seg[end] - seg[1] ≈ 5.0
true
```
"""
function discrete_segment(start::T, stop::T, n_points::Int) where {T <: Number}
    LinRange{T}(start, stop, n_points)
end
# function discrete_segment(extent::T, n_points::Int) where {T <: Number}
#     discrete_segment(-extent/2,extent/2,n_points)
# end
"""
    coordinates(extent, n_points)

Return an object containing `n_points` equidistant coordinates along each dimension of a grid of length `extent` along each dimension, centered at (0,0,...,0).
"""
# function discrete_lattice(extent::NTuple{N,T}, n_points::NTuple{N,Int}) where {N,T}
#     Iterators.product(discrete_segment.(extent, n_points)...) |> collect
# end
function discrete_lattice(start::Tup, stop::Tup, n_points::IntTup) where {Nminusone,T,Tup<:Tuple{T,Vararg{T,Nminusone}},IT<:Int,IntTup<:Tuple{IT,Vararg{IT,Nminusone}}}
    Iterators.product(discrete_segment.(start, stop, n_points)...) |> collect
end
coordinates(lattice::AbstractLattice) = lattice.arr
coordinate_axes(lattice::AbstractLattice) = (discrete_segment.(start(lattice), stop(lattice), size(lattice))...,)

using Statistics
# @recipe function f(lattice::AbstractLattice, values::AbstractArray{<:AbstractArray,1})
#     standard_deviations = std.(values)
#     means = mean.(values)
#     up_stds = means .+ standard_deviations
#     down_stds = means .- standard_deviations
#     seriestype := :line
#     @series begin
#         linestyle := [:dot :dot :solid]
#         (lattice, hcat(up_stds, down_stds, means))
#     end
# end


# @recipe function f(lattice::AbstractLattice{T,1}, values; val_lim=nothing) where T
#     x := coordinate_axes(lattice)[1] |> collect
#     y := values
#     seriestype := :line
#     if val_lim != nothing
#         ylim := val_lim
#     end
#     ()
# end
# 
# @recipe function f(lattice::AbstractLattice{T,2}, values::Array{T,2}; val_lim=nothing) where T
#     (x, y) = coordinate_axes(lattice) .|> collect
#     seriestype := :heatmap
#     if val_lim != nothing
#         clim := val_lim
#         zlim := val_lim
#     end
#     aspect_ratio := :equal
#     (x,y,values)
# end

    
