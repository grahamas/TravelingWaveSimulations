
abstract type AbstractAugmentedLattice{T,N_ARR,N_CDT,L} <: AbstractLattice{T,N_ARR,N_CDT} end
abstract type AbstractEmbeddedLattice{T,N_ARR,N_CDT,L} <: AbstractAugmentedLattice{T,N_ARR,N_CDT,L} end

function coordinates(lattice::AbstractEmbeddedLattice)
    lattice.coordinates
end
function size(lattice::AbstractEmbeddedLattice)
    size(lattice.lattice)
end

function difference(aug_lattice::AbstractEmbeddedLattice{T,N_ARR,N_CDT,L},
                    edge::Tuple{PT,PT}) where {T,N_ARR,N_CDT,L_N_CDT,
                                               L<:AbstractLattice{T,N_ARR,L_N_CDT},
                                               PT<:NTuple{N_CDT,T}
                                               }
    edge_first_dims = (edge[1][1:L_N_CDT], edge[2][1:L_N_CDT])
    edge_trailing_dims = (edge[1][L_N_CDT+1:end], edge[2][L_N_CDT+1:end])
    return (difference(aug_lattice.lattice, edge_first_dims)...,
        difference(aug_lattice.embedded_lattice, edge_trailing_dims)...)
end

function Base.step(aug_lattice::AbstractEmbeddedLattice)
    (step(aug_lattice.lattice)..., step(aug_lattice.embedded_lattice)...)
end

function unembed_values(lattice::AbstractEmbeddedLattice{T,N_ARR,N_CDT}, values::AbstractArray{T,N_ARR}) where {T,N_ARR,N_CDT}
    inner_coords = [coord[N_ARR+1:end] for coord in coordinates(lattice)]
    return [values[findall(map((x) -> all(isapprox.(embedded_coord, x)), inner_coords))] for embedded_coord in coordinates(lattice.embedded_lattice)]
end

linear_next(num::Int) = num + 1
#@recipe function f(lattice::AbstractEmbeddedLattice, values; layout=nothing, subplot=nothing)
#    @series begin
#        subplot := subplot
#        (lattice.lattice, values)
#    end
#    @series begin
#        subplot := linear_next(subplot)
#        (lattice.embedded_lattice, unembed_values(lattice, values))
#    end
#end


struct RandomlyEmbeddedLattice{T,N_ARR,N_CDT,L<:AbstractLattice{T,N_ARR},E<:AbstractSpace{T}} <: AbstractEmbeddedLattice{T,N_ARR,N_CDT,L}
    lattice::L
    embedded_lattice::E
    coordinates::Array{NTuple{N_CDT,T},N_ARR}
end
function RandomlyEmbeddedLattice(; lattice::L, embedded_lattice::E) where {T,N_ARR,L<:AbstractLattice{T,N_ARR},E<:AbstractSpace{T}}
    embedded_coordinates = embed_randomly(lattice, embedded_lattice)
    RandomlyEmbeddedLattice(lattice, embedded_lattice, embedded_coordinates)
end
function embed_randomly(lattice, embedded_lattice)
    [(lattice_coord..., sample(embedded_lattice)...) for lattice_coord in coordinates(lattice)]
end
function sample(lattice::AbstractLattice)
    (rand(length(extent(lattice))...) .* extent(lattice)) .- (extent(lattice) ./ 2)
end
