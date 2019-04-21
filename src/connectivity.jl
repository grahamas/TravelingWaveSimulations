
abstract type AbstractConnectivity{T,D} <: AbstractParameter{T} end

macro make_make_connectivity_mutator(num_dims)
    D = eval(num_dims)
    to_syms = [Symbol(:to,:_,i) for i in 1:D]
    from_syms = [Symbol(:from,:_,i) for i in 1:D]
    D_P = D + 1
    D_CONN_P = D + D + 2
    D_CONN = D + D
    tensor_prod_expr = @eval @macroexpand @tensor dA[$(to_syms...),i] = dA[$(to_syms...),i] + connectivity_tensor[$(to_syms...),$(from_syms...),i,j] * A[$(from_syms...),j]
    quote
        @memoize Dict function make_mutator(conn::AbstractArray{<:AbstractConnectivity{T,$D}}, space::AbstractSpace{T,$D}) where {T}
            connectivity_tensor::Array{T,$D_CONN_P} = tensor(conn, space)
            function connectivity!(dA::Array{T,$D_P}, A::Array{T,$D_P}, t::T) where T
                $tensor_prod_expr
            end
        end
    end |> esc
end

@make_make_connectivity_mutator 1
@make_make_connectivity_mutator 2

#
@with_kw struct ShollConnectivity{T} <: AbstractConnectivity{T,1}
    amplitude::T
    spread::T
end
function tensor(connectivity::ShollConnectivity{T}, space::AbstractSpace{T,1}) where T
    directed_weights(connectivity, space)
end
function tensor(arr::AbstractArray{<:AbstractConnectivity{T,N}}, space::Pops{P,T,N}) where {T,N,P}
    ret_tensor = Array{T,(N+N+2)}(undef, one_pop_size(space)..., one_pop_size(space)..., P, P)
    for dx in CartesianIndices(arr)
        view_slice_last(ret_tensor, dx) .= tensor(arr[dx], space)
    end
    return ret_tensor
end


@with_kw struct GaussianConnectivity{T} <: AbstractConnectivity{T,2}
    amplitude::T
    spread::Tuple{T,T}
end
function tensor(connectivity::GaussianConnectivity{T}, space::AbstractSpace{T,2}) where T
    directed_weights(connectivity, space)
end

###################

function exponential_decay_abs(distance::T, amplitude::T, spread::T, step_size::T) where {T <: Number}
    amplitude * step_size * exp(
        -abs(distance / spread)
    ) / (2 * spread)
end

function exponential_decay_gaussian(coord_distances::Tup, amplitude::T, spread::Tup, step_size::Tup) where {T,N, Tup<:NTuple{N,T}}
    amplitude * prod(step_size) * exp(
        -sum( (coord_distances ./ spread) .^ 2)
    ) / (2 * prod(spread))
end

# * Sholl connectivity

function directed_weights(connectivity::ShollConnectivity{T}, locations::AbstractSpace{T,1}) where {T}
    A = connectivity.amplitude
    σ = connectivity.spread
    distances = get_distances(locations)
    step_size = step(locations)
    return exponential_decay_abs.(distances, A, σ, step_size)
    # In comparing to Neuman, there is no division by 2 on the edges
    # but the edges are very small, so there's no difference
end

function directed_weights(connectivity::GaussianConnectivity{T}, locations::AbstractSpace{T,2}) where {T}
    distances = get_distances(locations)
    step_size = step(locations)
    return exponential_decay_gaussian.(distances, connectivity.amplitude, Ref(connectivity.spread), step_size)
end


#function calculate(connectivity::GaussianConnectivity, calc_dist_tensor::CalculatedDistanceTensor)

"""
We use an exponential connectivity function, inspired both by Sholl's
experimental work, and by certain theoretical considerations.

The interaction between two populations is characterized by this
function and its two parameters: the amplitude (weight) and the spread
(σ). The spatial step size is also a factor, but as a computational concern
rather than a fundamental one.
"""



# function meijer_matrix(connectivity::MeijerConnectivity, calc_dist_mx::CalculatedDistanceMatrix)
#     A = connectivity.amplitude
#     σ = connectivity.spread
#     dist_mx = calc_dist_mx.value
#     step_size = step(calc_dist_mx)
#     return meijer_matrix(A, σ, dist_mx, step_size)
#     # In comparing to Neuman, there is no division by 2 on the edges
#     # but the edges are very small, so there's no difference
# end
#


#endregion
