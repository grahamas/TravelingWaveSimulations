
abstract type AbstractConnectivity{T,D} <: AbstractParameter{T} end

function update(calc_arr::AA, new_arr::AbstractArray{C,2}, space::AbstractSpace{T}) where {T, C <: AbstractConnectivity{T}, CC <: CalculatedType{C}, AA<:AbstractArray{CC,2}}
    [calc_arr[i].connectivity != new_arr[i] ? Calculated(new_arr[i], space) : calc_arr[i] for i in CartesianIndices(calc_arr)]
end

#region ShollConnectivity

@calculated_type(struct ShollConnectivity{T} <: AbstractConnectivity{T,1}
    amplitude::T
    spread::T
end, function calculate(calculated_space::CalculatedType{<:AbstractSpace{T,1}}) where T
    directed_weights(source, calculated_space)
end
)

@calculated_type(struct GaussianConnectivity{T} <: AbstractConnectivity{T,2}
    amplitude::T
    spread::Tuple{T,T}
end, function calculate(calculated_space::CalculatedType{<:AbstractSpace{T,2}}) where T
    distances = get_distances(calculated_space)
    step_size = step(calculated_space)
    exponential_decay_gaussian.(distances, amplitude, Ref(spread), Ref(step_size))
end
)

###################

function exponential_decay_abs(distance::T, amplitude::T, spread::T, step_size::T) where {T <: Number}
    amplitude * step_size * exp(
        -abs(distance / spread)
    ) / (2 * spread)
end

function exponential_decay_gaussian(x::Tup, amplitude::T, spread::Tup, step_size::Tup) where {T, Tup<:Tuple{Vararg{T}}}
    amplitude * prod(step_size) * exp(
        -sum( (x ./ spread) .^ 2)
    ) / (2 * prod(spread))
end
# function exponential_decay_gaussian((x1, x2)::Tuple{Tup,Tup}, amplitude::T, spread::Tup, step_size::Tup) where {T, Tup<:Tuple{Vararg{T}}}
#     exponential_decay_gaussian(x1 .- x2, amplitude, spread, step_size)
# end


# * Sholl connectivity

function directed_weights(connectivity::ShollConnectivity{T}, locations::CalculatedType{<:AbstractSpace{T,1}}) where {T}
    A = connectivity.amplitude
    σ = connectivity.spread
    distances = get_distances(locations)
    step_size = step(locations)
    return exponential_decay_abs.(distances, A, σ, step_size)
    # In comparing to Neuman, there is no division by 2 on the edges
    # but the edges are very small, so there's no difference
end

function directed_weights(connectivity::GaussianConnectivity{T}, locations::CalculatedType{<:AbstractSpace{T,2}}) where {T}
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
