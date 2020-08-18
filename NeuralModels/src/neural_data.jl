
abstract type AbstractHomogeneousNeuralData{T,N} <: AbstractArray{T,N} end
const AbstractHeterogeneousNeuralData{T,N} = AbstractArray{T,N}
