module WCMNonlinearity
# * Types

using Parameters
using CalculatedParameters
import CalculatedParameters: Calculated, update

abstract type Nonlinearity{T} <: Parameter{T} end

function update(calc_arr::AA, new_arr::AbstractArray{L,1}) where {T, L <: Nonlinearity{T}, CL <: CalculatedParam{L}, AA<:AbstractArray{CL,1}}
    [calc_arr[i].nonlinearity != new_arr[i] ? Calculated(new_arr[i]) : calc_arr[i] for i in CartesianIndices(calc_arr)]
end

### Sigmoid ###

"""
The sigmoid function is defined

```math
\\begin{align}
\\mathcal{S}(x) = \\frac{1}{1 + \\exp(-a(x - θ))}
\\end{align}
```
where ``a`` describes the slope's steepness and ``θ`` describes translation of the slope's center away from zero.

This is "simple" because in practice we use the rectified sigmoid.
"""
function simple_sigmoid_fn(x, a, theta)
    return @. (1.0 / (1 + exp(-a * (x - theta))))
end

"""
A rectified version of `simple_sigmoid_fn`.

In practice, we use rectified sigmoid functions because firing rates cannot be negative.

TODO: Rename to rectified_sigmoid_fn.
"""
function rectified_sigmoid_fn(x, a, theta)
    return max.(0, simple_sigmoid_fn(x, a, theta) .- simple_sigmoid_fn(0, a, theta))
end

function sigmoid_diff_fn(x, a, θ, width)
    return max.(0,sigmoid_fn(x, a, θ) - sigmoid_fn(x, a, θ + width))
end

function neg_domain_sigmoid_diff_fn(input, a, θ, width)
    return max.(0,simple_sigmoid_fn(input, a, θ) - simple_sigmoid_fn(input, a, θ + width))
end

@with_kw struct SigmoidNonlinearity{T} <: Nonlinearity{T}
    a::T
    θ::T
end

mutable struct CalculatedSigmoidNonlinearity{T} <: CalculatedParam{SigmoidNonlinearity{T}}
    nonlinearity::SigmoidNonlinearity{T}
    a::T
    θ::T # This is nonsense, but fits with Calculated pattern
end

function Calculated(sigmoid::SigmoidNonlinearity{T}) where T
    CalculatedSigmoidNonlinearity{T}(sigmoid, sigmoid.a, sigmoid.θ)
end

function nonlinearity!(output::AT, csn::CalculatedSigmoidNonlinearity{T}) where {T, AT<:AbstractArray{T}}
    output .= rectified_sigmoid_fn.(output,csn.a,csn.θ)
end
function nonlinearity(csn::CalculatedSigmoidNonlinearity{T}, input_arr::AT) where {T, AT<:Array{T}}
    ret_arr = copy(input_arr)
    nonlinearity!(ret_arr, csn)
    return ret_arr
end

CalculatedParameters.get_value(csn::CalculatedSigmoidNonlinearity{T}) where T = csn

############

### Sech2 ###

function sech2_fn(x, a, θ)
    return @. 1 - tanh(a * (x - θ))^2
end

@with_kw struct Sech2Nonlinearity{T} <: Nonlinearity{T}
    a::T
    θ::T
end

mutable struct CalculatedSech2Nonlinearity{T} <: CalculatedParam{Sech2Nonlinearity{T}}
    nonlinearity::Sech2Nonlinearity{T}
    a::T
    θ::T # This is nonsense, but fits with Calculated pattern
end

function Calculated(sech::Sech2Nonlinearity{T}) where T
    CalculatedSech2Nonlinearity{T}(sech, sech.a, sech.θ)
end

function nonlinearity!(output::AT, csn::CalculatedSech2Nonlinearity{T}) where {T, AT<:AbstractArray{T}}
    output .= sech2_fn.(output,csn.a,csn.θ)
end
function nonlinearity(csn::CalculatedSech2Nonlinearity{T}, input_arr::AT) where {T, AT<:Array{T}}
    ret_arr = copy(input_arr)
    nonlinearity!(ret_arr, csn)
    return ret_arr
end

CalculatedParameters.get_value(csn::CalculatedSech2Nonlinearity{T}) where T = csn


export Nonlinearity, Calculated, update, nonlinearity!, nonlinearity
export SigmoidNonlinearity, CalculatedSigmoidNonlinearity
export Sech2Nonlinearity, CalculatedSech2Nonlinearity
# * end
end
