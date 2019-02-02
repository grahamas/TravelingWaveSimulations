abstract type AbstractNonlinearity{T} <: AbstractParameter{T} end

function update(calc_arr::AA, new_arr::AbstractArray{L,1}) where {T, L <: AbstractNonlinearity{T}, CL <: CalculatedType{L}, AA<:AbstractArray{CL,1}}
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

struct SigmoidNonlinearity{T} <: AbstractNonlinearity{T}
    a::T
    θ::T
end

SigmoidNonlinearity{T}(; a::T=nothing, θ::T=nothing) where T = SigmoidNonlinearity{T}(a,θ)

mutable struct CalculatedSigmoidNonlinearity{T} <: CalculatedType{SigmoidNonlinearity{T}}
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

CalculatedTypes.get_value(csn::CalculatedSigmoidNonlinearity{T}) where T = csn

############

### Sech2 ###

function sech2_fn(x, a, θ)
    return @. 1 - tanh(a * (x - θ))^2
end

struct Sech2Nonlinearity{T} <: AbstractNonlinearity{T}
    a::T
    θ::T
end

Sech2Nonlinearity{T}(; a::T=nothing, θ::T=nothing) where T = Sech2Nonlinearity{T}(a,θ) 

mutable struct CalculatedSech2Nonlinearity{T} <: CalculatedType{Sech2Nonlinearity{T}}
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

CalculatedTypes.get_value(csn::CalculatedSech2Nonlinearity{T}) where T = csn

############

### Gaussian ###

function gaussian_fn(x, sd, θ)
    return @. exp(-((x - θ) / sd)^2 ) - exp(-(-θ / sd))
end

struct GaussianNonlinearity{T} <: AbstractNonlinearity{T}
    sd::T
    θ::T
end

GaussianNonlinearity{T}(; sd::T=nothing, θ::T=nothing) where T = GaussianNonlinearity{T}(sd,θ) 

mutable struct CalculatedGaussianNonlinearity{T} <: CalculatedType{GaussianNonlinearity{T}}
    nonlinearity::GaussianNonlinearity{T}
    sd::T
    θ::T # This is nonsense, but fits with Calculated pattern
end

function Calculated(gaussian::GaussianNonlinearity{T}) where T
    CalculatedGaussianNonlinearity{T}(gaussian, gaussian.sd, gaussian.θ)
end

function nonlinearity!(output::AT, cgn::CalculatedGaussianNonlinearity{T}) where {T, AT<:AbstractArray{T}}
    output .= gaussian_fn.(output,cgn.sd,cgn.θ)
end
function nonlinearity(cgn::CalculatedGaussianNonlinearity{T}, input_arr::AT) where {T, AT<:Array{T}}
    ret_arr = copy(input_arr)
    nonlinearity!(ret_arr, cgn)
    return ret_arr
end

CalculatedTypes.get_value(cgn::CalculatedGaussianNonlinearity{T}) where T = cgn