abstract type AbstractNonlinearity{T} <: AbstractParameter{T} end

@memoize Dict function make_mutator(nonlinearity_arr::AbstractArray{<:AbstractNonlinearity{T}}) where T
    nonlinearity_mutators = [make_nonlinearity(nonl) for nonl in nonlinearity_arr]
    function nonlinearity_mutator!(dA::AbstractArray{T,D}, A::AbstractArray{T,D}, t::T) where {T,D}
        for (i, nonlinearity!) in enumerate(nonlinearity_mutators)
            nonlinearity!(view_slice_last(dA, i))
        end
    end
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
struct SigmoidNonlinearity{T} <: AbstractNonlinearity{T}
    a::T
    θ::T
end
SigmoidNonlinearity{T}(; a::T=nothing, θ::T=nothing) where T = SigmoidNonlinearity{T}(a,θ)

@memoize Dict function make_nonlinearity(sn::SigmoidNonlinearity{T}) where {T}
    (output::AbstractArray{T}) -> output .= rectified_sigmoid_fn.(output,sn.a,sn.θ)
end

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
@memoize Dict function make_nonlinearity(sn::Sech2Nonlinearity{T}) where {T}
    (output::AbstractArray{T}) -> output .= sech2_fn.(output,sn.a,sn.θ)
end

############

### Gaussian ###

function gaussian_fn(x, sd, θ)
    return @. exp(-((x - θ) / sd)^2 ) - exp(-(-θ / sd)^2)
end
struct GaussianNonlinearity{T} <: AbstractNonlinearity{T}
    sd::T
    θ::T
end
GaussianNonlinearity{T}(; sd::T=nothing, θ::T=nothing) where T = GaussianNonlinearity{T}(sd,θ)
@memoize Dict function make_nonlinearity(gn::GaussianNonlinearity{T}) where {T}
    (output::AbstractArray{T}) -> output .= gaussian_fn.(output,gn.sd,gn.θ)
end
