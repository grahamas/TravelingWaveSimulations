abstract type AbstractNonlinearity{T} <: AbstractAction{T} end

### Sigmoid ###
scalar_exp!(A) = A .= exp.(A)

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
    1.0 / (1 + exp(-a * (x - theta)))
end

function zeroed_sigmoid_fn(x, a, theta)
    simple_sigmoid_fn(x, a, theta) - simple_sigmoid_fn(0., a, theta)
end

"""
A rectified version of `simple_sigmoid_fn`.

In practice, we use rectified sigmoid functions because firing rates cannot be negative.

TODO: Rename to rectified_sigmoid_fn.
"""
function rectified_zeroed_sigmoid_fn(x, a, theta)
    max(0, zeroed_sigmoid_fn(x, a, theta))
end
struct ZeroedSigmoidNonlinearity{T} <: AbstractNonlinearity{T}
    a::T
    θ::T
    ZeroedSigmoidNonlinearity(a::T,θ::T) where T = new{T}(a,θ)
end
ZeroedSigmoidNonlinearity(; a, θ) = ZeroedSigmoidNonlinearity(a,θ)
(s::ZeroedSigmoidNonlinearity)(inplace::AbstractArray, ignored_source, ignored_t) = inplace .= rectified_zeroed_sigmoid_fn.(inplace, s.a, s.θ)

function rectified_unzeroed_sigmoid_fn(x, a, theta)
    max(0, simple_sigmoid_fn(x, a, theta))
end
near_zero_start(a,θ) = 0.0 <= rectified_unzeroed_sigmoid_fn(0.0,a,θ) < 0.05
struct SigmoidNonlinearity{T} <: AbstractNonlinearity{T}
    a::T
    θ::T
    function SigmoidNonlinearity(a::T,θ::T) where T
        if !near_zero_start(a, θ)
            return missing
        end
        new{T}(a,θ)
    end
end
SigmoidNonlinearity(; a, θ) = SigmoidNonlinearity(a,θ)
(s::SigmoidNonlinearity)(inplace::AbstractArray, ignored_source, ignored_t) = inplace .= rectified_unzeroed_sigmoid_fn.(inplace, s.a, s.θ)

struct UnrectifiedSigmoidNonlinearity{T} <: AbstractNonlinearity{T}
    a::T
    θ::T
end
UnrectifiedSigmoidNonlinearity(; a, θ) = UnrectifiedSigmoidNonlinearity(a,θ)
(s::UnrectifiedSigmoidNonlinearity)(inplace::AbstractArray, ignored_source, ignored_t) = inplace .= simple_sigmoid_fn.(inplace, s.a, s.θ)

############

### Sech2 ###

function sech2_fn(x, a, θ)
    @. 1 - tanh(a * (x - θ))^2
end
struct Sech2Nonlinearity{T} <: AbstractNonlinearity{T}
    a::T
    θ::T
    Sech2Nonlinearity(a::T,θ::T) where T = new{T}(a,θ)
end
Sech2Nonlinearity(; a, θ) = Sech2Nonlinearity(a,θ)
(sn::Sech2Nonlinearity)(output::AbstractArray, ignored_source, ignored_t) = output .= sech2_fn.(output,sn.a,sn.θ)

############

### Gaussian ###

function gaussian_fn(x, sd, θ)
    @. exp(-((x - θ) / sd)^2 ) - exp(-(-θ / sd)^2)
end
struct GaussianNonlinearity{T} <: AbstractNonlinearity{T}
    sd::T
    θ::T
    GaussianNonlinearity(sd::T,θ::T) where T = new{T}(sd,θ)
end
GaussianNonlinearity(; sd, θ) = GaussianNonlinearity(sd,θ)
(gaussian::GaussianNonlinearity)(output::AbstractArray, ignored_source, ignored_t) = output .= gaussian_fn.(output,gaussian.sd,gaussian.θ)

##############

function nonnegative_everywhere(fsig::SigmoidNonlinearity{T},bsig::SigmoidNonlinearity{T}) where T
    max_θ = max(fsig.θ, bsig.θ)
    min_a = min(fsig.a, fsig.a)
    max_test_val = max_θ + (10.0 / min_a)
    test_step = min_a / 10.0
    test_vals = 0.0:test_step:max_test_val |> collect
    dos_fn(fsig, bsig, test_vals)
    all(test_vals .>= 0)
end

### Difference of Sigmoids
struct DifferenceOfSigmoids{T} <: AbstractNonlinearity{T}
    firing_sigmoid::SigmoidNonlinearity{T}
    blocking_sigmoid::SigmoidNonlinearity{T}
    function DifferenceOfSigmoids(fsig::SigmoidNonlinearity{T},bsig::SigmoidNonlinearity{T}) where T
        if !nonnegative_everywhere(fsig, bsig)
            return missing
        end
        new{T}(fsig,bsig)
    end
end
DifferenceOfSigmoids(::Any, ::Missing) = missing
DifferenceOfSigmoids(::Missing, ::Any) = missing
DifferenceOfSigmoids(::Missing, ::Missing) = missing

DifferenceOfSigmoids(; firing_a, firing_θ, blocking_a, blocking_θ) = DifferenceOfSigmoids(SigmoidNonlinearity(; θ=firing_θ,a=firing_a), SigmoidNonlinearity(; θ=blocking_θ,a=blocking_a))
function dos_fn(up_sig::SigmoidNonlinearity, down_sig::SigmoidNonlinearity, output::AbstractArray)
    output .= rectified_unzeroed_sigmoid_fn.(output, up_sig.a, up_sig.θ) .- rectified_unzeroed_sigmoid_fn.(output, down_sig.a, down_sig.θ)
end
    
function (dos::DifferenceOfSigmoids)(output::AbstractArray, ignored_source, ignored_t)
    dos_fn(dos.firing_sigmoid, dos.blocking_sigmoid, output)
end
