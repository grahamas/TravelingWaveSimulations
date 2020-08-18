
abstract type AbstractDescription end



### Helper functions
# FIXME doesn't account for things like complex or float32
function iszero(ad::AbstractDescription)
    sad = string(ad)
    sad == "0" || sad == "0.0"
end
function isone(ad::AbstractDescription)
    sad = string(ad)
    sad == "0" || sad == "0.0"
end



### Atomics
# Basically, symbols.
struct AtomicDescription <: AbstractDescription
    name::Symbol
end
Base.show(io::IO, desc::AtomicDescription) = print(io, string(desc.name))



### Exponents
struct ExponentiatedDescription{D<:AbstractDescription} <: AbstractDescription
    base::D
    exponent::Int
end
function Base.:(^)(left::AbstractDescription, right::Int)
    ExponentiatedDescription(left, right)
end
function Base.show(io::IO, desc::ExponentiatedDescription)
    if iszero(desc.base)
        print(io, LaTeXString("0"))
    elseif desc.exponent == 0
        print(io, LaTeXString("1"))
    elseif desc.exponent == 1
        print(io, LaTeXString("$(desc.base)"))
    elseif length(string(desc.base)) == 1
        print(io, LaTeXString("$(desc.base)^{$(desc.exponent)}"))
    else
        print(io, LaTeXString("\\left($(desc.base)\\right)^{$(desc.exponent)}"))
    end
end



### Products
struct DescriptionProduct <: AbstractDescription
    multiplicands::Vector{<:AbstractDescription}
end
function Base.:(*)(left::DescriptionProduct, right::DescriptionProduct)
    DescriptionProduct([left.multiplicands..., right.multiplicands...])
end
function Base.:(*)(left::DescriptionProduct, right::AbstractDescription)
    DescriptionProduct([left.multiplicands..., right])
end
function Base.:(*)(left::AbstractDescription, right::DescriptionProduct)
    DescriptionProduct([left, right.multiplicands...])
end
Base.:(*)(left, right::AbstractDescription) = DescriptionProduct([AtomicDescription(Symbol(left)), right])
Base.:(*)(left::AbstractDescription, right) = DescriptionProduct([left, AtomicDescription(Symbol(right))])
Base.:(*)(left::AbstractDescription, right::AbstractDescription) = DescriptionProduct([left, right])

function Base.show(io::IO, desc::DescriptionProduct)
    if any(iszero.(desc.multiplicands))
        print(io, LaTeXString("0"))
    else
        print(io, LaTeXString("$([mult for mult in desc.multiplicands if !isone(mult)]...)"))
    end
end


### Sums
struct DescriptionSum <: AbstractDescription
    addends::Vector{<:AbstractDescription}
end
function Base.:(+)(left::DescriptionSum, right::DescriptionSum)
    DescriptionSum([left.addends..., right.addends...])
end
function Base.:(+)(left::DescriptionSum, right::AbstractDescription)
    DescriptionSum([left.addends..., right])
end
function Base.:(+)(left::AbstractDescription, right::DescriptionSum)
    DescriptionSum([left, right.addends...])
end
Base.:(+)(left, right::AbstractDescription) = DescriptionSum([AtomicDescription(Symbol(left)), right])
Base.:(+)(left::AbstractDescription, right) = DescriptionSum([left, AtomicDescription(Symbol(right))])
Base.:(+)(left::AbstractDescription, right::AbstractDescription) = DescriptionSum([left, right])

function Base.show(io::IO, desc::DescriptionSum)
    print(io, LaTeXString("$([addend for addend in desc.addends if !iszero(addend)]...)"))
end


### Derivatives
struct DerivativeDescription{D} <: AbstractDescription
    base::D
    derivatives::Dict{Symbol,Int}
end
function differentiate(desc::DerivativeDescription, sym::Symbol, order=1)
    derivatives = Dict(desc.derivatives)
    if haskey(derivatives, sym)
        derivatives[sym] += order
    else
        derivatives[sym] = order
    end
    return DerivativeDescription(desc.base, derivatives)
end
differentiate(desc::AbstractDescription, sym::Symbol, order=1) = DerivativeDescription(desc, Dict(sym => order))
function Base.show(io::IO, desc::DerivativeDescription)
    derivative_strs = [(val != 1 ? "D_{$(key)}^{$(val)}" : "D_{$(key)}") for (key, val) in pairs(desc.derivatives)]
    print(io, LaTeXString("$(derivative_strs...)$(desc.base)"))
end
