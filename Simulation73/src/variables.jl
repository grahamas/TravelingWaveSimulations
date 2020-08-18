

abstract type AbstractVariable{T<:Real} end

"A parameter to vary, without bounds."
struct UnboundedVariable{T} <: AbstractVariable{T}
    value::T
end

struct BoundedVariable{T} <: AbstractVariable{T}
    value::T
    bounds::Tuple{T,T}
end

function BoundedVariable(value::T) where {T <: Real}
	ten_pct = value / 10
	BoundedVariable{T}(value, (value-ten_pct,value+ten_pct))
end

const MaybeVariable{T} = Union{T,AbstractVariable{T}}

Base.zero(::Type{<:MaybeVariable{T}}) where T <: Number = zero(T)
default_value(var::V) where V <: AbstractVariable = var.value
bounds(var::V) where V <: BoundedVariable = var.bounds

"""
    deconstruct(ad::AbstractDeconstructor, param)

Deconstruct a parameter into nested tuples of the form `(type,value)`.

```
    Parameter{Varying{T}} |> deconstruct |> reconstruct # isa Parameter{Varying{T}}
```
"""
function deconstruct(m::M, ad::AbstractDeconstructor=Deconstructor()) where {M <: AbstractParameter}
    deconstruction = Tuple{Type,Any}[]
    for i_field in 1:nfields(m)
        substruct = ad(getfield(m, i_field))
        push!(deconstruction, substruct)
    end
    return (typeof(m), deconstruction)
end

"""
	deconstruct(..., DeconstructorFixingVariables())

Deconstructs as normal (i.e. [Deconstructor](@ref)) except on [AbstractVariable](@ref) parameters, which are fixed to their default value.
"""
struct DeconstructorFixingVariables <: AbstractDeconstructor end

function deconstruct(val::AbstractVariable, ad::AbstractDeconstructor=Deconstructor())
    return (typeof(val), val)
end
function deconstruct(v::AbstractVariable, deconstructor::DeconstructorFixingVariables)
    val = default_value(v)
    return (typeof(val), val)
end


function base_type(::Type{P}) where {P <: AbstractParameter}
    BTs = map(base_type, P.parameters)
    return (P.name.wrapper){BTs...}
end

function base_type(::Type{V}) where {T, V<: Union{AbstractVariable{T}, T}}
    return T
end

struct VariableDeconstruction{OBJ} <: AbstractDeconstruction{OBJ}
    source::OBJ
    deconstruction::Tuple{Type,Array}
    dxs_map
    default_values
    bounds
end

function VariableDeconstruction(obj)
	dxs_map, default_values, bounds = variables(obj)
	deconstruction = deconstruct(obj, DeconstructorFixingVariables())
	VariableDeconstruction(obj, deconstruction, dxs_map, default_values, bounds)
end

"""Takes model with variable parameters,
and returns default variable values and indices to those variables."""
function variables(variable_parameter::P) where {T,
											P <: AbstractParameter{<:MaybeVariable{T}}}
    deconstructed = deconstruct(variable_model)
    default_value, dxs_map, bounds = variables(deconstructed, T)
    return dxs_map, default_value, bounds
end

"Initialize all variables to their default value."
function variables(deconstructed::Tuple{Type,<:AbstractArray}, T::Type)
    dxs_map = []
    default_values = T[]
    bounds = Tuple{T,T}[]
    for dx in CartesianIndices(size(deconstructed[2]))
        (typ, val) = deconstructed[2][dx]
        if typ <: AbstractVariable
            push!(dxs_map, [dx])
            push!(default_values, default_value(val))
            push!(bounds, bounds(val))
        elseif val isa AbstractArray
            dxs, ps, bds = variables((typ, val), T)
            @assert length(ps) == length(dxs)
            push!(dxs_map, [vcat(dx, x) for x in dxs]...)
            push!(default_values, ps...)
            push!(bounds, bds...)
        end
    end
    return dxs_map, default_values, bounds
end

function param_from_vec(var_deconstruction::VariableDeconstruction, values_vec)
    target = deepcopy(var_deconstruction.deconstruction)
	for (dxs, value) in zip(var_deconstruction.dxs_map, values_vec)
		set_deep_dx!(target, dxs, value)
	end
	reconstruct(target)
end


function param_from_vec(varying_param::AbstractParameter{MaybeVariable{T}}, values_vec) where T
    varying_param_deconstruction = VariableDeconstruction(varying_param)
    param_from_vec(varying_param_deconstruction, values_vec)
end
