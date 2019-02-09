
abstract type AbstractExampleTarget{T} <: AbstractTarget{T} end
abstract type AbstractFunctionTarget{T} <: AbstractTarget{T} end

struct SpatioTemporalFnTarget{F<:Function} <: AbstractFunctionTarget{T} end


function make_target_fn(...)
	@fn (xt, p) -> begin
		subs(body, x=xt[1:end-1,:], t=xt[end,:])
	end
end

#region Sech2Target

#goal @optim_target sech2(x, t; amplitude, width, velocity) = amplitude * sech(width * (x - velocity * t)) ^ 2

function params_as_args_expr(params, type)
	[:($(param)::MaybeVariable{$type}=UnboundVariable()) for param in params]
end

function optim_fn_conditional_expr(; name, fn_args, varying_params, fixed_params, fn_body, generic_fn_dict)
	fn_name = gensym(name)
	fn_def = combine_def(merge(generic_fn_dict, Dict(:name => fn_name, :args => fn_args, :body => fn_body)))
	condition_exprs = [:(isnothing($param)) for param in varying_params]
	quote
		if all([$(condition_exprs...)])
			$fn_def
			return $fn_name
		end
	end
end

# Definitions
# kw-fn means function that takes kwargs
# arg-fn means function that only takes args
# optim-fn means function that optim takes (i.e. args xt and p)
# Method:
# Define kw-fn which returns an optim-fn
# kw-fn passes to inner arg-fn in canonical order
# An arg-fn is defined for every possible combination of input types
#	each arg being either Variable or Number
#	Variables get p[i] in returned optim-fn
#	Numbers get fixed

### Generic arg-fn returning optim-fn
# function arg-fn(arg1::Variable, arg2::NotVariable, arg3::Variable)
#	function fn-name-gensym(xt, p)
#		body with xt, p[i], and arg2
#	end
# end

### Generic kw-fn calling arg-fn and returning optim-fn
# function kw-fn(; arg1=UnboundVariable(), arg2=UnboundVariable(), arg3=UnboundVariable())
#	arg-fn(arg1, arg2, arg3)
# end

function optim_st_target(fn_expr)
	fn_dict = splitdef(fn_expr)
	completely_specified_subs_dict = Dict(:x=>:(xt[1:end-1,:]), :t=>:(xt[end,:]), [:($arg) => :($arg) for arg in fn_dict[:kwargs]]...)
	param_subset_dxs = subsets(1:length(fn_dict[:kwargs]))
	factory_definition_exprs = map(param_subset_dxs) do subset_dx
		varying_params = params[subset_dx]
		this_subs_dict = copy(subs_dict)
		for (i, varying_param) in enumerate(varying_params)
			this_subs_dict[varying_param] = :(p[$i])
		end
		this_body = subs(body; this_subs_dict...)
		optim_fn_conditional_expr(; fn_args=[:xt, :p], varying_params=varying_params, fixed_params, fn_body=this_body, fn_dict)
	end
	type = :T

	end
	constructor_expr = quote
		function SpatioTemporalFnTarget(fn::typeof($name); $(kwargs...))
			param_vec = [$(kwargs...)]
			target_fn = make_target_fn($(kwargs...))
			inital_p = default_value.(param_vec)
			lower, upper = bounds(param_vec)
	end
end


#endregion

#region MatchExample

struct MatchExample{T} <: AbstractExampleTarget{T}
	data::Array{T}
	x::Array{T,1}
	t::Array{T,1}
end
function MatchExample(; file_name::String="")
	@load file_name wave x t
	MatchExample(wave, x, t)
end

function (p::MatchExample{T})(soln::AbstractArray{T}) where T#, x_dxs::AbstractArray{Int,1}, pop_dxs::AbstractArray{Int,1}, t_dxs::AbstractArray{Int,1}) where {T}
	res = sum((soln .- p.data) .^ 2)
	return res
end

function target_loss(fn::MatchExample{T}, model::WCMSpatial1D{T}, solver::Solver{T}) where T
	t_target, x_target = fn.t, fn.x
	t_dxs = subsampling_time_idxs(t_target, solver)
	x_dxs = subsampling_space_idxs(x_target, model, solver)
	pop_dxs = 1
	(soln) -> (fn(soln[x_dxs, pop_dxs, t_dxs]))
end

#endregion

#region StretchExample

struct StretchExample{T} <: AbstractExampleTarget{T}
	data::Array{T}
	x::Array{T,1}
	t::Array{T,1}
	stretch_dx::Int
end
function StretchExample(; file_name::String="", stretch_dx=nothing)
	@load file_name wave x t
	StretchExample(wave, x, t, stretch_dx)
end

function (p::StretchExample{T})(soln::AbstractArray{T}, x_dxs::AbstractArray{Int,1}, pop_dxs::AbstractArray{Int,1}, t_dxs::AbstractArray{Int,1}) where {T}
	sum((soln[x_dxs, pop_dxs, t_dxs] .- p.data) .^ 2)
end

function target_loss(fn::StretchExample{T}, model::WCMSpatial1D{T}, solver::Solver{T}) where T
	t_target, x_target = fn.t, fn.x
	t_dxs = subsampling_time_idxs(t_target, solver)
	x_dxs = subsampling_space_idxs(x_target, model, solver)
	pop_dxs = [1]
	(soln) -> fn(soln, x_dxs, pop_dxs, t_dxs)
end

#endregion
