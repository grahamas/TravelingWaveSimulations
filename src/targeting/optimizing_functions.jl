
abstract type AbstractFunctionTarget{T,F<:Function} <: AbstractTarget{T} end

# Fitting function with spatiotemporal domain
struct SpatioTemporalFnTarget{T,F} <: AbstractFunctionTarget{T,F}
	fn::F
	p0::Array{T,1}
	lower::Array{T,1}
	upper::Array{T,1}
end

function spatiotemporal_input(x::AbstractArray{T,1},t::AbstractArray{T,1}) where T
	xs = repeat([x...,0], outer=(1,length(t)))
	xs[end,:] = t
	return xs
end

function target_loss(target::SpatioTemporalFnTarget{T,F}, initial_model::Model{T}, solver)
	t = time_arr(solver)
	x = space_arr(initial_model, solver)
	x_dx, pop_dx, t_dx = subsample_dxs(initial_model, solver, target.subsampler)
	x = x[x_dx]
	t = t[t_dx]
	target_data = target.fn.(spatiotemporal_input(x,t))
	return @fn (soln) -> sumsqdiff(soln[x_dx, 1, t_dx], target_data)
end

#goal: @optim_target sech2(x, t; amplitude, width, velocity) = amplitude * sech(width * (x - velocity * t)) ^ 2
# leads to:
#	target = SpatioTemporalFnTarget(sech2; ...)
#	target = Subsampled(target=target, time_subsampling = ..., space_subsampling = ...)

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

function make_arg_fns(optim_fn_default_body, optim_fn_args, param_syms, paramvec_sym, arg_fn_name, number_type)
	# Now keep in mind two dicts:
	# 	arg_fn_dict -- Outer fn, takes all three parameters, of type either T or Variable
	#	optim_fn_dict -- Inner fn, takes xt, p

	arg_fn_dict_skeleton = Dict(:name => arg_fn_name,
							   :kwargs => [],
							   :whereparams => (number_type,)
							   )
	optim_fn_default_dict = Dict(:name => Symbol(:_optim, arg_fn_name),
					   :args => optim_fn_args,
					   :kwargs => [],
					   :whereparams => (number_type,),
					   :body => optim_fn_default_body
					   )
	optim_fn_default_expr = MacroTools.combinedef(optim_fn_default_dict)
	l_varying_param_dxs = subsets(1:length(param_syms), length(param_syms)) #TODO: Remove second param when bug #30741 is resolved
	arg_fn_exprs = map(l_varying_param_dxs) do varying_param_dxs
		varying_params = param_syms[varying_param_dxs]
		free_variables_subs_dict = Dict(:($param) => :($(paramvec_sym)[$i]) for (i, param) in enumerate(varying_params))
		optim_fn_overwritten_expr = subs(optim_fn_default_expr, free_variables_subs_dict)
		arg_fn_args = map(enumerate(param_syms)) do (i, param)
			if i in varying_param_dxs
				:($(param)::AbstractVariable{$number_type})
			else
				:($(param)::$number_type)
			end
		end
		MacroTools.combinedef(merge(arg_fn_dict_skeleton, Dict(:args => arg_fn_args, :body => optim_fn_overwritten_expr)))
	end
	return arg_fn_exprs
end

function param_bounds(var::BoundedVariable)
	return (var.value, var.bounds...)
end
function param_bounds(var::UnboundedVariable)
	return (var.value, -Inf, Inf)
end
function optim_bounds(params_list...)
	varying_params_list = filter((x) -> x isa AbstractVariable, [params_list...])
	return zip(param_bounds.(varying_params_list)...)
end

macro function_target(input_fn_expr)
	input_fn_dict = splitdef(input_fn_expr)

	number_type = :T

	spacetime_sym = :xt
	paramvec_sym = :p
	optim_fn_args = [:($(sym)::AbstractArray{$number_type}) for sym in [spacetime_sym, paramvec_sym]]
	optim_fn_default_body = subs(input_fn_dict[:body], Dict(:x=>:($(spacetime_sym)[1:end-1,:]), :t=>:($(spacetime_sym)[end,:])))

	param_exprs = input_fn_dict[:kwargs]
	param_names = param_exprs

	# Define all arg-fns
	arg_fn_name = Symbol(:_, input_fn_dict[:name])
	arg_fn_exprs = make_arg_fns(optim_fn_default_body, optim_fn_args, param_names, paramvec_sym, arg_fn_name, number_type)
	# Define kw-fn
	params = input_fn_dict[:kwargs]
	kw_fn_dict = Dict(
		:name => :(SpatioTemporalFnTarget),
		:args => [:(::typeof($(input_fn_dict[:name])))],
		:kwargs => [Expr(:kw, sym, :(UnboundedVariable(0.0))) for sym in params],
		:whereparams => (),
		:body => quote
			p0, lower, upper = optim_bounds($(params...))
			SpatioTemporalFnTarget($(arg_fn_name)($(params...)), [p0...], [lower...], [upper...])
		end
	)
	return quote
		$(esc(input_fn_expr))
		$(esc.(arg_fn_exprs)...)
		$(esc(MacroTools.combinedef(kw_fn_dict)))
	end
end
