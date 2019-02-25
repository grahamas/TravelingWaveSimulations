
abstract type AbstractFunctionTarget{T,F<:Function} <: AbstractTarget{T} end

# Fitting function with spatiotemporal domain
@with_kw struct SpatioTemporalFnTarget{T,F} <: AbstractFunctionTarget{T,F}
	fn::F
	p0::Array{T,1}
	lower::Array{T,1}
	upper::Array{T,1}
	time_subsampler::Subsampler
	space_subsampler::Subsampler
end

function spatiotemporal_input(x::AbstractArray{T,1},t::AbstractArray{T,1}) where T
	xt = repeat([x...,0], outer=(1,length(t)))
	xt[end,:] = t
	return xt
end

function target_loss(target::SpatioTemporalFnTarget{T,F}, initial_model::Model{T}, solver) where {T,F}
	t = time_arr(solver)
	x = space_arr(initial_model, solver)
	x_dx, pop_dx, t_dx = subsampling_idxs(initial_model, solver, target.time_subsampler, target.space_subsampler)
	x = x[x_dx]
	t = t[t_dx]
	xt = spatiotemporal_input(x, t)
	return function optim_fit(soln)
		fit_loss(p) = sumsqdiff(soln[x_dx, 1, t_dx], target.fn(xt,p))
		result = optimize(fit_loss, target.lower, target.upper, target.p0, Fminbox(Optim.ConjugateGradient()))
		@show Optim.minimizer(result)
		return Optim.minimum(result)
	end
end

### REDUX:  Instead of calling the macro before, call the macro as part of the
### 		target call.

function parse_xt_args(exprs)
	@assert length(exprs) == 2
	@assert (@capture exprs[1] (x | (x ∈ (xlower_,xupper_,dx_)) | (x ∈ (xlower_,xupper_))))# | (x ∈ xwindow_ Δ xΔ_))
	@assert (@capture exprs[2] (t | (t ∈ (tlower_,tupper_,dt_)) | (t ∈ (tlower_,tupper_))))# | (t ∈ twindow_ Δ tΔ_))
	x_subsampler = :(Subsampler(; window = ($xlower,$xupper), Δ = $dx))#, Δ = $xΔ))
	t_subsampler = :(Subsampler(; window = ($tlower,$tupper), Δ = $dt))
	return (x_subsampler, t_subsampler)
end

function parse_params(exprs)
	free_param_names = []
	lower_bounds = []
	upper_bounds = []
	guesses = []
	fixed_param_mapping = []

	for expr in exprs
		if @capture(expr, freename_ ∈ (lower_, upper_) : guess_)
			push!(free_param_names, freename)
			push!(lower_bounds, lower)
			push!(upper_bounds, upper)
			push!(guesses, guess)
		elseif @capture(expr, fixedname_ = value_)
			push!(fixed_param_mapping, fixedname => value)
		elseif expr isa Symbol
			push!(free_param_names, expr)
			push!(lower_bounds, :(-Inf))
			push!(upper_bounds, :(Inf))
			push!(guesses, :(0.0))
		else
			error("invalid parameter argument structure: $expr")
		end
	end
	return (free_param_names, lower_bounds, upper_bounds, guesses, fixed_param_mapping)
end

macro function_target(input_fn_expr)
	input_fn_dict = MacroTools.splitdef(input_fn_expr)

	number_type = :T

	spacetime_sym = :xt
	paramvec_sym = :p
	fn_args = [:($(sym)::AbstractArray{$number_type}) for sym in [spacetime_sym, paramvec_sym]]

	x_subsampler, t_subsampler = parse_xt_args(input_fn_dict[:args])

	free_param_names, lower_bounds, upper_bounds, guesses, fixed_param_mapping = parse_params(input_fn_dict[:kwargs])
	free_param_mapping = Dict(:($param) => :($(paramvec_sym)[$i]) for (i, param) in enumerate(free_param_names))

	fn_body = subs(input_fn_dict[:body], Dict(:x=>:($(spacetime_sym)[1:end-1,:]), :t=>:($(spacetime_sym)[end,:]), free_param_mapping..., fixed_param_mapping...))

	fn_def = MacroTools.combinedef(Dict(
		:name => input_fn_dict[:name],
		:args => fn_args,
		:kwargs => [],
		:whereparams => (number_type,),
		:body => fn_body
	))

	return quote
		SpatioTemporalFnTarget(
			fn=$fn_def,
			p0=[$(guesses...)],
			lower=[$(lower_bounds...)],
			upper=[$(upper_bounds...)],
			space_subsampler=$x_subsampler,
			time_subsampler=$t_subsampler
		)
	end
end
