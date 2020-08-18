
"""
For `obj isa AbstractSearch{T,P}`, the method `search(obj)` finds the realization of obj.parameter
"""
@with_kw struct Search{T, PRM <: AbstractParameter{<:MaybeVariable{T}}, TGT <: AbstractTarget{T}}
    parameter::PRM
    target::TGT
end

struct SearchExecution
    search::Search
    result::Execution
end
function execute(search::Search{T,<:Simulation}) where T
    varying_simulation = search.parameter
    result = search(varying_simulation, target; kwargs...)
    result_simulation = param_from_p(varying_simulation, minimizer(result))
    SearchExecution(search, execute(result_simulation))
end

minimizer(res::BlackBoxOptim.OptimizationResults) = best_candidate(res)
minimizer(res::Optim.OptimizationResults) = Optim.minimizer(res)

result_simulation(search_execution::SearchExecution) = search_execution.result.simulation

# * Run the search

make_problem_generator() = error("undefined.")

function _build_loss_objective(initial_problem, solver::Solver{T,Euler}, space, loss_fn, prob_generator) where T
    build_loss_objective(initial_problem, Euler(), loss_fn;
        prob_generator=prob_generator, dt=solver.simulated_dt,
        save_at=saved_dt(solver), save_idxs=save_idxs(solver, space))
end

function _build_loss_objective(initial_problem, solver::Solver{T,Nothing}, space, loss_fn, prob_generator) where T
    initial_problem = problem_generator(nothing, initial_values)
    build_loss_objective(initial_problem, Tsit5(), loss_fn;
        prob_generator=prob_generator,
        save_at=saved_dt(solver),
        timeseries_steps=solver.time_save_every,
        save_idxs=save_idxs(solver, space),
        alg_hints=[solver.stiffness])
end

function search(varying_simulation::Simulation{MaybeVariable{T}}, target::AbstractTarget{T}; MaxSteps=11e3) where T
    varying_deconstruction = VariableDeconstruction(varying_simulation)
    problem_generator = make_problem_generator(varying_deconstruction)
    initial_problem = problem_generator(nothing, initial_values)
    loss_fn = target_loss(target, varying_simulation)
    loss_obj = _build_loss_objective(initial_problem, varying_simulation.solver, loss_fn, problem_generator)
    result = bboptimize(loss_obj; NumDimensions=length(varying_deconstruction.initial_values),
        MaxSteps=MaxSteps, SearchRange=varying_deconstruction.bounds)
    return result
end

# function run_search_optim(varying_model, variable_map, solver, target, initial_p, p_bounds::Array{<:Tuple}; MaxSteps=11e3)
#     @show MaxSteps
#     initial_model = model_from_p(varying_model, variable_map, initial_p)
#     problem_generator = make_problem_generator(initial_model, solver, variable_map)
#     initial_problem = problem_generator(nothing, initial_p)
#     loss_fn = target_loss(target, initial_model, solver)
#     loss_obj = _build_loss_objective(initial_problem, solver, initial_model.space, loss_fn,
#                                                 problem_generator)
#     lower = [bounds[1] for bounds in p_bounds]
#     upper = [bounds[2] for bounds in p_bounds]
#     result = optimize(loss_obj, lower, upper, initial_p, Fminbox(ParticleSwarm()),
#                         Optim.Options(
#                          iterations = 10,
#                          show_every = 1))
#     return result
# end

function run_search(jl_filename::AbstractString)
    include(jl_filename)
    filecopy(p_search.output, jl_filename, basename(jl_filename))
    return p_search
end

function make_problem_generator(deconstructed_simulation::VariableDeconstruction{SIM}) where {SIM <: Simulation}
    function problem_generator(prob, values_vec::AbstractArray)
        simulation = reconstruct_with_vec(deconstructed_simulation, values_vec)
        system_function = make_system_function(simulation)
        ode_fn = convert(ODEFunction{true}, system_function)
        return ODEProblem(ode_fn, initial_value(simulation), time_span(simulation), values_vec)
    end

    return problem_generator
end
