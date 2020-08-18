"A AbstractModel specifies all parameters of a system."
abstract type AbstractModel{T,N,P} <: AbstractParameter{T} end
abstract type AbstractODEModel{T,N,P} <: AbstractModel{T,N,P} end
abstract type AbstractModelwithDelay{T,N,P} <: AbstractModel{T,N,P} end
abstract type AbstractNoisyModel{T,N,P} <: AbstractModel{T,N,P} end
abstract type AbstractNoisySpaceAction{T,N} <: AbstractSpaceAction{T,N} end
abstract type AbstractSimulation{T} <: AbstractParameter{T} end 
abstract type AbstractExecution{T,SIM} end 
# FIXME shouldn't the second type of full execution be the type of the solution?
abstract type AbstractFullExecution{T,SIM} <: AbstractExecution{T,SIM} end
export AbstractFullExecution


#function (am::Type{<:AbstractModel})(fallback_args...; fallback_kwargs...)
#    @warn """Model $am undefined!
#    ---------------------
#    $fallback_args
#    
#    $fallback_kwargs
#    ---------------------   
#    """
#    missing
#end

struct FailedSimulation{T} <: AbstractSimulation{T} end
struct FailedExecution{T,S<:FailedSimulation{T}} <: AbstractExecution{T,S}
    sim::S
end

export AbstractNoisyModel, NoisyInputModel, NoisyInputAction, AbstractODEModel

struct NoisyInputModel{T,N,P,M<:AbstractModel{T,N,P}} <: AbstractNoisyModel{T,N,P}
    model::M
    noise_process::NoiseProcess
end
function Base.getproperty(wnm::NoisyInputModel, sym::Symbol)
    if sym ∈ [:noise_process, :model]
        return getfield(wnm, sym)
    else
        return getproperty(getfield(wnm, :model), sym)
    end
end

struct NoisyInputAction{T,N,SA<:AbstractSpaceAction{T,N}} <: AbstractNoisySpaceAction{T,N}
    space_action::SA
    noise_process::NoiseProcess
end
function (act::NoisyInputAction)(du, u, p, t, W)
    du .= W
    act.space_action(du, u, p, t)
end


function (nim::NoisyInputModel)(args...) 
    inner_fn = nim.model(args...)
    NoisyInputAction(inner_fn, nim.noise_process)
end
    

n_populations(::AbstractModel{T,N,P}) where {T,N,P} = P

function initial_value(::AbstractModel{T,N,P}, space::AbstractSpace{T,N}) where {T,N,P}
    init_val = AxisArray(population_repeat(zeros(space), P), (coordinate_axes(space)..., 1:P))
    return init_val
end
    

"A Simulation holds an AbstractModel to be solved, the space on which to solve it, the time for which to solve it, the initial value, and various solver options."
struct Simulation{
        T,
        M<:AbstractModel{T},
        S<:AbstractSpace{T},
        IV<:AbstractArray{T},
        ALG,
        DT<:Union{T,Nothing},
        SV_IDX<:Union{AbstractArray,Nothing},
        SR<:Union{Tuple{Function,DataType},Nothing},
        GR<:Function} <: AbstractSimulation{T}
    model::M
    space::S
    tspan::Tuple{T,T}
    initial_value::IV
    algorithm::ALG
    dt::DT
    save_idxs::SV_IDX
    step_reduction::SR
    global_reduction::GR
    solver_options
end
function Simulation(model::M; space::S, tspan, initial_value::IV=initial_value(model,space), 
                    algorithm::ALG=nothing, dt::DT=nothing, save_idxs::SV_IDX=nothing, step_reduction::SR=nothing, global_reduction::GR=identity, 
        opts...) where {T,N,P,M<:AbstractModel{T,N,P}, S<:AbstractSpace{T,N},IV,ALG,DT,SV_IDX,SR,GR}
    save_idxs = parse_save_idxs(space, P, save_idxs)
    return Simulation{T,M,S,IV,ALG,DT,typeof(save_idxs),SR,GR}(model, space, tspan, initial_value, algorithm, dt, save_idxs, step_reduction, global_reduction, opts)
end
function Simulation(model::Missing; kwargs...)
    return FailedSimulation{Missing}()
end

"An Execution holds a Simulation and the solution obtained by running the Simulation."
struct Execution{T,S<:AbstractSimulation{T},D<:DESolution} <: AbstractFullExecution{T,S}
    simulation::S
    solution::D
end
struct ReducedExecution{T,ST,S<:AbstractSimulation{T},SV<:SavedValues{T,ST}} <: AbstractExecution{T,S}
    simulation::S
    saved_values::SV
end
struct AugmentedExecution{T,ST,S<:AbstractSimulation{T},D<:DESolution,SV<:SavedValues{T,ST}} <: AbstractFullExecution{T,S}
    simulation::S
    solution::D
    saved_values::SV
end
make_execution(s::S, sol::DESolution) where {T,S <: Simulation{T}} = Execution(s,sol)
make_execution(s::S, sv::SavedValues) where {T,S <: Simulation{T}} = ReducedExecution(s,sv)
make_execution(s::S, (sol,sv)::Tuple{<:DESolution,<:SavedValues}) where {T,S <: Simulation{T}} = AugmentedExecution(s,sol,sv)
export ReducedExecution, AugmentedExecution,make_execution
function execute(s::Simulation)
    return make_execution(s, solve(s))
end
function execute(s::FailedSimulation)
    return FailedExecution(s)
end
            
coordinates(sim::Simulation) = coordinates(space(sim))
coordinates(ex::AbstractExecution) = coordinates(space(ex))
timepoints(ex::AbstractFullExecution) = ex.solution.t
reduced_space(sim::Simulation) = subsample(sim.space, sim.save_idxs)
reduced_space(ex::AbstractExecution) = reduced_space(ex.simulation)
frame_xs(exec::AbstractExecution{T,<:Simulation{T,<:M}}) where {T,M<:AbstractModel{T,1}} = [x[1] for x in reduced_space(exec).arr]
origin_idx(sim::Simulation) = origin_idx(sim.space)
origin_idx(ex::AbstractExecution) = origin_idx(ex.simulation)
extrema(ex::AbstractFullExecution) = extrema(ex.solution.u)
stimulus_center(mod::AbstractModel) = center(mod.stimulus)
stimulus_center(sim::Simulation) = stimulus_center(sim.model)
export stimulus_center

"""
    make_system_mutator(simulation)

Construct the differential function to be provided to the ODE solver.
"""
function make_system_mutator(sim::SIM) where {SIM <: Simulation}
    sim.model(sim.space)
end

"""
    generate_problem(simulation)

Return an ODEProblem of the `simulation.model` with time span specified by `simulation.solver`.
"""
function generate_problem(simulation::Simulation{T,<:AbstractODEModel}; callback=nothing) where {T}
    system_fn! = make_system_mutator(simulation)# simulation.model(simulation.space)
    ode_fn = convert(ODEFunction{true}, system_fn!)
    return ODEProblem(ode_fn, simulation.initial_value, simulation.tspan, callback=callback)
end

function generate_problem(simulation::Simulation{T,<:AbstractNoisyModel}; callback=nothing) where {T}
    system_fn! = make_system_mutator(simulation)# simulation.model(simulation.space)
    return RODEProblem(system_fn!, simulation.initial_value, simulation.tspan, noise=simulation.model.noise_process, noise_prototype=zeros(size(simulation.initial_value)...), callback=callback)
end

# TODO: Add history functionality
# function generate_problem(simulation::Simulation{T,MwD}) where {T, MwD<:AbstractModelwithDelay}
#     system_mutator! = make_system_mutator(simulation)
#     return DDEProblem(system_mutator!, simulation.initial_value, history(simulation), simulation.tspan)
# end
parse_save_idxs(::Any, P, save_idx) = save_idx
function parse_save_idxs(space::AbstractSpace, P, subsampler::Union{AbstractSubsampler,AbstractArray{<:AbstractSubsampler}})
	one_pop_coordinates = coordinate_indices(space, subsampler)
	population_coordinates(one_pop_coordinates, P)
end

# TODO this could probably all be better served by dispatch; but that gets difficult
function _solve(simulation::Simulation, alg; callback=nothing, dt=nothing, solver_options...)
    if dt != nothing
        solver_options = (solver_options..., dt=dt)
    end
    problem = generate_problem(simulation; callback=callback)
    solve(problem, alg; solver_options...)
end
function _solve(simulation::Simulation, alg, ensemble_alg; prob_func::Function, output_func=(sol, i) -> (simulation.global_reduction(sol), false), reduction=(u,data,i) -> (append!(u,data), false), callback=nothing, dt=nothing, u_init=[], solver_options...)
    if dt != nothing
        solver_options = (solver_options..., dt=dt)
    end
    initial_problem = generate_problem(simulation; callback=callback)
    ensemble_problem = EnsembleProblem(initial_problem;
                    output_func=output_func,
                    prob_func=prob_func,
                    reduction=reduction,
                    u_init=u_init)
    solve(ensemble_problem, alg, ensemble_alg; solver_options...)
end

function solve(simulation::Simulation{T,M,S,IV,ALG,DT,SV_IDX,Nothing,GR}, args...; kwargs...) where {T,M,S,IV,ALG,DT,SV_IDX,GR}
    sol = _solve(simulation, simulation.algorithm, args...; dt=simulation.dt, save_idxs=simulation.save_idxs, simulation.solver_options..., kwargs...)
    return sol
end

function solve(simulation::Simulation{T,M,S,IV,ALG,DT,SV_IDX,<:Tuple,GR}, args...; kwargs...) where {T,M,S,IV,ALG,DT,SV_IDX,GR}
    save_func, save_type = simulation.step_reduction
    sv = SavedValues(T,save_type)
    all_callbacks = SavingCallback(save_func, sv)
        
    if :callback ∈ keys(simulation.solver_options)
        all_callbacks = CallbackSet(simulation.solver_options[:callback], all_callbacks)
    end
    # save_callback has all callbacks
    
    sol = _solve(simulation, simulation.algorithm, args...; dt=simulation.dt, save_idxs=simulation.save_idxs, simulation.solver_options..., callback=all_callbacks, kwargs...)
    return (sol, sv)
end

