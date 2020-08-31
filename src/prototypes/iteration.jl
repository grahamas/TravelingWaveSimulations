
function execute_single_modification(prototype, modification)
    mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=";")
    if mod_name == ""
        mod_name = "no_mod"
    end
    simulation = prototype(; modification...)
    execution = execute(simulation)
    if execution isa FailedExecution
        #@warn "$mod_name failed!"
        return (mod_name, missing)
    end
    if execution.solution.retcode != :Success && execution.solution.retcode != :Default && execution.solution.retcode != :Terminated
        return (mod_name, missing)
    end
    return (mod_name, execution)
end


"""
    iterate_prototype(prototype_name::AbstractString, 
                      modifications; 
                       data_root=datadir(),
                       prototype_name,
                       no_save_raw=false,
                       modifications=[]
                       no_save_raw::Bool=false,
                       max_sims_in_mem::Int=floor(Int, Sys.free_memory() / 2^20),
                       backup_paths=[],
                       delete_original::Bool=false)

Run simulation based on prototype named `prototype_name` described in `src/prototypes/prototypes.jl`.

`data_root` defines the root of the data saving directory tree.
"""
function iterate_prototype(prototype_name::AbstractString, args...; kwargs...)
    prototype = get_prototype(prototype_name)
    iterate_prototype(prototype, args...; kwargs..., prototype_name=prototype_name)
end
function iterate_prototype(prototype::Function, 
                           modifications_specification; 
                           kwargs...)
    modifications, modification_strs = parse_modifications_argument(modifications_specification)
    iterate_prototype(prototype, modifications, modification_strs; kwargs...)
end
function iterate_prototype(prototype::Function,
                           modifications::AbstractArray,
                           modification_strs::Vector{<:AbstractString}; 
                           prototype_name,
                           experiment_name="",
                           data_root::AbstractString=datadir(), 
                           max_sims_in_mem::Int=floor(Int,Sys.free_memory() / 2^20),
                           save=true, parallelize_if_possible=true) 

    @warn "# of mods: $(length(modifications))"

    # Initialize saving paths
    data_path = init_data_path(modification_strs; data_root=data_root, 
                                                  prototype_name=prototype_name,
                                                  experiment_name=experiment_name)

    # Initialize prototype
    prototype = get_prototype(prototype_name)
    
    pkeys = mods_to_pkeys(modifications)
   
    function prob_function(prob, i, repeat)
        these_mods = modifications[i]
        new_sim = prototype(; these_mods...)        
        pkeys_nt = NamedTuple{Tuple(pkeys)}([these_mods[key] for key in pkeys])
        remake(prob, p=pkeys_nt, f=convert(ODEFunction{true},make_system_mutator(new_sim)))
    end

    # Run simulation for every modification
    initial_simulation = prototype()
    sample_execution = execute(initial_simulation)
    sample_data = initial_simulation.global_reduction(sample_execution.solution)
    u_init = init_results_tuple(pkeys, sample_data)
    
    (ensemble_solver, batch_size) = if parallelize_if_possible
        if nprocs() > 1
            @warn "Parallelizing over $(nworkers()) workers"
            (EnsembleDistributed(),
             min(max_sims_in_mem, ceil(Int, length(modifications) / (nprocs() - 1))))
        elseif Threads.nthreads() > 1
            @warn "Parallelizing over $(Threads.nthreads()) threads"
            (EnsembleThreads(),
             min(max_sims_in_mem, ceil(Int, length(modifications) / (Threads.nthreads() - 1))))
        else
            @warn "Not parallelizing..."
            (EnsembleSerial(), length(modifications))
        end
    else
        (EnsembleSerial(), length(modifications))
    end
    ensemble_solution = solve(initial_simulation, ensemble_solver; 
                              prob_func=prob_function, 
                              # FIXME: consider trying setindex instead of append
                              reduction=(u,data,i) -> (append_namedtuple_arr!(u,data), false),
                              u_init=u_init, 
                              trajectories=length(modifications), 
                              batch_size=batch_size)
    (ensemble_solution, savedvals) = if ensemble_solution isa Tuple
        ensemble_solution
    else
        (ensemble_solution, nothing)
    end

    # TODO: only works for namedtuple u
    if save == true
        @warn """saving $(joinpath(data_path, "ensemble_solution.jdb"))""" 
        JuliaDB.save(table(ensemble_solution.u, pkey=pkeys), joinpath(data_path, "ensemble_solution.jdb"))
    end
    return ensemble_solution
end
