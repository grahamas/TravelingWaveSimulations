
export based_on_example, based_on_example_serial

function execute_single_modification(example, modification)
    mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=";")
    if mod_name == ""
        mod_name = "no_mod"
    end
    simulation = example(; modification...)
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
    based_on_example(; data_root=datadir(),
                       no_save_raw=false,
                       example_name=nothing,
                       modifications=[])

Run simulation based on example named `example_name` described in `src/examples`.

`modifications` is an array of *strings* modifying the example's built-in parameters. These strings can 1) name a modifications file in `scripts/modifications`, 2) modify a keyword value using the syntax `name=value` or 3) describe a range of values for a keyword, using colon syntax (`name=start:step:end`). Valid keyword targets are *any* keywords in the signature or body of the example, so long as that keyword is *unique*.

`no_save_raw` overrides the default behavior of saving the raw simulation solution when set to `true`.

`data_root` defines the root of the data saving directory tree.

# Example
```jldoctest
julia> using TravelingWaveSimulations
julia> based_on_example(; example_name="sigmoid_normal", modifications=["iiS=0.7"])
```
"""
function based_on_example(; data_root::AbstractString=datadir(), 
        no_save_raw::Bool=false,
        example_name::AbstractString,
        modifications::AbstractArray=[],
        max_sims_in_mem=floor(Int,Sys.free_memory() / 2^20), #Assuming a sim will never be larger than a MiB
        backup_paths=[],
        progress=false,
        delete_original=false)

    modifications, modifications_prefix = parse_modifications_argument(modifications)
    @warn "# of mods: $(length(modifications))"

    # Initialize saving paths
    data_path = joinpath(data_root, example_name, "$(modifications_prefix)$(Dates.now())_$(gitdescribe())")
    if !no_save_raw
        raw_path = joinpath(data_path, "raw")
        run(`mkdir -p $raw_path`)
    else
        raw_path = nothing
    end

    # Initialize example
    example = get_example(example_name)
    
    # Initialize parallelism
    parallel_batch_size = if nprocs() == 1
        @warn "MUST HAVE MORE THAN ONE CORE"
        throw(InvalidStateException)
    else
        min(max_sims_in_mem, ceil(Int, length(modifications) / (nprocs() - 1)))
    end
    @show parallel_batch_size
    pkeys = mods_to_pkeys(modifications)
   
    function prob_func(prob, i, repeat)
        these_mods = modifications[i]
        new_sim = example(; these_mods...)        
        pkeys_nt = NamedTuple{Tuple(pkeys)}([these_mods[key] for key in pkeys])
        remake(prob, p=pkeys_nt, f=convert(ODEFunction{true},make_system_mutator(new_sim)))
    end

    # Run simulation for every modification
    initial_simulation = example()
    sample_execution = execute(initial_simulation)
    sample_data = initial_simulation.global_reduction(sample_execution.solution)
    u_init = init_results_tuple(pkeys, sample_data)
    ensemble_solution = solve(initial_simulation, EnsembleDistributed(); 
                              prob_func=prob_func, 
                              reduction=(u,data,i) -> (append_namedtuple_arr!(u,data), false),
                              u_init=u_init, 
                              trajectories=length(modifications), 
                              batch_size=parallel_batch_size)

    return ensemble_solution
    # for backup_path in backup_paths
    #     run(`scp -r $data_path $backup_path`)
    # end
    # if delete_original
    #     run(`rm -rf $data_path`)
    # end
end

function based_on_example_serial(; data_root::AbstractString=datadir(), 
        no_save_raw::Bool=false,
        example_name::AbstractString,
        modifications::AbstractArray=[],
        backup_paths=[],
        progress=false,
        delete_original=false)

    modifications, modifications_prefix = parse_modifications_argument(modifications)
    @warn "# of mods: $(length(modifications))"

    # Initialize saving paths
    data_path = joinpath(data_root, example_name, "$(modifications_prefix)$(Dates.now())_$(gitdescribe())")
    if !no_save_raw
        raw_path = joinpath(data_path, "raw")
        run(`mkdir -p $raw_path`)
    else
        raw_path = nothing
    end

    # Initialize example
    example = get_example(example_name)
    
    pkeys = mods_to_pkeys(modifications)
   
    function prob_func(prob, i, repeat)
        these_mods = modifications[i]
        new_sim = example(; these_mods...)        
        pkeys_nt = NamedTuple{Tuple(pkeys)}([these_mods[key] for key in pkeys])
        remake(prob, p=pkeys_nt, f=convert(ODEFunction{true},make_system_mutator(new_sim)))
    end

    # Run simulation for every modification
    initial_simulation = example()
    sample_execution = execute(initial_simulation)
    sample_data = initial_simulation.global_reduction(sample_execution.solution)
    u_init = init_results_tuple(pkeys, sample_data, length(modifications))
    @show typeof(u_init)
    ensemble_solution = solve(initial_simulation, EnsembleSerial(); 
                              prob_func=prob_func, 
                              reduction=(u,data,i) -> (setindex_namedtuple!(u, data, i), false),
                              u_init=u_init, 
                              trajectories=length(modifications))

    return ensemble_solution
    # for backup_path in backup_paths
    #     run(`scp -r $data_path $backup_path`)
    # end
    # if delete_original
    #     run(`rm -rf $data_path`)
    # end
end
