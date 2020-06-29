
export based_on_example, based_on_example_serial

function catch_for_saving(results_channel, data_path, pkeys, n_batches, progress=false)
    @show data_path
    mkpath(data_path)
    n_remaining_batches = n_batches
    counter = 1
    while n_remaining_batches > 0
        these_data = take!(results_channel)
        JuliaDB.save(table(these_data, pkey=pkeys), joinpath(data_path,"$(counter).jdb"))
        these_data = nothing
        GC.gc()
        if progress
            println("batches completed: $(counter) / $(n_batches)")
        end
        counter += 1
        n_remaining_batches -= 1
    end
end

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
        example_name::AbstractString=nothing,
        modifications::AbstractArray=[],
        max_sims_in_mem=floor(Int,Sys.free_memory() / 2^20), #Assuming a sim will never be larger than a MiB
        backup_paths=[],
        progress=false,
        delete_original=false)

    modifications, modifications_prefix = parse_modifications_argument(modifications)

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
    max_held_batches = floor(Int, max_sims_in_mem / parallel_batch_size)
    @assert max_held_batches > 0
    @warn "max_held_batches = $max_held_batches"
    pkeys = mods_to_pkeys(modifications)
    batches = Iterators.partition(modifications, parallel_batch_size)
    n_batches = length(batches)
    @show n_batches
    @show length(modifications)
    
    # Sample data to initialize results
    sample_execution = execute(example())
    sample_data = reduce_to_namedtuple(sample_execution)
    sample_results = init_results_tuple(pkeys, sample_data)
    missing_data = init_missing_data(sample_data)
    
    # Initialize results channel to receive and process output
    results_channel = RemoteChannel(() -> Channel{typeof(sample_results)}(max_held_batches))
    
    # Run simulation for every modification
    task = @distributed for modifications_batch in collect(batches)
        my_results = init_results_tuple(pkeys, sample_data)
        GC.gc()
        for modification in modifications_batch
            mod_name, execution = execute_single_modification(example, modification)
            these_params = extract_params_tuple(modification, pkeys)
            if execution !== missing #is success
                these_data = reduce_to_namedtuple(execution)
            else
                these_data = missing_data
            end
            push_namedtuple!(my_results, merge(these_params, these_data))

        end
        put!(results_channel, my_results)
    end
    success = catch_for_saving(results_channel, data_path, pkeys, n_batches, progress)
    @show fetch(task)
        
    for backup_path in backup_paths
        run(`scp -r $data_path $backup_path`)
    end
    if delete_original
        run(`rm -rf $data_path`)
    end
end

function based_on_example_serial(; data_root::AbstractString=datadir(),
        no_save_raw::Bool=false,
        example_name::AbstractString=nothing,
        modifications::AbstractArray=[],
        backup_paths=[],
        progress=false,
        delete_original=false)

    modifications, modifications_prefix = parse_modifications_argument(modifications)

    # Initialize saving paths
    data_path = joinpath(data_root, example_name, "$(modifications_prefix)$(Dates.now())_$(gitdescribe())")
    @show data_path
    if !no_save_raw
        raw_path = joinpath(data_path, "raw")
        run(`mkdir -p $raw_path`)
    else
        raw_path = nothing
    end

    # Initialize example
    example = get_example(example_name)
    
    pkeys = mods_to_pkeys(modifications)
    @show length(modifications)

    # Sample data to initialize results
    sample_execution = execute(example())
    sample_data = reduce_to_namedtuple(sample_execution)
    
    results = init_results_tuple(pkeys, sample_data)
    missing_data = init_missing_data(sample_data)
    for modification in modifications
        mod_name, execution = execute_single_modification(example, modification)
        these_params = extract_params_tuple(modification, pkeys)
        if execution !== missing #is success
            these_data = reduce_to_namedtuple(execution)
        else
            these_data = missing_data
        end
        push_namedtuple!(results, merge(these_params, these_data))
    end
    
    JuliaDB.save(table(results, pkey=pkeys), joinpath(data_path, "results.jdb"))
    if progress
        println("batches completed: $(counter) / $(n_batches)")
    end
        
    for backup_path in backup_paths
        run(`scp -r $data_path $backup_path`)
    end
    if length(backup_paths) > 0 && delete_original
        run(`rm -rf $data_path`)
    end
end
