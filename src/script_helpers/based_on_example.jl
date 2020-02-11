
export based_on_example

DEFAULT_SAVE_BATCH_SIZE = 10000
"""
    based_on_example(; data_root=datadir(),
                       no_save_raw=false,
                       example_name=nothing,
                       modifications=[],
                       analyses=[])

Run simulation based on example named `example_name` described in `src/examples`.

`modifications` is an array of *strings* modifying the example's built-in parameters. These strings can 1) name a modifications file in `scripts/modifications`, 2) modify a keyword value using the syntax `name=value` or 3) describe a range of values for a keyword, using colon syntax (`name=start:step:end`). Valid keyword targets are *any* keywords in the signature or body of the example, so long as that keyword is *unique*.

`analyses` is an array of *strings* naming files (sans `.jl`) found in `scripts/analyses`.

`no_save_raw` overrides the default behavior of saving the raw simulation solution when set to `true`.

`data_root` defines the root of the data saving directory tree.

# Example
```jldoctest
julia> using TravelingWaveSimulations
julia> based_on_example(; example_name="sigmoid_normal", analyses=["radial_slice"], modifications=["iiS=0.7"])
```
"""
function based_on_example(; data_root::AbstractString=datadir(), no_save_raw::Bool=false,
        example_name::AbstractString=nothing,
        modifications::AbstractArray=[],
        analyses::AbstractArray=[],
        batch=DEFAULT_SAVE_BATCH_SIZE,
        max_sims_in_mem=nothing,
        backup_paths=[],
        progress=false,
        force_parallel=false)

    modifications, modifications_prefix = parse_modifications_argument(modifications)
    analyses, analyses_prefix = parse_analyses_argument(analyses)
    @show analyses

    # Initialize saving paths
    if length(analyses) > 0
        analyses_path = joinpath(plotsdir(), example_name, "$(modifications_prefix)$(analyses_prefix)_$(Dates.now())_$(gitdescribe())")
        mkpath(analyses_path)
    else
        analyses_path = ""
    end
    if !no_save_raw
        data_path = joinpath(data_root, example_name, "$(modifications_prefix)$(Dates.now())_$(gitdescribe())")
    else
        data_path = nothing
    end

    example = get_example(example_name)
    if !force_parallel && (!(@isdefined nprocs) || nprocs() == 1)
        @warn "Not parallelizing parameter sweep."
        failures = execute_modifications_serial(example, modifications, analyses, data_path, analyses_path, no_save_raw, progress)
        @show failures
    else
        @warn "Parallelizing parameter sweep."
        if no_save_raw
            failures = execute_modifications_parallel_nosaving(example, modifications, analyses, data_path, analyses_path, batch, progress)
        else
            failures = execute_modifications_parallel_saving(example, modifications, analyses, data_path, analyses_path, batch, max_sims_in_mem, progress)
        end
        @show failures
    end
    for backup_path in backup_paths
        run(`scp -r $data_path $backup_path`)
    end
end

function getkeys(d, keys)
    [d[key] for key in keys]
end

function extract_data_namedtuple(execution::AbstractExecution)
    soln = execution.solution
    u = soln.u
    t = soln.t
    x = coordinates(space(execution)) |> collect
    return (u=u, t=t, x=x)
end
function extract_data_namedtuple(execution::Missing)
    return (u=missing, t=missing, x=missing)
end


function extract_params_tuple(modification, pkeys)
    return NamedTuple{Tuple(pkeys)}(getkeys(modification, pkeys))
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
    if execution.solution.retcode == :Unstable
        @warn "$mod_name unstable!"
        return (mod_name, missing)
    end
    return (mod_name, execution)
end
function mods_to_pkeys(modifications)
    pkeys = keys(modifications[1]) |> collect
    disallowed_keys = [:algorithm, :u, :x, :t, :n, :n_points, :extent, :save_idxs]
    pkeys = filter((x) -> !(x in disallowed_keys), pkeys)
    pkeys = filter((x) -> !(modifications[1][x] === nothing), pkeys)
end
function execute_modifications_serial(example, modifications::Array{<:Dict}, analyses,
        data_path::String, analyses_path::String, no_save_raw, progress=false)
    # FIXME: Needs to incorporate batching
    pkeys = mods_to_pkeys(modifications)
    data_namedtuple = nothing
    failures = nothing
    for modification in modifications
        mod_name, execution = execute_single_modification(example, modification)
        these_params = extract_params_tuple(modification, pkeys)
        if execution === missing
            push_namedtuple!(failures, these_params)
            continue
        end
        if !no_save_raw
            these_data_namedtuple = extract_data_namedtuple(execution)
            data_namedtuple = push_namedtuple!(data_namedtuple, merge(these_data_namedtuple, these_params))
        end
        analyse.(analyses, Ref(execution), analyses_path, mod_name)
    end
    if !no_save_raw      
        mkpath(data_path)
        results_table = table(data_namedtuple, pkey=pkeys)
        JuliaDB.save(results_table, joinpath(data_path, "raw_data.csv"))
        failures !== nothing && JuliaDB.save(table(failures), joinpath(data_path, "failures.csv"))
    end
    return failures
end
function catch_for_saving(results_channel, data_path, pkeys, n_batches, progress=false)
    @show data_path
    mkpath(data_path)
    n_remaining_batches = n_batches
    all_failures = nothing
    counter = 1
    while n_remaining_batches > 0
        (these_failures, these_data) = take!(results_channel)
        @show counter
        all_failures = push_namedtuple!(all_failures, these_failures)
        JuliaDB.save(table(these_data, pkey=pkeys), joinpath(data_path,"$(counter).jdb"))
        these_data = nothing; these_failures = nothing
        GC.gc()
        if progress
            println("batches completed: $(counter) / $(n_batches)")
        end
        counter += 1
        n_remaining_batches -= 1
    end
    all_failures !== nothing && JuliaDB.save(table(all_failures), joinpath(data_path, "failures.jdb"))
end

function init_results_tuple(pkeys::Array{Symbol}, results::NamedTuple{NAMES,TYPES}) where {NAMES,TYPES}
    N = length(pkeys)
    ALL_NAMES = Tuple([pkeys...,NAMES...])
    ARR_TYPES = Tuple{[Array{Float64,1} for _ in 1:N]..., [Array{Union{Missing,TYPE},1} for TYPE in TYPES.parameters]...}
    empty_arrs = [[Float64[] for _ in 1:N]..., [Union{Missing,TYPE}[] for TYPE in TYPES.parameters]...]
    NamedTuple{ALL_NAMES,ARR_TYPES}(empty_arrs)
end

function init_failures_tuple(pkeys::Array{Symbol})
    N = length(pkeys)
    NamedTuple{Tuple(pkeys),NTuple{N,Array{Float64,1}}}([Float64[] for _ in 1:N])
end

Base.getindex(nt::NamedTuple, dx::Array{Symbol}) = getindex.(Ref(nt), dx)
function execute_modifications_parallel_saving(example, modifications::Array{<:Dict}, analyses,
        data_path::String, analyses_path::String, max_batch_size, max_sims_in_mem::Int, progress=false)
    parallel_batch_size = nprocs() == 1 ? max_batch_size : min(max_batch_size, ceil(Int, length(modifications) / (nprocs() - 1)))
    pkeys = mods_to_pkeys(modifications)
    max_held_batches = floor(Int, max_sims_in_mem / parallel_batch_size)
    @assert max_held_batches > 0
    @warn "max_held_batches = $max_held_batches"
    batches = Iterators.partition(modifications, parallel_batch_size)
    n_batches = length(batches)
    @show n_batches
    results_channel = RemoteChannel(() -> Channel{Tuple}(max_held_batches))
    failures = @spawnat :any catch_for_saving(results_channel, data_path, pkeys, n_batches, progress)
    sample_execution = execute(example())
    sample_data = extract_data_namedtuple(sample_execution)
    task = @distributed for modifications_batch in collect(batches)
        results = init_results_tuple(pkeys, sample_data)
        failures = init_failures_tuple(pkeys)
        GC.gc()
        for modification in modifications_batch
            mod_name, execution = execute_single_modification(example, modification)
            @show mod_name
            these_params = extract_params_tuple(modification, pkeys)
            if execution !== nothing #is success
                these_data = extract_data_namedtuple(execution)
                @show typeof(results)
                results = push_namedtuple!(results, merge(these_params, these_data))
                analyse.(analyses, Ref(execution), analyses_path, mod_name)
            else
                failures = push_namedtuple!(failures, these_params)
            end
        end
        put!(results_channel, (failures, results))
    end
    @show fetch(task)
    return fetch(failures)
end
function execute_modifications_parallel_nosaving(example, modifications::Array{<:Dict}, analyses,
        data_path::String, analyses_path::String, parallel_batch_size, progress=false)
    pkeys = mods_to_pkeys(modifications)
    results_channel = RemoteChannel(() -> Channel{NamedTuple}(10 * nworkers()))
    batches = Iterators.partition(modifications, parallel_batch_size)
    n_batches = length(batches)
    task = @async @distributed for modifications_batch in collect(batches)
        failures = nothing
        for modification in modifications_batch
            mod_name, execution = execute_single_modification(example, modification)
            these_params = extract_params_tuple(modification, pkeys)
            if execution !== nothing #is success
                analyse.(analyses, Ref(execution), analyses_path, mod_name)
            else
                failures = push_namedtuple!(failures, these_params)
            end
        end
        put!(results_channel, failures)
    end
    n_remaining_batches = n_batches
    all_failures = nothing
    while n_remaining_batches > 0
        these_failures = take!(results_channel)
        all_failures = push_namedtuple!(all_failures, these_failures)
        n_remaining_batches -= 1
        if progress
            done_batches = n_batches - n_remaining_batches
            println("batches completed: $(done_batches) / $(n_batches)")
        end
    end
    @show fetch(task)
    return all_failures
end
function merge_ddb(::Nothing, tbl)
    return distribute(tbl, nworkers())
end
function merge_ddb(ddb::JuliaDB.DNDSparse, tbl::NDSparse)
    @error "Sparse unsupported"
    return merge(ddb, tbl)
end
function merge_ddb(ddb::JuliaDB.DIndexedTable, tbl::IndexedTable)
    return merge(ddb, tbl)
end

function push_namedtuple!(::Nothing, mods::NamedTuple{NAMES,TYPES}) where {NAMES,TYPES}
    arrd_TYPES = Tuple{[Array{T,1} for T in TYPES.parameters]...}
    NamedTuple{NAMES, arrd_TYPES}([[val] for val in values(mods)])
end
function push_namedtuple!(tup::NamedTuple, mods)
    for mod in pairs(mods)
        push!(tup[mod[1]], mod[2])
    end
end
push_namedtuple!(::Nothing, ::Nothing) = nothing
push_namedtuple!(nt::NamedTuple, ::Nothing) = nt
