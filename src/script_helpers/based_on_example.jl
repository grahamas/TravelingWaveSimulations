
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
        batch=DEFAULT_SAVE_BATCH_SIZE)

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
        mkpath(data_path)
    else
        data_path = nothing
    end

    example = get_example(example_name)
    if !(@isdefined nworkers) || nworkers() == 1
        @warn "Not parallelizing parameter sweep."
        results, failures = execute_modifications(example, modifications, analyses, data_path, analyses_path, no_save_raw)
        @show failures
    else
        @warn "Parallelizing parameter sweep."
        if no_save_raw
            results, failures = execute_modifications_parallel_nosaving(example, modifications, analyses, data_path, analyses_path, batch)
        else
            results, failures = execute_modifications_parallel_saving(example, modifications, analyses, data_path, analyses_path, batch)
        end
        @show failures
    end
end

function getkeys(d, keys)
    [d[key] for key in keys]
end

function extract_data_namedtuple(execution)
    soln = execution.solution
    u = soln.u
    t = soln.t
    x = coordinates(space(execution)) |> collect
    return (u=u, t=t, x=x)
end


function extract_params_tuple(execution, modification, pkeys)
    return NamedTuple{Tuple(pkeys)}(getkeys(modification, pkeys))
end

function execute_single_modification(example, modification)
    mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=";")
    if mod_name == ""
        mod_name = "no_mod"
    end
    @show mod_name
    simulation = example(; modification...)
    execution = execute(simulation)
    if execution.solution.retcode == :Unstable
        @warn "$mod_name unstable!"
        return (mod_name, nothing)
    end
    return (mod_name, execution)
end
function mods_to_pkeys(modifications)
    pkeys = keys(modifications[1]) |> collect
    disallowed_keys = [:algorithm, :u, :x, :t, :n, :n_points, :extent]
    pkeys = filter((x) -> !(x in disallowed_keys), pkeys)
end
function execute_modifications(example, modifications::Array{<:Dict}, analyses,
        data_path::String, analyses_path::String, no_save_raw)
    # FIXME
    pkeys = mods_to_pkeys(modifications)
    data_namedtuple = nothing
    for modification in modifications
        mod_name, execution = execute_single_modification(example, modification)
        execution === nothing && push!(failures, mod_name)
        if !no_save_raw
            these_data = expand_soln_and_modification(execution, modification, pkeys)
            data_namedtuple = push_namedtuple!(data_namedtuple, these_data)
        end
        analyse.(analyses, Ref(execution), analyses_path, mod_name)
    end
    if !no_save_raw        
        results_table = table(data_namedtuple, pkey=pkeys)
        JuliaDB.save(results_table, joinpath(data_path, "raw_data.csv"))
    else
        results_table = nothing
    end
    return results_table
end
using Nullables
Base.convert(::Type{T}, x::Nullable{T}) where T = x.value
Base.getindex(nt::NamedTuple, dx::Array{Symbol}) = getindex.(Ref(nt), dx)
function execute_modifications_parallel_saving(example, modifications::Array{<:Dict}, analyses,
        data_path::String, analyses_path::String, parallel_batch_size)
    pkeys = mods_to_pkeys(modifications)
    results_channel = RemoteChannel(() -> Channel{Tuple}(2 * nworkers()))
    @distributed for modifications_batch in Iterators.partition(modifications, parallel_batch_size)
        results = nothing
        failures = nothing
        for modification in modifications_batch
            mod_name, execution = execute_single_modification(example, modification)
            these_params = extract_params_tuple(modification)
            if execution !== nothing #is success
                these_data = extract_data_namedtuple(execution)
                results = push_results(results, merge(these_params, these_data))
                analyse.(analyses, Ref(execution), analyses_path, mod_name)
            else
                failures = push_namedtuple!(these_params, failures)
            end
        end
        put!(results_channel, (failures, results))
    end
    n = length(modifications)
    all_failures = nothing
    counter = 1
    while n > 0
        (these_failures, these_data) = take!(results_channel)
        all_failures = push_namedtuple(all_failures, these_failures)
        save(these_successes, joinpath(data_path,"$(counter).jdb"))
        counter += 1
        n -= 1
    end
    save(table(all_failures), joinpath(data_path, "failures.jdb"))
    return failures
end
function execute_modifications_parallel_nosaving(example, modifications::Array{<:Dict}, analyses,
        data_path::String, analyses_path::String, parallel_batch_size)
    pkeys = mods_to_pkeys(modifications)
    results_channel = RemoteChannel(() -> Channel{NamedTuple}(2 * nworkers()))
    batches = Iterators.partition(modifications, parallel_batch_size)
    n_batches = length(batches)
    @distributed for modifications_batch in Iterators.partition(modifications, parallel_batch_size)
        failures = nothing
        for modification in modifications_batch
            mod_name, execution = execute_single_modification(example, modification)
            these_params = extract_params_tuple(modification)
            if execution !== nothing #is success
                analyse.(analyses, Ref(execution), analyses_path, mod_name)
            else
                failures = push_namedtuple!(these_params, failures)
            end
        end
        put!(results_channel, failures)
    end
    n_remaining_batches = n
    all_failures = nothing
    while n_remaining_batches > 0
        these_failures = take!(results_channel)
        all_failures = push_namedtuple(all_failures, these_failures)
        n_remaining_batches -= 1
    end
    return failures
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
    return tup
end
