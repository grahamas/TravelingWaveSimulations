
export based_on_example

DEFAULT_SAVE_BATCH_SIZE = 10000
"""
    based_on_example(; data_root=datadir(),
                       no_save_raw=false,
                       example_name=nothing,
                       modifications=[],
                       analyses=[])

Run simulation based on example named `example_name` described in the `examples` subdirectory of `src`.

`modifications` is an array of *strings* modifying the example's built-in parameters. These strings can 1) name a modifications file in the `modifications` subdirectory of `scripts`, 2) modify a keyword value using the syntax `name=value` or 3) describe a range of values for a keyword, using colon syntax (`name=start:step:end`). Valid keyword targets are *any* keywords in the signature or body of the example, so long as that keyword is *unique*.

`analyses` is an array of *strings* naming files (sans `.jl`) found in the `analyses` subdirectory of `scripts`.

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
    if !(@isdefined nworkers)
        @warn "Not parallelizing parameter sweep."
        execute_modifications(example, modifications, analyses, data_path, analyses_path, no_save_raw)
    else
        @warn "Parallelizing parameter sweep."
        execute_modifications_parallel(example, modifications, analyses, data_path, analyses_path, no_save_raw, batch)
    end
end

function getkeys(d, keys)
    [d[key] for key in keys]
end

function expand_soln_and_modification(execution, modification, pkeys)
   # u = soln.u[:]
   # u1 = soln(0.0)
   # dims = size(u1)
   # n_pops = dims[end]
   # t = repeat(soln.t, inner=prod(dims))
   # x = repeat(coordinates(space(execution))[:], inner=n_pops, outer=length(soln.t))
   # mod = NamedTuple{Tuple(pkeys)}(getkeys(modification, pkeys))
   # repeated_modification = NamedTuple{keys(mod)}([repeat([mod_val], inner=(prod(dims) * length(soln.t))) for mod_val in values(mod)])
   # return (u=u, t=t, x=x, repeated_modification...)
   soln = execution.solution
   u = soln.u
   t = soln.t
   x = coordinates(space(execution)) |> collect
   mod = NamedTuple{Tuple(pkeys)}(getkeys(modification, pkeys))
   return [mod, (u=u, t=t, x=x)]
end

function execute_modifications(example, modifications::Array{<:Dict}, analyses,
        data_path::String, analyses_path::String, no_save_raw)
    pkeys = keys(modifications[1]) |> collect
    disallowed_keys = [:algorithm, :u, :x, :t, :n, :n_points, :extent]
    pkeys = filter((x) -> !(x in disallowed_keys), pkeys)
    results = nothing; mods = nothing
    for modification in modifications
        mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=";")
        if mod_name == ""
            mod_name = "no_mod"
        end
        @show mod_name
        simulation = example(; modification...)
        execution = execute(simulation, )
        if execution.solution.retcode == :Unstable
            @warn "$mod_name unstable!"
            continue
        end
        if !no_save_raw
            these_mods, these_results = expand_soln_and_modification(execution, modification, 
                                                                            pkeys)
            results = push_results(results, these_results)
            mods = push_results(mods, these_mods)
        end
        analyse.(analyses, Ref(execution), analyses_path, mod_name)
    end
    if !no_save_raw
        JuliaDB.save(ndsparse(mods, results), joinpath(data_path, "raw_data.csv"))
    end
end
using Nullables
Base.convert(::Type{T}, x::Nullable{T}) where T = x.value
function execute_modifications_parallel(example, modifications::Array{<:Dict}, analyses,
        data_path::String, analyses_path::String, no_save_raw, parallel_batch_size)
    pkeys = keys(modifications[1]) |> collect
    disallowed_keys = [:algorithm, :u, :x, :t, :n, :n_points, :extent]
    pkeys = filter((x) -> !(x in disallowed_keys), pkeys)
    @show pkeys
    results_channel = RemoteChannel(() -> Channel{Array}(2 * nworkers()))
    @distributed for modification in modifications
        mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=";")
        if mod_name == ""
            mod_name = "no_mod"
        end
        @show mod_name
        simulation = example(; modification...)
        execution = execute(simulation)
        if execution.solution.retcode == :Unstable
            @warn "$mod_name unstable!"
            continue
        end
        if !no_save_raw
            put!(results_channel, expand_soln_and_modification(execution, modification, pkeys))
        end
        analyse.(analyses, Ref(execution), analyses_path, mod_name)
    end
    n = length(modifications)
    results_batch = nothing
    mods_batch = nothing
    ddb = nothing
    all_results = nothing
    @show n
    @show parallel_batch_size
    while n > 0
        these_mods, these_results = take!(results_channel)
        #results_batch = push_results(results_batch, these_results)
        #mods_batch = push_results(mods_batch, these_mods)
        all_results = push_results(all_results, NamedTuple{(pkeys...,:u,:t,:x),typeof((these_mods...,these_results...))}((these_mods..., these_results...)))
        if ((length(all_results[1]) >= parallel_batch_size) || (n == 1))
            println("writing! ($n)")
            #ddb = join_ddb(ddb, table(all_results, pkey=pkeys))# table((mods_batch..., results_batch...), pkey=pkeys))
            #@show ddb
            JuliaDB.save(table(all_results, pkey=pkeys), joinpath(data_path, "$n.jdb"))
            results_batch = nothing
            mods_batch = nothing
            all_results = nothing
        end
        n -= 1
    end
    #JuliaDB.save(ddb, data_path)
end
function join_ddb(::Nothing, tbl)
    #return distribute(tbl, 1)
    return tbl
end
function join_ddb(ddb::JuliaDB.DNDSparse, tbl::JuliaDB.NDSparse)
    return join(ddb, tbl)
end


function push_results(::Nothing, mods::NamedTuple{NAMES,TYPES}) where {NAMES,TYPES}
    arrd_TYPES = Tuple{[Array{T,1} for T in TYPES.parameters]...}
    NamedTuple{NAMES, arrd_TYPES}([[val] for val in values(mods)])
end
function push_results(tup::NamedTuple, mods)
    for mod in pairs(mods)
        push!(tup[mod[1]], mod[2])
    end
    return tup
end
