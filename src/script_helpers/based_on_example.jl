
export based_on_example

MIN_FOR_PARALLEL = 100
PARALLEL_BATCH_SIZE = 100

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
        batch=10)

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
    if length(modifications) < (batch * nworkers())
        @warn "Not parallelizing parameter sweep."
        execute_modifications(modifications, analyses, data_path, analyses_path, no_save_raw)
    else
        @warn "Parallelizing parameter sweep."
        execute_modifications_parallel(modifications, analyses, data_path, analyses_path, no_save_raw, batch)
    end
end

function execute_modifications(modifications::Array{<:Dict}, analyses,
        data_path::String, analyses_path::String, no_save_raw)
    for modification in modifications
        mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=";")
        @show mod_name
        simulation = example(; modification...)
        execution = execute(simulation)
        if execution.solution.retcode == :Unstable
            @warn "$mod_name unstable!"
            continue
        end
        if mod_name == ""
            mod_name = "no_mod"
        end
        if !no_save_raw
            #save_csv("$(joinpath(data_path, mod_name)).csv", execution, modification)
            save(data_path, execution, modification)
        end
        analyse.(analyses, Ref(execution), analyses_path, mod_name)
    end
end

function execute_modifications_parallel(modifications::Array{<:Dict}, analyses,
        data_path::String, analyses_path::String, no_save_raw, parallel_batch_size)
    results = nothing
    @distributed for modifications_subset in collect(IterTools.partition(modifications, parallel_batch_size))
        for modification in modifications_subset
            mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=";")
            simulation = example(; modification...)
            execution = execute(simulation)
            if execution.solution.retcode == :Unstable
                @warn "$mod_name unstable!"
                continue
            end
            if mod_name == ""
                mod_name = "no_mod"
            end
            if !no_save_raw
                soln = execution.solution
                push!(results, (u=soln.u,t=soln.t,x=coordinates(space(execution)),modification...))
            end
            analyse.(analyses, Ref(execution), analyses_path, mod_name)
        end
        open(joinpath(data_path, "$(myid()).csv"), "w") do f
            table(results, pkey=keys(modifications_subset[1])) |> CSV.write(f)
        end
    end
end
push_results(::Nothing, mods...) = NamedTuple{keys(mods)}([[val] for val in values(mods)])
function push_results(tup::NamedTuple, mods...)
    for mod in mods
        push!(tup[mod[1]], mod[2])
    end
end
