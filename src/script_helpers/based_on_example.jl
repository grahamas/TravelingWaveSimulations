
export based_on_example

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
        analyses::AbstractArray=[])

    modifications, modifications_prefix = parse_modifications_argument(modifications)
    analyses, analyses_prefix = parse_analyses_argument(analyses)

    # Initialize saving paths
    if length(analyses) > 0
        analyses_path = joinpath(plotsdir(), example_name, "$(modifications_prefix)$(analyses_prefix)_$(Dates.now())_$(gitdescribe())")
        mkpath(analyses_path)
    else
        analyses_path = ""
    end
    if !no_save_raw
        sim_output_path = joinpath(data_root, "sim", example_name, "$(modifications_prefix)$(Dates.now())_$(gitdescribe())")
        mkpath(sim_output_path)
    else
        sim_output_path = nothing
    end

    example = get_example(example_name)
        for modification in modifications
            simulation = example(; modification...)
            execution = execute(simulation)
            mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=";")
            if execution.solution.retcode == :Unstable
                @warn "$mod_name unstable!"
                continue
            end
            if mod_name == ""
                mod_name = "no_mod"
            end
            if !no_save_raw
                execution_dict = @dict execution
                @warn("not currently saving raw")
                #@tagsave("$(joinpath(sim_output_path, mod_name)).bson", execution_dict, true)
            end
            analyse.(analyses, Ref(execution), analyses_path, mod_name)
        end
#    end
end
