
export based_on_example

function based_on_example(; data_root::AbstractString=datadir(), no_save_raw::Bool=false,
        example_name::AbstractString=nothing,
        modifications::AbstractArray=[],
        analyses::AbstractArray=[])

    modifications, modifications_prefix = parse_modifications_argument(modifications)
    analyses, analyses_prefix = parse_analyses_argument(analyses)

    # Initialize saving paths
    if length(analyses) > 0
        analyses_path = joinpath(plotsdir(), example_name, "$(modifications_prefix)$(analyses_prefix)_$(Dates.now())_$(current_commit())")
        mkpath(analyses_path)
    else
        analyses_path = ""
    end
    if !no_save_raw
        sim_output_path = joinpath(data_root, "sim", example_name, "$(modifications_prefix)$(Dates.now())_$(current_commit())")
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
