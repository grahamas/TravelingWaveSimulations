
export based_on_example

function based_on_example(; data_root::AbstractString=datadir(), no_save_raw::Bool=false,
        example_name::AbstractString=nothing,
        modifications::AbstractArray=[],
        plot_specs::AbstractArray=[])

    modifications, modifications_prefix = parse_modifications_argument(modifications)
    plot_specs, plot_specs_prefix = parse_plot_specs_argument(plot_specs)

    # Initialize saving paths
    if length(plot_specs) > 0
        plots_path = joinpath(plotsdir(), example_name, "$(modifications_prefix)$(plot_specs_prefix)_$(Dates.now())_$(current_commit())")
        mkpath(plots_path)
    else
        plots_path = ""
    end
    if !no_save_raw
        sim_output_path = joinpath(data_root, "sim", example_name, "$(modifications_prefix)$(Dates.now())_$(current_commit())")
        mkpath(sim_output_path)
    else
        sim_output_path = nothing
    end

    example = get_example(example_name)
    if length(modifications) > 2
        plot_lock = Threads.Mutex()
        @warn "Parallelizing with $(Threads.nthreads()) threads."
        Threads.@threads for modification in modifications
            simulation = example(; modification...)
            execution = execute(simulation)
            mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=";")
            if mod_name == ""
                mod_name = "no_mod"
            end
            if !no_save_raw
                execution_dict = @dict execution
                @warn("not currently saving raw")
                #@tagsave("$(joinpath(sim_output_path, mod_name)).bson", execution_dict, true)
            end
            Threads.lock(plot_lock)
            plot_and_save.(plot_specs, Ref(execution), plots_path, mod_name)
            Threads.unlock(plot_lock)
        end
    else
        for modification in modifications
            simulation = example(; modification...)
            execution = execute(simulation)
            mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=";")
            if mod_name == ""
                mod_name = "no_mod"
            end
            if !no_save_raw
                execution_dict = @dict execution
                @warn("not currently saving raw")
                #@tagsave("$(joinpath(sim_output_path, mod_name)).bson", execution_dict, true)
            end
            plot_and_save.(plot_specs, Ref(execution), plots_path, mod_name)
        end
    end
end
