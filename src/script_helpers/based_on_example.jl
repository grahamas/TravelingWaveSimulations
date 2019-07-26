
export based_on_example

function init_path(args...)
    path = joinpath(args...) |> make_path_windows_safe
    mkpath(path)
    return path
end

function based_on_example(ARGS)
    arg_settings = ArgParseSettings(; autofix_names = true)
    @add_arg_table arg_settings begin
        "--data-root", "-d"
            help = "Location in which to store output data"
            default = datadir()
        "--example-name"
            help = "Name of example defined examples.jl"
        "--modifications", "--mod"
            nargs = '*'
            help = "Name of file specifying dict of modifications"
        "--plot-specs", "--plot"
            nargs = '*'
            help = "Name of file specifying plots"
        "--no-save-raw"
            help = "Don't save raw simulation"
            action = :store_true
    end
    kwargs = parse_args(ARGS, arg_settings; as_symbols=true)
    based_on_example(; kwargs...)
end



function based_on_example(; data_root::AbstractString=datadir(), no_save_raw::Bool=false,
        example_name::AbstractString=nothing,
        modifications::AbstractArray=[],
        plot_specs::AbstractArray=[])

    modifications, modifications_prefix = parse_modifications_argument(modifications)
    plot_specs, plot_specs_prefix = parse_plot_specs_argument(plot_specs)

    # Initialize saving paths
    if length(plot_specs) > 0
        plots_path = init_path(plotsdir(), example_name, "$(modifications_prefix)$(plot_specs_prefix)_$(Dates.now())_$(current_commit())")
    else
        plots_path = ""
    end
    if !no_save_raw
        sim_output_path = init_path(data_root, "sim", example_name, "$(modifications_prefix)$(Dates.now())_$(current_commit())")
    else
        sim_output_path = nothing
    end

    example = get_example(example_name)
#    if length(modifications) > 2
#        plot_lock = Threads.Mutex()
#        @warn "Parallelizing with $(Threads.nthreads()) threads."
#        Threads.@threads for modification in modifications
#            simulation = example(; modification...)
#            execution = execute(simulation)
#            mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=MOD_SEP)
#            if mod_name == ""
#                mod_name = "no_mod"
#            end
#            if !no_save_raw
#                execution_dict = @dict execution
#                @warn("not currently saving raw")
#                #@tagsave("$(joinpath(sim_output_path, mod_name)).bson", execution_dict, true)
#            end
#            Threads.lock(plot_lock)
#            plot_and_save.(plot_specs, Ref(execution), plots_path, mod_name)
#            Threads.unlock(plot_lock)
#        end
#    else
        for modification in modifications
            simulation = example(; modification...)
            execution = execute(simulation)
            mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=MOD_SEP)
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
#    end
end
