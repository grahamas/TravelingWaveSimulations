
fcmb_monotonic_A_fpath = "/home/graham/data/dos_ring/Aee=40.0:12.0:250.0;Aei=20.0:16.0:250.0;Aie=15.0:16.0:250.0;Aii=1.0:20.0:250.0;blocking_θE=25.0;blocking_θI=25.0;firing_θE=6.0;firing_θI=7.0;step_reduction=nothing_2020-07-08T21:24:38.395_v1.0-208-gc225329_dirty"
fcmb_blocking_A_fpath = "/home/graham/data/dos_ring/Aee=40.0:12.0:250.0;Aei=20.0:16.0:250.0;Aie=15.0:16.0:250.0;Aii=1.0:20.0:250.0;blocking_θE=25.0;blocking_θI=10.0;firing_θE=6.0;firing_θI=7.0;step_reduction=nothing_2020-07-06T11:27:29.306_v1.0-205-g3f85b45_dirty"

fcmb_monotonic_S_fpath = "/home/graham/data/dos_ring/See=14.0:5.0:100.0;Sei=14.0:5.0:100.0;Sie=14.0:5.0:100.0;Sii=14.0:5.0:100.0;blocking_θE=25.0;blocking_θI=25.0;firing_θE=6.0;firing_θI=7.0;save_everystep=false_2020-07-09T08:35:00.708_v1.0-210-g00df749_dirty"
fcmb_blocking_S_fpath = "/home/graham/data/dos_ring/See=14.0:5.0:100.0;Sei=14.0:5.0:100.0;Sie=14.0:5.0:100.0;Sii=14.0:5.0:100.0;blocking_θE=25.0;blocking_θI=10.0;firing_θE=6.0;firing_θI=7.0;save_everystep=false_2020-07-09T08:04:38.231_v1.0-210-g00df749_dirty"

using DrWatson
include(projectdir("repl", "setup", "basic.jl"))
include(projectdir("repl", "setup", "plot.jl"))


function calc_binary_segmentation(arr)
    never = sum(arr .== 0)
    always = sum(arr .== 1)
    total = prod(size(arr))
    sometimes = total - (always + never)
    return (never = never / total,
            sometimes = sometimes / total,
            always = always / total)
end
function _fpath_names(fpath)
    path_components = splitpath(fpath)
    @assert path_components[end-2] == "data"
    prototype_name, sim_name = path_components[end-1:end]
    return (prototype_name, sim_name)
end
function _fpath_params(fpath)
    prototype_name, sim_name = _fpath_names(fpath)
    sim_params = TravelingWaveSimulations.parse_modifications_filename(sim_name)
    return prototype_name, sim_params
end

_namedaxisarray_names(naa::NamedAxisArray{names}) where names = names
_getaxis(naa::NamedAxisArray, dims) = naa.data.axes[[dim(naa, dims)...]]
function save_figure_contrast_monotonic_blocking((x_sym, y_sym)::Tuple{Symbol,Symbol}, 
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol,
                                     unique_id::String="")
    scene, _ = figure_contrast_monotonic_blocking((x_sym, y_sym),
                                               monotonic_fpath,
                                               blocking_fpath,
                                               property_sym)
    fname = "figure_contrast_monotonic_blocking_$(x_sym)_$(y_sym)_$(property_sym).png"
    mkpath(plotsdir(unique_id))
    Makie.save(plotsdir(unique_id,fname), scene)
end
function figure_contrast_monotonic_blocking((x_sym, y_sym)::Tuple{Symbol,Symbol}, 
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol;
                                     scene_resolution=(1200,1200))
    monotonic_prototype_name, monotonic_sim_name = _fpath_names(monotonic_fpath)
    blocking_prototype_name, blocking_sim_name = _fpath_names(blocking_fpath)
    @show monotonic_sim_name, monotonic_prototype_name
    @show blocking_sim_name, blocking_prototype_name
    @assert monotonic_prototype_name == blocking_prototype_name

    monotonic_prototype, blocking_prototype = get_prototype.((monotonic_prototype_name, blocking_prototype_name))


    (monotonic_x,monotonic_y,monotonic_data), (blocking_x,blocking_y,blocking_data) = map((monotonic_fpath, blocking_fpath)) do fpath
        A = TravelingWaveSimulations.load_ExecutionClassifications(AbstractArray, fpath)[property_sym]
        name_syms = _namedaxisarray_names(A)
        collapsed_syms = Tuple(setdiff(name_syms, (y_sym, x_sym)))
        @assert length(collapsed_syms) == length(name_syms) - 2
        reduced_data = Simulation73.avg_across_dims(A, collapsed_syms)
        x, y = _getaxis(A, (x_sym, y_sym)) .|> ax -> ax.keys
        if findfirst(name_syms .== y_sym) < findfirst(name_syms .== x_sym)
            return (x,y,reduced_data')
        else
            return (x,y,reduced_data)
        end
    end

    scene, layout = layoutscene(resolution=scene_resolution)
    
    title_size = 20
    monotonic_title = layout[1,1] = LText(scene, "Monotonic nonl.", textsize=20)
    blocking_title = layout[1,2] = LText(scene, "Blocking nonl.", textsize=20)

    monotonic_sweep_ax = LAxis(scene) 
    monotonic_heatmap = heatmap!(monotonic_sweep_ax, monotonic_x,monotonic_y,monotonic_data, colorrange=(0,1))
    blocking_sweep_ax = LAxis(scene) 
    blocking_heatmap = heatmap!(blocking_sweep_ax,blocking_x,blocking_y, blocking_data, colorrange=(0,1))
    tightlimits!.([monotonic_sweep_ax, blocking_sweep_ax])
    linkaxes!(monotonic_sweep_ax, blocking_sweep_ax)
    hideydecorations!(blocking_sweep_ax)
    monotonic_sweep_ax.ylabel = string(y_sym)
    monotonic_sweep_ax.xlabel = blocking_sweep_ax.xlabel = string(x_sym)
    sweep_sublayout  = GridLayout()
    sweep_sublayout[:h] = [monotonic_sweep_ax, blocking_sweep_ax]
    layout[2,1:2] = sweep_sublayout

    monotonic_histogram = StatsMakie.histogram(monotonic_data[:], closed=:right)
    blocking_histogram = StatsMakie.histogram(blocking_data[:], closed=:right)
    histogram_sublayout = GridLayout()

    monotonic_histogram_ax = LAxis(scene)
    plot!(monotonic_histogram_ax, monotonic_histogram)
    xlims!(monotonic_histogram_ax, 0, 1)
    blocking_histogram_ax = LAxis(scene)
    plot!(blocking_histogram_ax, blocking_histogram)
    xlims!(blocking_histogram_ax, 0, 1)
    linkaxes!(monotonic_histogram_ax, blocking_histogram_ax)
    tightlimits!.([monotonic_histogram_ax, blocking_histogram_ax])
    hideydecorations!(blocking_histogram_ax)
    monotonic_histogram_ax.ylabel = "count"
    monotonic_histogram_ax.xlabel = blocking_histogram_ax.xlabel = "proportion"

    histogram_sublayout[:h] = [monotonic_histogram_ax, blocking_histogram_ax]
    layout[3,3:4] = histogram_sublayout

    monotonic_title.tellwidth = false
    blocking_title.tellwidth = false
    
    return scene, layout
    
end

function save_figure_example_contrast_monotonic_blocking((x_sym, y_sym)::Tuple{Symbol,Symbol}, 
                                                         example_specs::Array{<:NamedTuple},
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol,
                                     unique_id::String="")
    scene, _ = figure_example_contrast_monotonic_blocking((x_sym, y_sym),
                                               example_specs,
                                               monotonic_fpath,
                                               blocking_fpath,
                                               property_sym)
    fname = "figure_examples_contrast_monotonic_blocking_$(x_sym)_$(y_sym)_$(property_sym).png"
    mkpath(plotsdir(unique_id))
    Makie.save(plotsdir(unique_id,fname), scene)
end
function figure_example_contrast_monotonic_blocking((x_sym, y_sym)::Tuple{Symbol,Symbol}, 
                                                    example_specs::Array{<:NamedTuple},
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol)
    monotonic_prototype_name, monotonic_spec = _fpath_params(monotonic_fpath)
    blocking_prototype_name, blocking_spec = _fpath_params(blocking_fpath)
    @assert monotonic_prototype_name == blocking_prototype_name

    @show monotonic_spec
    monotonic_prototype, blocking_prototype = get_prototype.((monotonic_prototype_name, blocking_prototype_name))

    scene_height = 450 * (2 + length(example_specs))
    scene_width = 450 * 4
    scene, layout = figure_contrast_monotonic_blocking((x_sym, y_sym), monotonic_fpath, blocking_fpath, property_sym; scene_resolution=(scene_width, scene_height))

    for spec in example_specs
        _, monotonic_example_exec = execute_single_modification(monotonic_prototype, merge(monotonic_spec, pairs(spec)))
        _, blocking_example_exec = execute_single_modification(blocking_prototype, merge(blocking_spec, pairs(spec)))
        
        examples_layout = GridLayout()
        examples_layout[2,1] = monotonic_ax = exec_heatmap!(scene, monotonic_example_exec)
        examples_layout[2,2] = blocking_ax = exec_heatmap!(scene, blocking_example_exec)
        examples_layout[1,:] = LText(scene, "$(spec)", tellwidth=false)
        
        layout[end+1, 1:4] = examples_layout
    end
    return scene, layout


end

