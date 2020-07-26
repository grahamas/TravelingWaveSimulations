
#fcmb_monotonic_A_fpath = "/home/graham/data/dos_ring/Aee=40.0:12.0:250.0;Aei=20.0:16.0:250.0;Aie=15.0:16.0:250.0;Aii=1.0:20.0:250.0;blocking_θE=25.0;blocking_θI=25.0;firing_θE=6.0;firing_θI=7.0;step_reduction=nothing_2020-07-08T21:24:38.395_v1.0-208-gc225329_dirty"
#fcmb_blocking_A_fpath = "/home/graham/data/dos_ring/Aee=40.0:12.0:250.0;Aei=20.0:16.0:250.0;Aie=15.0:16.0:250.0;Aii=1.0:20.0:250.0;blocking_θE=25.0;blocking_θI=10.0;firing_θE=6.0;firing_θI=7.0;step_reduction=nothing_2020-07-06T11:27:29.306_v1.0-205-g3f85b45_dirty"
#
#fcmb_monotonic_S_fpath = "/home/graham/data/dos_ring/See=14.0:5.0:100.0;Sei=14.0:5.0:100.0;Sie=14.0:5.0:100.0;Sii=14.0:5.0:100.0;blocking_θE=25.0;blocking_θI=25.0;firing_θE=6.0;firing_θI=7.0;save_everystep=false_2020-07-09T08:35:00.708_v1.0-210-g00df749_dirty"
#fcmb_blocking_S_fpath = "/home/graham/data/dos_ring/See=14.0:5.0:100.0;Sei=14.0:5.0:100.0;Sie=14.0:5.0:100.0;Sii=14.0:5.0:100.0;blocking_θE=25.0;blocking_θI=10.0;firing_θE=6.0;firing_θI=7.0;save_everystep=false_2020-07-09T08:04:38.231_v1.0-210-g00df749_dirty"

fcmb_monotonic_A_fpath = joinpath(homedir(), "data/ring_monotonic/report2/2020-07-26T15:53:11.486_v1.0-239-ga150dd3_dirty")
fcmb_blocking_A_fpath = joinpath(homedir(), "data/ring_blocking/report2/2020-07-26T16:29:43.140_v1.0-239-ga150dd3_dirty")

fcmb_monotonic_S_fpath = joinpath(homedir(), "data/ring_monotonic/report2/2020-07-26T17:16:51.735_v1.0-241-g8392d24")
fcmb_blocking_S_fpath = joinpath(homedir(), "data/ring_blocking/report2/2020-07-26T17:17:24.569_v1.0-241-g8392d24")

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
function _fpath_params(fpath)
    sim_params = read_modifications_from_data_path(fpath)
    path_components = splitpath(fpath)
    data_dx = findfirst(path_components .== "data")
    prototype_name = path_components[data_dx + 1]
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

function slice_2d_and_steepest_line_and_histogram!(scene::Scene, 
                                                   (x_sym, y_sym)::Tuple{Symbol,Symbol},
                                                   fpath::String,
                                                   property_sym::Symbol; 
                                                   title, titlesize=20, hide_y=false,
                                                   colorbar_width=nothing)
    prototype_name, sim_params = _fpath_params(fpath)

    A = TravelingWaveSimulations.load_ExecutionClassifications(AbstractArray, fpath)[property_sym]
    name_syms = _namedaxisarray_names(A)
    collapsed_syms = Tuple(setdiff(name_syms, (y_sym, x_sym)))
    @assert length(collapsed_syms) == length(name_syms) - 2
    reduced_data = Simulation73.avg_across_dims(A, collapsed_syms)
    x, y = _getaxis(A, (x_sym, y_sym)) .|> ax -> ax.keys
    data = if findfirst(name_syms .== y_sym) < findfirst(name_syms .== x_sym)
        reduced_data'
    else
        reduced_data
    end

    layout = GridLayout()
    title_facet = layout[1,1] = LText(scene, title, textsize=titlesize, tellwidth=false)

    sweep_ax = LAxis(scene) 
    heatmap = heatmap!(sweep_ax, x,y,data, colorrange=(0,1))
    tightlimits!(sweep_ax)
    sweep_ax.xlabel = string(x_sym)
    if hide_y
        hideydecorations!(sweep_ax)
    else
        sweep_ax.ylabel = string(y_sym)
    end

    histogram = StatsMakie.histogram(data[:], closed=:right)
    histogram_ax = LAxis(scene)
    plot!(histogram_ax, histogram)
    xlims!(histogram_ax, 0, 1)
    histogram_ax.xlabel = "proportion"
    if hide_y
        #hideydecorations!(histogram_ax)
    else
        histogram_ax.ylabel = "count"
    end

    segmented = calc_binary_segmentation(data)
    segmented_ax = LAxis(scene)
    segment_names = string.([keys(segmented)...])
    barplot!(segmented_ax, [values(segmented)...])
    ylims!(segmented_ax, 0, 1)
    segmented_ax.xticks =( 1:3, segment_names)
    #segmented_ax.xticklabelrotation = (pi,)

    # TODO Make diagonal slice

    if colorbar_width === nothing
        layout[2,1] = sweep_ax
    else
        sweep_layout = GridLayout()
        sweep_layout[:h] = [sweep_ax, LColorbar(scene, heatmap, width=colorbar_width)]
        layout[2,1] = sweep_layout
    end
    
    # TODO put the diagonal slice and the histogram next to each other
    # in a nested layout
    summary_layout = GridLayout()
    summary_layout[:h] = [histogram_ax, segmented_ax]
    layout[3,1] = summary_layout

    return layout

end


function figure_contrast_monotonic_blocking((x_sym, y_sym)::Tuple{Symbol,Symbol}, 
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol;
                                     scene_resolution=(1200,1200))
    
    scene, layout = layoutscene(resolution=scene_resolution)

    layout[1,1] = monotonic_sweep_ax = slice_2d_and_steepest_line_and_histogram!(scene, (x_sym, y_sym),
                                                            monotonic_fpath,
                                                            property_sym; title="Monotonic")
    layout[1,2] = blocking_sweep_ax = slice_2d_and_steepest_line_and_histogram!(scene, (x_sym, y_sym),
                                                            blocking_fpath,
                                                            property_sym; title="Blocking", hide_y=true)

    return scene, layout
    
end


function figure_contrast_monotonic_blocking_all((x_sym, y_sym)::Tuple{Symbol,Symbol}, (other_x_sym, other_y_sym)::Tuple{Symbol,Symbol},
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol;
                                     scene_resolution=(1200,1200))
    
    scene, layout = layoutscene(resolution=scene_resolution)

    layout[1,1] = monotonic_sweep_ax = slice_2d_and_steepest_line_and_histogram!(
                                                            scene, 
                                                            (x_sym, y_sym),
                                                            monotonic_fpath,
                                                            property_sym; 
                                                            title="Monotonic")
    layout[1,2] = blocking_sweep_ax = slice_2d_and_steepest_line_and_histogram!(
                                                            scene, 
                                                            (x_sym, y_sym),
                                                            blocking_fpath,
                                                            property_sym; 
                                                            title="Blocking", 
                                                            hide_y=true)
    layout[1,3] = other_monotonic_sweep_ax = slice_2d_and_steepest_line_and_histogram!(
                                                           scene, 
                                                           (other_x_sym, other_y_sym),
                                                           monotonic_fpath,
                                                           property_sym; 
                                                           title="Monotonic")
    layout[1,4] = other_blocking_sweep_ax = slice_2d_and_steepest_line_and_histogram!(
                                                          scene, 
                                                          (other_x_sym, other_y_sym),
                                                          blocking_fpath,
                                                          property_sym; 
                                                          title="Blocking",
                                                          hide_y=true,
                                                          colorbar_width=25)

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
function save_figure_example_contrast_monotonic_blocking_all((x_sym, y_sym)::Tuple{Symbol,Symbol}, other_syms::Tuple{Symbol,Symbol}, 
                                                         example_specs::Array{<:NamedTuple},
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol,
                                     unique_id::String="")
    scene, _ = figure_example_contrast_monotonic_blocking_all((x_sym, y_sym), other_syms,
                                               example_specs,
                                               monotonic_fpath,
                                               blocking_fpath,
                                               property_sym)
    fname = "figure_examples_contrast_monotonic_blocking_all_$(x_sym)_$(y_sym)_$(property_sym).png"
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

function figure_example_contrast_monotonic_blocking_all((x_sym, y_sym)::Tuple{Symbol,Symbol}, 
                                                        other_syms::Tuple{Symbol,Symbol},
                                                        example_specs::Array{<:NamedTuple},
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol)
    monotonic_prototype_name, monotonic_spec = _fpath_params(monotonic_fpath)
    blocking_prototype_name, blocking_spec = _fpath_params(blocking_fpath)

    @show monotonic_spec
    monotonic_prototype, blocking_prototype = get_prototype.((monotonic_prototype_name, blocking_prototype_name))

    scene_height = 450 * (2 + length(example_specs))
    scene_width = 450 * 4
    scene, layout = figure_contrast_monotonic_blocking_all((x_sym, y_sym), other_syms, monotonic_fpath, blocking_fpath, property_sym; scene_resolution=(scene_width, scene_height))

    for spec in example_specs
        _, monotonic_example_exec = execute_single_modification(monotonic_prototype, merge(monotonic_spec, pairs(spec), Dict(:save_everystep => true)))
        @show ExecutionClassifications(monotonic_example_exec).has_propagation
        _, blocking_example_exec = execute_single_modification(blocking_prototype, merge(blocking_spec, pairs(spec), Dict(:save_everystep => true)))
        @show ExecutionClassifications(blocking_example_exec).has_propagation
        
        examples_layout = GridLayout()
        examples_layout[2,1] = monotonic_ax = exec_heatmap!(scene, monotonic_example_exec)
        examples_layout[2,2] = blocking_ax = exec_heatmap!(scene, blocking_example_exec)
        examples_layout[1,:] = LText(scene, "$(spec)", tellwidth=false)
        
        layout[end+1, 1:4] = examples_layout
    end
    return scene, layout


end

