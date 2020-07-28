
#fcmb_monotonic_A_fpath = joinpath(homedir(), "data/ring_monotonic/report2/2020-07-26T15:53:11.486_v1.0-239-ga150dd3_dirty")
#fcmb_blocking_A_fpath = joinpath(homedir(), "data/ring_blocking/report2/2020-07-26T16:29:43.140_v1.0-239-ga150dd3_dirty")
#
#fcmb_monotonic_S_fpath = joinpath(homedir(), "data/ring_monotonic/report2/2020-07-26T17:16:51.735_v1.0-241-g8392d24")
#fcmb_blocking_S_fpath = joinpath(homedir(), "data/ring_blocking/report2/2020-07-26T17:17:24.569_v1.0-241-g8392d24")

fcmb_monotonic_A_fpath = get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_monotonic", "report2_A_sweep"))
fcmb_blocking_A_fpath = get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_blocking", "report2_A_sweep"))

fcmb_monotonic_S_fpath = get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_monotonic", "report2_S_sweep"))
fcmb_blocking_S_fpath = get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_blocking", "report2_S_sweep"))

@show fcmb_monotonic_S_fpath

using DrWatson
include(projectdir("repl", "setup", "basic.jl"))
include(projectdir("repl", "setup", "plot.jl"))


function calc_binary_segmentation(arr)
    never = sum(arr .== 0)
    always = sum(arr .== 1)
    total = prod(size(arr))
    sometimes = total - (always + never)
    return (none = never / total,
            some = sometimes / total,
            all = always / total)
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

l2(x,y) = sqrt(sum((x .- y) .^ 2))
using Interpolations
function diagonal_slice(x_axis, y_axis, data::Matrix, y_intercept, slope, dx=1)
    interp = LinearInterpolation((x_axis, y_axis), data')
    calc_y(x) = slope * x + y_intercept
    sample_points = [(x, calc_y(x)) for x in x_axis if y_axis[begin] <= calc_y(x) <= y_axis[end]]
    @show first(sample_points)
    sample_values = sample_points .|> (x) -> interp(x...)
    @show (sample_points[begin], sample_points[end])
    distance_along_line = l2.(Ref(sample_points[begin]), sample_points)
    return (distance_along_line, sample_points, sample_values)
end

function slice_2d_and_steepest_line_and_histogram!(scene::Scene, 
                                                   (x_sym, y_sym)::Tuple{Symbol,Symbol},
                                                   fpath::String,
                                                   property_sym::Symbol; 
                                                   title, titlesize=20, hide_y=false,
                                                   colorbar_width=nothing,
                                                   slice_y_intersect,
                                                   slice_slope)
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

    segmented = calc_binary_segmentation(data)
    segmented_ax = LAxis(scene)
    segment_names = string.([keys(segmented)...])
    barplot!(segmented_ax, [values(segmented)...])
    ylims!(segmented_ax, 0, 1)
    segmented_ax.xticks =( 1:3, segment_names)
    segmented_ax.xticklabelrotation = 0.0

    slice_x, slice_sample_points, slice_val = diagonal_slice(x, y, reduced_data, slice_y_intersect,
                                                            slice_slope)
    diagonal_slice_ax = LAxis(scene)
    plot!(diagonal_slice_ax, slice_x, slice_val)
    diagonal_slice_ax.xticks = ([slice_x[begin], slice_x[end-1]], 
                                string.([floor.(Ref(Int), slice_sample_points[begin]), 
                                         floor.(Int, slice_sample_points[end])]))
    tightlimits!(diagonal_slice_ax)
    ylims!(diagonal_slice_ax, 0, 1)


    if colorbar_width === nothing
        layout[2,1] = sweep_ax
    else
        sweep_layout = GridLayout()
        sweep_layout[:h] = [sweep_ax, LColorbar(scene, heatmap, width=colorbar_width)]
        layout[2,1] = sweep_layout
    end
    
    summary_layout = GridLayout()
    summary_layout[:v] = [segmented_ax, diagonal_slice_ax]
    layout[3,1] = summary_layout

    return layout

end


function figure_contrast_monotonic_blocking((x_sym, y_sym)::Tuple{Symbol,Symbol}, 
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol;
                                     scene_resolution=(800,1200),
                                     slice_y_intersect, slice_slope)
    
    scene, layout = layoutscene(resolution=scene_resolution)

    layout[1,1] = monotonic_sweep_ax = slice_2d_and_steepest_line_and_histogram!(scene, 
                                                            (x_sym, y_sym),
                                                            monotonic_fpath,
                                                            property_sym; 
                                                            slice_y_intersect=slice_y_intersect,
                                                            slice_slope=slice_slope,
                                                            title="Monotonic")
    layout[1,2] = blocking_sweep_ax = slice_2d_and_steepest_line_and_histogram!(scene, 
                                                            (x_sym, y_sym),
                                                            blocking_fpath,
                                                            property_sym; 
                                                            slice_y_intersect=slice_y_intersect,
                                                            slice_slope=slice_slope,
                                                            title="Blocking", hide_y=true)

    return scene, layout
    
end


function figure_contrast_monotonic_blocking_all((x_sym, y_sym)::Tuple{Symbol,Symbol}, (other_x_sym, other_y_sym)::Tuple{Symbol,Symbol},
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol;
                                     scene_resolution=(1200,1200),
                                     primary_slice_y_intersect, primary_slice_slope,
                                     secondary_slice_y_intersect, secondary_slice_slope)
    
    scene, layout = layoutscene(resolution=scene_resolution)


    layout[1:3,1] = monotonic_sweep_ax = slice_2d_and_steepest_line_and_histogram!(
                                                            scene, 
                                                            (x_sym, y_sym),
                                                            monotonic_fpath,
                                                            property_sym; 
                                                            slice_y_intersect = primary_slice_y_intersect,
                                                            slice_slope= primary_slice_slope,
                                                            title="Monotonic")
    layout[1:3,2] = blocking_sweep_ax = slice_2d_and_steepest_line_and_histogram!(
                                                            scene, 
                                                            (x_sym, y_sym),
                                                            blocking_fpath,
                                                            property_sym; 
                                                            slice_y_intersect = primary_slice_y_intersect,
                                                            slice_slope= primary_slice_slope,
                                                            title="Blocking", 
                                                            hide_y=true)
    layout[1:3,3] = other_monotonic_sweep_ax = slice_2d_and_steepest_line_and_histogram!(
                                                           scene, 
                                                           (other_x_sym, other_y_sym),
                                                           monotonic_fpath,
                                                           property_sym; 
                                                           slice_y_intersect = secondary_slice_y_intersect,
                                                           slice_slope= secondary_slice_slope,
                                                           title="Monotonic")
    layout[1:3,4] = other_blocking_sweep_ax = slice_2d_and_steepest_line_and_histogram!(
                                                          scene, 
                                                          (other_x_sym, other_y_sym),
                                                          blocking_fpath,
                                                          property_sym; 
                                                          slice_y_intersect = secondary_slice_y_intersect,
                                                          slice_slope= secondary_slice_slope,
                                                          title="Blocking",
                                                          hide_y=true,
                                                          colorbar_width=25)
    layout[end+1, 1] = LText(scene, "($(x_sym), $(y_sym))", tellwidth=false)
    layout[end+1, 2] = LText(scene, "($(x_sym), $(y_sym))", tellwidth=false)
    layout[end, 3] = LText(scene, "($(other_x_sym), $(other_y_sym))", tellwidth=false)
    layout[end, 4] = LText(scene, "($(other_x_sym), $(other_y_sym))", tellwidth=false)

    layout[end-1:end-2, 0] = LText(scene, "prop. sims", rotation=pi/2, tellheight=false)

    return scene, layout
    
end

function save_figure_example_contrast_monotonic_blocking((x_sym, y_sym)::Tuple{Symbol,Symbol}, 
                                                         example_specs::Array{<:NamedTuple},
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol,
                                     unique_id::String=""; kwargs...)
    scene, _ = figure_example_contrast_monotonic_blocking((x_sym, y_sym),
                                               example_specs,
                                               monotonic_fpath,
                                               blocking_fpath,
                                               property_sym; kwargs...)
    fname = "figure_examples_contrast_monotonic_blocking_$(x_sym)_$(y_sym)_$(property_sym).png"
    mkpath(plotsdir(unique_id))
    Makie.save(plotsdir(unique_id,fname), scene)
end
function save_figure_example_contrast_monotonic_blocking_all((x_sym, y_sym)::Tuple{Symbol,Symbol}, other_syms::Tuple{Symbol,Symbol}, 
                                                         example_specs::Array{<:NamedTuple},
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol,
                                     unique_id::String=""; kwargs...)
    scene, _ = figure_example_contrast_monotonic_blocking_all((x_sym, y_sym), other_syms,
                                               example_specs,
                                               monotonic_fpath,
                                               blocking_fpath,
                                               property_sym; kwargs...)
    fname = "figure_examples_contrast_monotonic_blocking_all_$(x_sym)_$(y_sym)_$(property_sym).png"
    mkpath(plotsdir(unique_id))
    Makie.save(plotsdir(unique_id,fname), scene)
end
#function figure_example_contrast_monotonic_blocking((x_sym, y_sym)::Tuple{Symbol,Symbol}, 
#                                                    example_specs::Array{<:NamedTuple},
#                                     monotonic_fpath::AbstractString, 
#                                     blocking_fpath::AbstractString,
#                                     property_sym::Symbol)
#    monotonic_prototype_name, monotonic_spec = _fpath_params(monotonic_fpath)
#    blocking_prototype_name, blocking_spec = _fpath_params(blocking_fpath)
#    @assert monotonic_prototype_name == blocking_prototype_name
#
#    @show monotonic_spec
#    monotonic_prototype, blocking_prototype = get_prototype.((monotonic_prototype_name, blocking_prototype_name))
#
#    scene_height = 450 * (2 + length(example_specs))
#    scene_width = 450 * 4
#    scene, layout = figure_contrast_monotonic_blocking((x_sym, y_sym), monotonic_fpath, blocking_fpath, property_sym; scene_resolution=(scene_width, scene_height))
#
#    count = 1
#    for spec in example_specs
#        _, monotonic_example_exec = execute_single_modification(monotonic_prototype, merge(monotonic_spec, pairs(spec)))
#        _, blocking_example_exec = execute_single_modification(blocking_prototype, merge(blocking_spec, pairs(spec)))
#        
#        examples_layout = GridLayout()
#        no_labels = count < length(example_spec)
#        @show "NO_LABELS: $no_labels"
#        examples_layout[2,1] = monotonic_ax = exec_heatmap!(scene, monotonic_example_exec; clims=(0.,0.5), no_labels=no_labels)
#        examples_layout[2,2] = blocking_ax = exec_heatmap!(scene, blocking_example_exec; clims=(0.,0.5), no_labels=true)
#        examples_layout[1,:] = LText(scene, "$(spec)", tellwidth=false)
#        
#        layout[end+1, 1:4] = examples_layout
#        count += 1
#    end
#    return scene, layout
#
#
#end


function figure_example_contrast_monotonic_blocking_all((x_sym, y_sym)::Tuple{Symbol,Symbol}, 
                                                        other_syms::Tuple{Symbol,Symbol},
                                                        example_specs::Array{<:NamedTuple},
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol; kwargs...)
    monotonic_prototype_name, monotonic_spec = _fpath_params(monotonic_fpath)
    blocking_prototype_name, blocking_spec = _fpath_params(blocking_fpath)
    @show blocking_spec

    @show monotonic_spec
    monotonic_prototype, blocking_prototype = get_prototype.((monotonic_prototype_name, blocking_prototype_name))

    scene_height = 450 * 3# (2 + length(example_specs))
    scene_width = 450 * 8
    scene, layout = figure_contrast_monotonic_blocking_all((x_sym, y_sym), other_syms, monotonic_fpath, blocking_fpath, property_sym; scene_resolution=(scene_width, scene_height), kwargs...)
    scene.attributes.attributes[:fontsize] = 20

    @show merge(blocking_spec, pairs(example_specs[2]))
    examples_layout = GridLayout()
    count = 0
    for spec in example_specs
        count += 1
        _, monotonic_example_exec = execute_single_modification(monotonic_prototype, merge(monotonic_spec, pairs(spec), Dict(:save_everystep => true)))
        @show ExecutionClassifications(monotonic_example_exec, n_traveling_frames_threshold=60).has_propagation
        _, blocking_example_exec = execute_single_modification(blocking_prototype, merge(blocking_spec, pairs(spec), Dict(:save_everystep => true)))
        @show ExecutionClassifications(blocking_example_exec, n_traveling_frames_threshold=60).has_propagation
        
        examples_layout[count,1] = monotonic_ax = exec_heatmap!(scene, monotonic_example_exec; clims=(0.0,0.5), no_labels=true)
        examples_layout[count,2] = blocking_ax = exec_heatmap!(scene, blocking_example_exec;
                                                           clims=(0.0,0.5), no_labels=true)
        examples_layout[count,3] = spec_text = LText(scene, 
                                                     join(["$k = $v" for (k,v) in pairs(spec)], "\n"), tellheight=false)
    end

    examples_layout[0,1] = LText(scene, "Monotonic", tellwidth=false)
    examples_layout[1,2] = LText(scene, "Blocking", tellwidth=false)
    examples_layout[end,0] = LText(scene, "space (μm)", rotation=pi/2, tellheight=false)
    examples_layout[end+1,2] = LText(scene, "time (ms)", tellwidth=false)
    
    ex_x, ex_y = size(examples_layout)
    @show size(examples_layout)
    next_col = size(layout)[2] + 1
    layout[1:4,next_col:next_col+3] = examples_layout
    
	label_a = layout[1, 2, TopLeft()] = LText(scene, "A", textsize = 35,
    	font = "Noto Sans Bold", halign = :right)
	label_b = layout[1, 4, TopLeft()] = LText(scene, "B", textsize = 35,
    	font = "Noto Sans Bold", halign = :right)
	label_c = layout[1, 6, TopLeft()] = LText(scene, "C", textsize = 35,
    	font = "Noto Sans Bold", halign = :right)

    return scene, layout



end

