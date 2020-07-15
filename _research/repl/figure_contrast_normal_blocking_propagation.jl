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
_fpath_names(fpath) = splitpath(fpath)[end-1:end]
_namedaxisarray_names(naa::NamedAxisArray{names}) where names = names
function figure_contrast_monotonic_blocking((y_sym, x_sym)::Tuple{Symbol,Symbol}, 
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol)
    monotonic_prototype_name, monotonic_sim_name = _fpath_names(monotonic_fpath)
    blocking_prototype_name, blocking_sim_name = _fpath_names(blocking_fpath)
    @show monotonic_sim_name, monotonic_prototype_name
    @show blocking_sim_name, blocking_prototype_name
    @assert monotonic_prototype_name == blocking_prototype_name

    monotonic_prototype, blocking_prototype = get_prototype.((monotonic_prototype_name, blocking_prototype_name))


    monotonic_data, blocking_data = map((monotonic_fpath, blocking_fpath)) do fpath
        A = TravelingWaveSimulations.load_ExecutionClassifications(AbstractArray, fpath)[property_sym]
        name_syms = _namedaxisarray_names(A)
        collapsed_syms = Tuple(setdiff(name_syms, (y_sym, x_sym)))
        @assert length(collapsed_syms) == length(name_syms) - 2
        Simulation73.avg_across_dims(A, collapsed_syms)
    end

    scene, layout = layoutscene(resolution=(1200,1200))
    
    title_size = 20
    monotonic_title = layout[1,1] = LText(scene, "Monotonic nonl.", textsize=20)
    blocking_title = layout[1,2] = LText(scene, "Blocking nonl.", textsize=20)

    monotonic_sweep_ax = LAxis(scene) 
    monotonic_heatmap = heatmap!(monotonic_sweep_ax, monotonic_data, colorrange=(0,1))
    blocking_sweep_ax = LAxis(scene) 
    blocking_heatmap = heatmap!(blocking_sweep_ax, blocking_data, colorrange=(0,1))
    tightlimits!.([monotonic_sweep_ax, blocking_sweep_ax])
    linkaxes!(monotonic_sweep_ax, blocking_sweep_ax)
    hideydecorations!(blocking_sweep_ax)
    monotonic_sweep_ax.ylabel = string(y_sym)
    monotonic_sweep_ax.xlabel = blocking_sweep_ax.xlabel = string(x_sym)
    sweep_sublayout  = GridLayout()
    sweep_sublayout[:h] = [monotonic_sweep_ax, blocking_sweep_ax]
    layout[2,1:2] = sweep_sublayout

    monotonic_slice_ax = LAxis(scene)
    monotonic_slice = lines!(monotonic_slice_ax, monotonic_data[:,end÷2])
    blocking_slice_ax = LAxis(scene)
    blocking_slice = lines!(blocking_slice_ax, blocking_data[:,end÷2])
    linkaxes!(monotonic_slice_ax, blocking_slice_ax)
    hideydecorations!(blocking_slice_ax)
    monotonic_slice_ax.ylabel = "proportion"
    monotonic_slice_ax.xlabel = monotonic_slice_ax.xlabel = string(x_sym)
    slice_sublayout = GridLayout()
    slice_sublayout[:h] = [monotonic_slice_ax, blocking_slice_ax]
    layout[3,1:2] = slice_sublayout


    monotonic_histogram = StatsMakie.histogram(monotonic_data[:])
    blocking_histogram = StatsMakie.histogram(blocking_data[:])
    histogram_sublayout = GridLayout()

    monotonic_histogram_ax = LAxis(scene)
    plot!(monotonic_histogram_ax, monotonic_histogram)
    #ylims!(monotonic_histogram_ax, 0, 1)
    blocking_histogram_ax = LAxis(scene)
    plot!(blocking_histogram_ax, blocking_histogram)
    #ylims!(blocking_histogram_ax, 0, 1)
    linkaxes!(monotonic_histogram_ax, blocking_histogram_ax)
    hideydecorations!(blocking_histogram_ax)
    monotonic_histogram_ax.ylabel = "count"
    monotonic_histogram_ax.xlabel = blocking_histogram_ax.xlabel = "proportion"

    histogram_sublayout[:h] = [monotonic_histogram_ax, blocking_histogram_ax]
    layout[4,1:2] = histogram_sublayout

    monotonic_title.tellwidth = false
    blocking_title.tellwidth = false
    
    return scene
    
end

