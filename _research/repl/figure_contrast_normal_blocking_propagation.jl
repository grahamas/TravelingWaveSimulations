include(projectdir("repl", "setup", "basic.jl"))
include(projectdir("repl", "setup", "plot.jl"))

function single_2D_slice_heatmap(scene, A::NamedAxisArray{name_syms}, 
                                 y_axis_sym::Symbol,
                                 x_axis_sym::Symbol;
                                 show_colorbar::Bool=true) where {name_syms}
    collapsed_syms = Tuple(setdiff(name_syms, (y_axis_sym, x_axis_sym)))
    mean_values = Simulation73.avg_across_dims(A, collapsed_syms)
    ax = LAxis(scene)
    heatmap = heatmap!(ax, mean_values, colorrange=(0,1), color=:magma)
    tightlimits!(ax)
    if show_colorbar
        layout = GridLayout()
        layout[1,1] = ax
        layout[1,2] = cbar = LColorbar(scene, heatmap)
        cbar.width = 25
        return layout
    else
        return ax
    end
end
_fpath_names(fpath) = splitpath(fpath)[end-1:end]
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
        TravelingWaveSimulations.load_ExecutionClassifications(AbstractArray, fpath)[property_sym]
    end

    scene, layout = layoutscene()
    
    blocking_sweep_panel = layout[1,1] = single_2D_slice_heatmap(scene, blocking_data, y_sym, x_sym)
    monotonic_sweep_panel = layout[1,2] = single_2D_slice_heatmap(scene, monotonic_data, y_sym, x_sym)
    
end

