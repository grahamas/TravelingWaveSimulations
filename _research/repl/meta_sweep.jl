
using DrWatson
include(projectdir("repl", "setup", "basic.jl"))
include(projectdir("repl", "setup", "plot.jl"))
include(projectdir("repl", "phase_space_analysis.jl"))

function get_subspace_coordinates(na::NamedAxisArray{NAMES}, subspace_dims) where NAMES
    subspace_axes = axes_keys(na)[NAMES .∈ subspace_dims]
    map(NamedTuple{NAMES}, product(subspace_axes...))
end
function meta_sweep_from_large_run(meta_params,
                                   phase_space_axis_names,
                                   fpath;
                                   property_sym=:has_propagation)
    whole_run_result = TravelingWaveSimulations.load_ExecutionClassifications(
                                                     AbstractArray, fpath
                                                    )[property_sym]
    sigmoid_fit_results = NamedAxisArray{meta_params}(map(get_subspace_coordinates(whole_run_result, meta_params)) do meta_param_nt
        param_sweep_results = whole_run_result[meta_param_nt...] # get all results for given meta_param
        phase_space_sigmoid_fit(param_sweep_results, phase_space_axis_names)
    end, get_subspace_axes(whole_run_result, meta_params))
end

function heatmap_meta_sweep_from_large_run!(scene, meta_params,
                                           phase_space_axis_names,
                                           args...; kwargs...)
    sigmoid_fit_results = meta_sweep_from_large_run(meta_params, phase_space_axis_names, args...; kwargs...)

    layout = GridLayout()
    change_ax = layout[1,1] = LAxis(scene)
    slope_ax = layout[1,2] = LAxis(scene)
    threshold_ax = layout[2,1] = LAxis(scene)
    error_ax = layout[2,2] = LAxis(scene)

    
    change_heatmap = heatmap!(change_ax, map(s -> s.change, sigmoid_fit_results),
                              colorrange=(0,1))
    slope_heatmap = heatmap!(slope_ax, map(s -> s.slope, sigmoid_fit_results),
                              colorrange=(0,1))
    threshold_heatmap = heatmap!(threshold_ax, map(s -> s.threshold, sigmoid_fit_results),
                              colorrange=(0,1))
    error_heatmap = heatmap!(error_ax, map(s -> s.error, sigmoid_fit_results),
                              colorrange=(0,1))

    return (scene, layout)
end


function save_heatmap_meta_sweep_from_large_run(meta_params, phase_space_axis_names,
                                           args...; unique_id="", kwargs...)
    fname = "metasweep_META_$(join(meta_params, '_'))_PS_$(join(phase_space_axis_names, '_')).png"
    dir = plotsdir(unique_id)
    mkpath(dir)

    scene, _ = heatmap_meta_sweep_from_large_run(meta_params,
                    phase_space_axis_names, args...; kwargs...)

    path = plotsdir(unique_id, fname)
    @info "saving $(path)"
    Makie.save(path, scene)
end