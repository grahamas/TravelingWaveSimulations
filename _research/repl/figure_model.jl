using DrWatson
include(projectdir("repl", "setup", "basic.jl"))
include(projectdir("repl", "setup", "plot.jl"))
include(joinpath(@__DIR__, "figure_contrast_monotonic_blocking_propagation.jl"))
include(joinpath(@__DIR__, "figure_connectivity.jl"))
include(joinpath(@__DIR__, "figure_nonlinearity.jl"))
include(joinpath(@__DIR__, "figure_stimulus.jl"))

DEFAULT_MODELS_RES=(1200,1200)
function figure_models!(scene::Scene, monotonic_simulation::Simulation,
                                     blocking_simulation::Simulation; 
                                     monotonic_prototype_name,
                                     blocking_prototype_name,
                                     monotonic_modifications, 
                                     blocking_modifications,
                                     save_dir=nothing,
                                     scene_resolution=DEFAULT_MODELS_RES)
    
    connectivity_pane = connectivity_plot!(scene, monotonic_simulation; 
            title = "Connectivity kernels (example)", save_dir=save_dir)
    monotonic_nonlinearity_pane = nonlinearity_plot!(scene, monotonic_simulation; 
            title = "Monotonic nonlinearity", save_dir=save_dir)
    blocking_nonlinearity_pane = nonlinearity_plot!(scene, blocking_simulation; 
            title = "Blocking nonlinearity", save_dir = save_dir)
    stimulus_pane = stimulus_plot!(scene, monotonic_simulation; save_dir=save_dir)

    layout = GridLayout()
    layout[1,1] = connectivity_pane
    layout[1,2] = monotonic_nonlinearity_pane
    layout[2,2] = blocking_nonlinearity_pane
    layout[2,1] = stimulus_pane
	label_a = layout[1, 1, TopLeft()] = LText(scene, "A", textsize = 35,
    	font = "Noto Sans Bold", halign = :right)
	label_b = layout[1, 2, TopLeft()] = LText(scene, "B", textsize = 35,
    	font = "Noto Sans Bold", halign = :right)
	label_c = layout[2, 1, TopLeft()] = LText(scene, "C", textsize = 35,
    	font = "Noto Sans Bold", halign = :right)
	label_d = layout[2, 2, TopLeft()] = LText(scene, "D", textsize = 35,
    	font = "Noto Sans Bold", halign = :right)
    return layout

end

function figure_models(monotonic_simulation::Simulation,
                       blocking_simulation::Simulation; 
                       scene_resolution=DEFAULT_MODELS_RES,
                                                      save=false,
                                                      unique_id="",
                                                      kwargs...) 
    scene, layout = layoutscene(resolution=scene_resolution)
    save_dir = if save
        _save_dir = plotsdir("figure_model",unique_id)
        mkpath(_save_dir)
        _save_dir
    else
        nothing
    end

    layout[1,1] = figure_models!(scene, monotonic_simulation, blocking_simulation; save_dir=save_dir, kwargs...)

    if save
        Makie.save(joinpath(_save_dir, "full_model.png"), scene)
    end

    return layout

end

function figure_models(monotonic_prototype::Function, monotonic_modifications,
                      blocking_prototype::Function, blocking_modifications; kwargs...)
    figure_models(monotonic_prototype(; monotonic_modifications...),
                 blocking_prototype(; blocking_modifications...);
                 monotonic_modifications=monotonic_modifications,
                 blocking_modifications=blocking_modifications,
                 kwargs...)
end

function figure_models(mono_fpath::String, blk_fpath, fixing_mods; kwargs...)
    mono_prototype_name, mono_sim_params = _fpath_params(mono_fpath)
    mono_prototype = get_prototype(mono_prototype_name)
    blk_prototype_name, blk_sim_params = _fpath_params(blk_fpath)
    blk_prototype = get_prototype(blk_prototype_name)
    figure_models(mono_prototype, 
        merge(mono_sim_params, Dict(pairs(fixing_mods)...)), 
        blk_prototype,
        merge(blk_sim_params, Dict(pairs(fixing_mods)...)); 
        monotonic_prototype_name=mono_prototype_name, 
        blocking_prototype_name=blk_prototype_name,
        kwargs...)
end


