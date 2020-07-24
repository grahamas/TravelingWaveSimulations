using DrWatson
include(projectdir("repl", "setup", "basic.jl"))
include(projectdir("repl", "setup", "plot.jl"))
include(joinpath(@__DIR__, "figure_contrast_monotonic_blocking_propagation.jl"))
include(joinpath(@__DIR__, "figure_connectivity.jl"))
include(joinpath(@__DIR__, "figure_nonlinearity.jl"))

DEFAULT_MODELS_RES=(800,800)
function figure_model!(scene::Scene, simulation::Simulation; prototype_name, 
                                                      modifications, 
                                                      save_dir=nothing,
                                                      scene_resolution=DEFAULT_MODELS_RES)
    
    connectivity_pane = connectivity_plot!(scene, simulation; save_dir=save_dir)
    nonlinearity_pane = nonlinearity_plot!(scene, simulation; save_dir=save_dir)
    #stimulus_pane = stimulus_plot!(scene, simulation; save_dir=save_dir)

    layout = GridLayout()
    layout[1,1] = connectivity_pane
    layout[2,1] = nonlinearity_pane
    return layout

end

function figure_model(simulation::Simulation; scene_resolution=DEFAULT_MODELS_RES,
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

    layout[1,1] = figure_model!(scene, simulation; save_dir=save_dir, kwargs...)

    if save
        Makie.save(joinpath(_save_dir, "full_model.png"), scene)
    end

    return layout

end

function figure_model(prototype::Function, modifications=Dict(); kwargs...)
    figure_model(prototype(; modifications...);
                 modifications=modifications,
                 kwargs...)
end

function figure_model(fpath::String, fixing_mods; kwargs...)
    prototype_name, sim_params = _fpath_params(fpath)
    prototype = get_prototype(prototype_name)
    figure_model(prototype, merge(sim_params, Dict(pairs(fixing_mods)...));prototype_name=prototype_name, kwargs...)
end


