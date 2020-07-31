
function stimulus_plot!(scene, simulation::Simulation; save_dir)
    layout = GridLayout()
    layout[1,1] = ax = LAxis(scene,
            ylabel="space (μm)",
            xlabel="time (ms)",
            title="Stimulation")

    xs = coordinate_axes(Simulation73.reduced_space(simulation))[1] |> collect
    stim = simulation.model.stimulus.p1
    stim_action! = stim(simulation.space)

    time_windows = stim.time_windows
    end_time = time_windows[end][end] * 1.1
    
    n_t = 100
    ts = range(0, stop=end_time, length=n_t)
    stim_arr = zeros(Float64, length(xs), n_t) 

    for (i_t,t) in enumerate(ts)
        @views stim_action!(stim_arr[:,i_t], nothing, t)
    end
    
    htmp = heatmap!(ax, ts, xs, stim_arr')
    tightlimits!(ax)
    cbar = layout[1,2] = LColorbar(scene, htmp)
    cbar.width = 25
    return layout
end
