using Colors

function connectivity_plot!(scene::Scene, 
                            simulation::Simulation73.Simulation; 
                            xlimits=nothing, save_dir=nothing,
                            title="Connectivity kernels")
    layout = GridLayout()
    ax = LAxis(scene, 
               xlabel="distance (μm)", 
               ylabel="connectivity strength (a.u.)",
               title=title)
    
    pop_names = simulation.model.pop_names
    n_pops = length(pop_names)
    colors = distinguishable_colors(n_pops, 
                                    parse(Colorant, ax.attributes[:backgroundcolor][]), 
                                    dropseed=true)
    
    # Connectivity must be calculated before subsampling
    raw_space = simulation.space
    midpoint_idx = NeuralModels.fft_center_idx(raw_space)
    midpoint = coordinates(raw_space)[midpoint_idx]
    dists = differences(raw_space, midpoint) .|> x -> x[1]
    dists[CartesianIndex(1):midpoint_idx] .*= -1
    
    x_plot_idxs = if xlimits !== nothing
        findfirst(dists .>= xlimits[1]):findlast(dists .<= xlimits[2])
    else
        Colon()
    end
    onto_idx_idx = 1
    connectivity_array = Simulation73.array(simulation.model.connectivity)
    plots = map(pairs(IndexCartesian(), connectivity_array) |> collect) do (idx, conn)
        kern = NeuralModels.kernel(conn, raw_space)
        plot!(ax, dists[x_plot_idxs], kern[x_plot_idxs], color=colors[Tuple(idx)[onto_idx_idx]])
    end
    ontos = [Tuple(idx)[onto_idx_idx] for idx in CartesianIndices(connectivity_array)]
    legend_names = ["onto $pop" for pop in pop_names]
    leg = LLegend(scene, plots[1:n_pops], legend_names,
                  tellheight=false, tellwidth=false,
                  halign=:right, valign=:top, orientation=:vertical)
    tightlimits!(ax)
    layout[1,1] = ax 
    layout[1,1] = leg
    
    layout
end
