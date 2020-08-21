
function heatmap_sweep_with_target(sweep::AbstractArray,
			target_mods_nt::NamedTuple{mod_names}, 
			prototype_name,
            sim_name;
            fixed_mods=Dict(),
            title,
			plot_color=:magma,
    		plot_side_size = 350 * (length(mod_names) - 1)
		) where {mod_names}
    mod_values = keys.(axes(sweep))
    mod_names_str = [string(name) for name in mod_names]
    line_prototype = get_prototype(prototype_name)
    all_dims = 1:length(mod_names)
    slices_2d = IterTools.subsets(all_dims, Val{2}())
    plot_size = (plot_side_size, plot_side_size)
    scene, layout = layoutscene(resolution=plot_size)

    heatmaps = map(slices_2d) do (x,y)
        (x,y) = x < y ? (x,y) : (y,x)
        my = mod_values[y] |> collect
        mx = mod_values[x] |> collect
        collapsed_dims = Tuple(setdiff(all_dims, (x,y)))
        sweep_2d_mean = Simulation73.avg_across_dims(sweep, collapsed_dims)
        
        @assert size(sweep_2d_mean) == length.((mod_values[x], mod_values[y]))
        
        layout[x,y] = ax = LAxis(scene); 
        tightlimits!(ax)
        
        heatmap = Makie.heatmap!(ax, my, mx, sweep_2d_mean', colorrange=(0,1))
        
        Makie.scatter!(ax, [target_mods_nt[y]], [target_mods_nt[x]], color=:red, markersize=5)
        
        heatmap
    end
    savedir(bn) = plotsdir(prototype_name, sim_name, mods_filename(; fixed_mods..., target_mods_nt...), bn)
    mkpath(savedir("") |> dirname)
    
    layout[:,1] = LText.(scene, mod_names_str[1:end-1], tellheight=false, rotation=pi/2)
    layout[end+1,2:end] = LText.(scene, mod_names_str[2:end], tellwidth=false)
    layout[0, :] = LText(scene, title, textsize = 30)
    cbar = layout[2:end-1, end+1] = LColorbar(scene, heatmaps[1], label = "Proportion")
    cbar.width = 25
    path = savedir("slices.png")
    Makie.save(path, scene)
    
    these_mods =  (save_idxs=nothing, other_opts=Dict(), fixed_mods..., target_mods_nt...)
    @show these_mods
    prototype = get_prototype(prototype_name)
    (mod_name, exec) = execute_single_modification(prototype, these_mods)
    #exec = execute(line_prototype(;mods..., other_opts=Dict()))
    #wp = ExecutionClassifications(exec.solution)
    #@show wp.has_propagation
    @warn "Not validating has_propagation"
    anim_path = savedir("sim_animation.mp4")
    animate_execution(anim_path, exec);
    
    scene, layout = Simulation73.exec_heatmap_slices(exec, 5, (1400,1000))
    layout[0,:] = LText(scene, "Simulation: $(target_mods_nt)\n Inhibition blocking threshold: $(these_mods[:blocking_θI])") 
    sim_heatmap_path = savedir("sim_heatmap.png")
    Makie.save( sim_heatmap_path, scene)
    
    scene, layout = Simulation73.exec_heatmap(exec)
    layout[0,:] = LText(scene, "Simulation: $(target_mods_nt)\n Inhibition blocking threshold: $(these_mods[:blocking_θI])") 
    sim_heatmap_path = savedir("sim_E_heatmap.png")
    Makie.save( sim_heatmap_path, scene)

end
