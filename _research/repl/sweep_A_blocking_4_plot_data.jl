
let classifications = classifications_A
    for name in keys(classifications)
        scene = sweep_2D_slice_heatmaps(classifications[name], title = "$(name) proportion")
        path = plotsdir("$(example_name)/$(sim_name)/$(name)_slices.png")
        mkpath(path |> dirname)
        Makie.save(path, scene)
    end
end

# +
# One plot
#using Makie, MakieLayout, DrWatson
#
#let mod_names = string.(mod_names)
#    plot_color = :magma
#    all_dims = 1:length(mod_names)
#    slices_2d = IterTools.subsets(all_dims, Val{2}())
#    plot_side_size = 350 * (length(all_dims) - 1)
#    plot_size = (plot_side_size, plot_side_size)
#    epileptic_scene, epileptic_layout = layoutscene(resolution=plot_size)
#    tw_scene, tw_layout = layoutscene(resolution=plot_size, title="Traveling Solitary Waves")
#
#    heatmap_pairs = map(slices_2d) do (x,y)
#        (x,y) = x < y ? (x,y) : (y,x)
#        collapsed_dims = Tuple(setdiff(all_dims, (x,y)))
#        mean_is_traveling_solitary = dropdims(mean_skip_missing(A_is_traveling_solitary, dims=collapsed_dims), dims=collapsed_dims) |> collect
#        mean_is_epileptic = dropdims(mean_skip_missing(A_is_epileptic, dims=collapsed_dims), dims=collapsed_dims) |> collect
#        my = mod_values[y] |> collect
#        mx = mod_values[x] |> collect
#        
#        @assert size(mean_is_epileptic) == length.((mod_values[x], mod_values[y]))
#        
#        epileptic_layout[x,y] = epileptic_ax = LAxis(epileptic_scene); 
#        tw_layout[x,y] = tw_ax = LAxis(tw_scene)
#        tightlimits!.([epileptic_ax, tw_ax])
#        (Makie.heatmap!(epileptic_ax, my, mx, mean_is_epileptic', colorrange=(0,1)),
#            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop epileptic"),
#        Makie.heatmap!(tw_ax, my, mx, mean_is_traveling_solitary', colorrange=(0,0.100001)))
#            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop traveling solitary")
#    end
#    ep_heatmaps, tw_heatmaps = zip(heatmap_pairs...)
#    epileptic_layout[:,1] = LText.(epileptic_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
#    epileptic_layout[end+1,2:end] = LText.(epileptic_scene, mod_names[2:end], tellwidth=false)
#    epileptic_layout[0, :] = LText(epileptic_scene, "Parameter sweep:\n Prop. simulations with traveling front", textsize = 30)
#    ep_cbar = epileptic_layout[2:end-1, end+1] = LColorbar(epileptic_scene, ep_heatmaps[1], label = "Proportion traveling fronts")
#    ep_cbar.width = 25
#    ep_path = plotsdir("$(example_name)/$(sim_name)/epileptic_slices.png")
#    mkpath(ep_path |> dirname)
#    Makie.save(ep_path, epileptic_scene)
#    
#    tw_layout[:,1] = LText.(tw_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
#    tw_layout[end+1,2:end] = LText.(tw_scene, mod_names[2:end], tellwidth=false)
#    tw_layout[0, :] = LText(tw_scene, "Parameter sweep:\n Prop. simulations with traveling wave", textsize = 30)
#    tw_cbar = tw_layout[2:end-1, end+1] = LColorbar(tw_scene, tw_heatmaps[1], label = "Proportion traveling waves")
#    tw_cbar.width = 25
#    tw_path = plotsdir("$(example_name)/$(sim_name)/tw_slices.png")
#    mkpath(tw_path |> dirname)
#    Makie.save(tw_path, tw_scene)
#end
