
function mean_skip_missing(A::AbstractArray; dims)
    missings = ismissing.(A)
    zeroed = copy(A)
    zeroed[missings] .= 0
    nonmissingsum = sum(zeroed; dims=dims)
    nonmissingmean = nonmissingsum ./ sum(.!missings; dims=dims)
    return nonmissingmean
end

function avg_across_dims(arr, dims)
    avgd = mean_skip_missing(arr, dims=dims)
    squeezed = dropdims(avgd, dims=dims)
    return collect(squeezed)
end

let classifications = classifications_A, mod_names = string.(mod_names)
    all_dims = 1:length(mod_names)
    for (x,y) in IterTools.subsets(all_dims, Val{2}())
        classification_names = Tuple(keys(classifications))
        collapsed_dims = Tuple(setdiff(all_dims, (x,y)))
        # FIXME: don't include in calculation if not sane
        mean_class = NamedTuple{classification_names}([avg_across_dims(a, collapsed_dims) for a in classifications])
        
        @assert size(mean_class[1]) == length.((mod_values[x], mod_values[y]))
        
        for name in classification_names
            scene, layout = layoutscene(resolution=(600, 600))
            #layout = GridLayout()
            layout[1,2] = ax = LAxis(scene); 
            heatmap = Makie.heatmap!(ax, mod_values[x], mod_values[y], mean_class[name])#, ylabel=mod_names[y], xlabel=mod_names[x], title=string(name))
            min_val, max_val = extrema(mean_class[name])
            layout[1,1] = LText(scene, mod_names[y], tellheight=false, rotation=pi/2)
            layout[2,2] = LText(scene, mod_names[x], tellwidth=false)
            if min_val != max_val
                layout[1,3] = cbar = LColorbar(scene, heatmap)
                cbar.width = 25
            end
            layout[0,:] = LText(scene, string(name))
        
            tightlimits!.([ax])#, tw_ax])
        
            path = plotsdir("$(example_name)/$(sim_name)/single_slices/$(mod_names[x])_$(mod_names[y])_$(name)_slice.png")
            mkpath(dirname(path))
            Makie.save(path, scene)
        end
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
