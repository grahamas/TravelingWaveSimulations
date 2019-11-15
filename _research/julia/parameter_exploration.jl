# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all
#     notebook_metadata_filter: all,-language_info
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Julia 1.3.0-rc4
#     language: julia
#     name: julia-1.3
# ---

# %%
using TravelingWaveSimulations, Simulation73, Plots, TSne, JuliaDB, Clustering

example=TravelingWaveSimulations.examples_dict["sigmoid_normal_fft"]


# %%
mdl = example(n=512, x=700.0, amplitude=([25.0 -25.2; 35.0 -4.0]), stop_time=30.0); execution=execute(mdl); ucat=cat([u[:,1] for u in execution.solution.u]..., dims=2); 
heatmap(ucat, xlab="time", ylab="space (1D slice from center at 0deg)") |> display
savefig("oscillating_512_slice.png")

# %% {"collapsed": true, "jupyter": {"outputs_hidden": true}}
mp4(custom_animate(execution), "oscillating_512_slice.mp4", fps=10)

# %%
mdl = example(n=512, x=700.0, amplitude=([25.0 -25.2; 35.0 -4.0]), stop_time=30.0, save_idxs=nothing); execution=execute(mdl)
mp4(custom_animate(execution), "oscillating_512_full.mp4", fps=10)

# %%
mdl = example(n=512, x=700.0, amplitude=([29.0 -13.0; 32.0 -1.0]), stop_time=30.0); execution=execute(mdl); ucat=cat([u[:,1] for u in execution.solution.u]..., dims=2); 
heatmap(ucat, xlab="time", ylab="space (1D slice from center at 0deg)") |> display
savefig("wavefront_512_slice.png")
mp4(custom_animate(execution), "wavefront_512_slice.mp4", fps=10)

# %%
mdl = example(n=512, x=700.0, amplitude=([29.0 -13.0; 32.0 -1.0]), stop_time=30.0, save_idxs=nothing); execution=execute(mdl)
mp4(custom_animate(execution), "wavefront_512_full.mp4", fps=10)

# %%
full_mdl_odd = example(n=512, x=700.0, amplitude=([25.0 -25.2; 35.0 -4.0]), stop_time=30.0, save_idxs=nothing); full_execution_odd=execute(full_mdl_odd)
mp4(custom_animate(full_execution_odd), "oscillating_wave_longer.mp4")

# %% {"collapsed": true, "jupyter": {"outputs_hidden": true}}
full_mdl_bad = example(n=128, x=700.0, amplitude=([29.0 -13.0; 32.0 -1.0]), stop_time=30.0, dt=1.0, save_idxs=nothing); full_execution_bad=execute(full_mdl_bad);  
mp4(custom_animate(full_execution_bad), "traveling_wavefront.mp4")

# %% {"collapsed": true, "jupyter": {"outputs_hidden": true}}
full_mdl_bad = example(n=127, x=700.0, amplitude=([29.0 -13.0; 32.0 -1.0]), stop_time=30.0, dt=1.0, save_idxs=nothing); full_execution_bad=execute(full_mdl_bad);  
mp4(custom_animate(full_execution_bad), "traveling_wavefront_127.mp4")

# %% {"collapsed": true, "jupyter": {"outputs_hidden": true}}
full_mdl_bad = example(n=512, x=700.0, amplitude=([29.0 -13.0; 32.0 -1.0]), stop_time=30.0, dt=1.0); full_execution_bad=execute(full_mdl_bad);  
mp4(custom_animate(full_execution_bad), "traveling_wavefront_slice_512.mp4")

# %%
mx = maximum(execution.solution)
for frame in execution.solution.u
    IJulia.clear_output(true)
    plot(frame, ylim=[0,mx]) |> display
end

# %%
mx = maximum(execution.solution)
for frame in execution.solution.u
    IJulia.clear_output(true)
    plot(frame, ylim=[0,mx]) |> display
end

# %%
using JuliaDB
# amp_path = "/home/grahams/git/TravelingWaveSimulations/plots/Aee=19.0:29.0;Aei=13.0:23.0;Aie=22.0:32.0;Aii=1.0:1.0:6.0_2019-10-17T22:33:44.392_v1.0-21-g16bf0a7_dirty/1.jdb"
# amp_data = load(amp_path)
# spr_path = "/home/grahams/git/TravelingWaveSimulations/plots/See=22.0:28.0;Sei=22.0:28.0;Sie=22.0:28.0;Sii=22.0:28.0_2019-10-17T23:16:46.175_v1.0-21-g16bf0a7_dirty/1.jdb"
# spr_data = load(spr_path)
broad_amp_path = "/home/grahams/git/TravelingWaveSimulations/plots/Aee=14.0:2.0:34.0;Aei=8.0:2.0:28.0;Aie=17.0:2.0:37.0;Aii=1.0:2.0:12.0_2019-10-22T09:16:30.843_v1.0-22-g56dee0a_dirty/1.jdb"
broad_amp_data = load(broad_amp_path)

# %%
using Lazy
db_to_E_mat(db) = @as x db begin
    select(x, :u)
    [cat([timepoint[:,1] for timepoint in timeseries]..., dims=1) for timeseries in x]
    cat(x..., dims=2)
    transpose(x)
end
function tsne_db(db)
    mat = db_to_E_mat(db)
    out_tsne = tsne(mat, 2)
    scatter(out_tsne[:,1], out_tsne[:,2]) |> display
    return out_tsne
end
function tsne_db(db, n_clusters)
    mat = db_to_E_mat(db)
    out_tsne = tsne(mat, 2)
    kmeans_clusters = kmeans(out_tsne, n_clusters)
    scatter(out_tsne[:,1], out_tsne[:,2], color=kmeans_clusters.assignments)
end
    

# %%
tsne_db(amp_data[1:5:end])

# %% {"scrolled": true}
tsne_db(spr_data[1:1800])


# %%
broad_amp_tsne = tsne_db(broad_amp_data[1:5:end])
