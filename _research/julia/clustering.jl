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
using JuliaDB, OnlineStats, Plots, TSne, TravelingWaveSimulations, Random


# %%
plots_dir(x) = joinpath(@__DIR__, "..", "plots", x)

amplitude_dir_path = plots_dir("Aee=19.0:29.0;Aei=13.0:23.0;Aie=22.0:32.0;Aii=1.0:1.0:6.0_2019-10-17T22:33:44.392_v1.0-21-g16bf0a7_dirty")
amplitude_file_path = joinpath(amplitude_dir_path, "1.jdb")
spread_file_path = "/home/grahams/git/TravelingWaveSimulations/plots/See=7.0:7.0:77.0;Sei=7.0:7.0:77.0;Sie=7.0:7.0:77.0;Sii=7.0:14.0:77.0;n=128;x=700.0_2019-10-25T09:25:26.251_v1.0-22-g56dee0a_dirty/1.jdb"
broad_amplitude_file_path = "/home/grahams/git/TravelingWaveSimulations/plots/Aee=14.0:2.0:34.0;Aei=8.0:2.0:28.0;Aie=17.0:2.0:37.0;Aii=1.0:2.0:12.0_2019-10-22T09:16:30.843_v1.0-22-g56dee0a_dirty/1.jdb"
amplitude_sweep = load(amplitude_file_path)
spread_sweep = load(spread_file_path)
@show spread_sweep
broad_amplitude_sweep = load(broad_amplitude_file_path)

# %%
n_clusters = 4
amplitude_series = reduce(FTSeries(KMeans(2015,n_clusters), transform=(x) -> cat([timepoint[:,1] for timepoint in x]..., dims=1)), amplitude_sweep, select=5);
centers = [reshape(amplitude_series.stats[1].value[i].value, (65, 31)) for i in 1:n_clusters];

plot([heatmap(centers[i], title="($i)") for i in 1:n_clusters]..., layout=n_clusters, xlab="time", ylab="space (radial slice)") |> display


# %%
savefig("amplitude_sweep_4_cluster_centers.png")

# %%
n_clusters = 4
spread_series = reduce(FTSeries(KMeans(2015,n_clusters), transform=(x) -> cat([timepoint[:,1] for timepoint in x]..., dims=1)), spread_sweep, select=5);
spread_sweep_centers = [reshape(spread_series.stats[1].value[i].value, (65, 31)) for i in 1:n_clusters];
plot([heatmap(spread_sweep_centers[i], title="($i)") for i in 1:n_clusters]..., layout=n_clusters, xlab="time", ylab="space (radial slice)") |> display


# %%
savefig("spread_sweep_4_cluster_centers_WIDER_other.png")

# %%
using LinearAlgebra
function which_cluster(o::KMeans, x)
    fill!(o.buffer, 0.0)
    for k in eachindex(o.buffer)
        cluster = o.value[k]
        o.buffer[k] = norm(x .- cluster.value[k])
    end
    k_star = argmin(o.buffer)
    return k_star
end

using Lazy
db_to_E_mat(db) = cat(db_to_E_vecs(db)..., dims=2) # COLUMNS ARE DATA POINTS
db_to_E_vecs(db) = @as x db begin
    select(x, :u)
    @view x[randperm(length(x))]
    [cat([timepoint[:,1] for timepoint in timeseries]..., dims=1) for timeseries in x]
end

using LinearAlgebra
function tsne_db(db, save_name::Union{Nothing,String}=nothing; kwargs...)
    mat = db_to_E_mat(db) |> transpose
    out_tsne = tsne(mat, 2)
    scatter(out_tsne[:,1], out_tsne[:,2]; kwargs...) |> display    
    save_name != nothing && savefig(save_name)
    return out_tsne
end
function tsne_db_nt(db, save_name::Union{Nothing,String}=nothing; kwargs...)
    mat = db_to_E_mat(db)
    out_tsne = tsne(mat, 2)
    scatter(out_tsne[:,1], out_tsne[:,2]; kwargs...) |> display    
    save_name != nothing && savefig(save_name)
    return out_tsne
end
# function tsne_db(db, n_clusters::Int)
#     mat = db_to_E_mat(db)
#     out_tsne = tsne(mat, 2)
#     kmeans_clusters = kmeans(out_tsne, n_clusters)
#     scatter(out_tsne[:,1], out_tsne[:,2], color=kmeans_clusters.assignments)
#     save_name != nothing && savefig(save_name)
# end
# function tsne_db(db, clusters::Vector{Int}, save_name::Union{Nothing,String}=nothing)
#     mat = db_to_E_mat(db)
#     out_tsne = tsne(mat, 2)
#     kmeans_clusters = kmeans(out_tsne, n_clusters)
#     scatter(out_tsne[:,1], out_tsne[:,2], color=clusters)
#     save_name != nothing && savefig(save_name)
# end
# function tsne_db(db, clusters, save_name::Union{Nothing,String}=nothing; kwargs...)
#     vecs = db_to_E_vecs(db)
#     mat = cat(vecs..., dims=2)
#     out_tsne = tsne(mat, 2)
#     cluster_ids = which_cluster.(Ref(clusters), vecs)
#     scatter(out_tsne[:,1], out_tsne[:,2]; color=cluster_ids, kwargs...) |> display
#     save_name != nothing && savefig(save_name)
#     return (out_tsne, cluster_ids)
# end

# %%
spread_clusters = spread_series.stats[1];
amp_clusters = amplitude_series.stats[1];

# %%
amp_tsne = tsne_db(amplitude_sweep[1:7:end], "test_tsne.png"; title="Amplitude Sweep tSNE (19-29,13-23;-(22-32),-(1-6))")

# %%
spread_sweep

# %%
spr_tsne = tsne_db(spread_sweep[1:7:end], "spread_sweep_7_tsne_WIDER.png";  title="Spread Sweep tSNE (7-35um)")

# %% {"collapsed": true, "jupyter": {"outputs_hidden": true}}
broad_amp_tsne = tsne_db(broad_amplitude_sweep[1:7:end], "broad_amplitude_sweep_7_tsne.png"; title="Broader Amplitude Sweep tSNE (14:2:34,8:2:28;-(17:2:37),-(1:2:12))")

# %%
sum(ids .== 1)
