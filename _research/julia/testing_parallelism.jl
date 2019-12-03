# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all
#     notebook_metadata_filter: all,-language_info
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Julia 1.3.0
#     language: julia
#     name: julia-1.3
# ---

# %%
using TravelingWaveSimulations, Simulation73, Plots, JuliaDB

# %%
example=TravelingWaveSimulations.examples_dict["sigmoid_normal_fft"]

# %%
execution = execute(example(n=512, x=700.0, amplitude=([25.0 -25.2; 35.0 -4.0]), stop_time=30.0));

# %%
anim = TravelingWaveSimulations.custom_animate(execution)
mp4(anim, "tmp/_tmp.mp4")

# %%
failures = based_on_example(; example_name="sigmoid_normal_fft", data_root="data", modifications=["See=1.0:10.0"], batch=5)

# %%
using Distributed;  @show nworkers(); addprocs(1); @show nworkers()
@everywhere using Pkg
@everywhere Pkg.activate("..")
@everywhere using TravelingWaveSimulations
@show nworkers()


# %%
# small parallel example
failures = based_on_example(; example_name="sigmoid_normal_fft", data_root="data", modifications=["See=1.0:10.0"], batch=5)

# %%
function most_recent_subdir(datadir)
    subdirs = joinpath.(Ref(datadir), readdir(datadir))
    sort(subdirs, by=mtime, rev=true)[1]
end

# %%
recent_dir = most_recent_subdir("data/sigmoid_normal_fft")
d1 = load(joinpath(recent_dir, "1.jdb"))
d2 = load(joinpath(recent_dir, "2.jdb"))

# %%
recent_dir = most_recent_subdir("$(ENV["HOME"])/sims")
recent_db = load(joinpath(recent_dir, "1.jdb"))

# %%
recent_db[1].u
