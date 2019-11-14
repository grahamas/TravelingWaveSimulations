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
using TravelingWaveSimulations, NeuralModels, WilsonCowanModel, Simulation73, Plots

# %%
space = PeriodicLattice{Float64,2}(; n_points=(128,128), extent=(700.0, 700.0))
stimulus_p = SharpBumpStimulusParameter(; strength=1.2, width=28.1, time_windows=[(0.0,5.0)])
stimulus = stimulus_p(space)
connectivity_p = GaussianConnectivityParameter(; amplitude=24.0, spread=(25.0, 25.0))

heatmap(stimulus.bump_frame, xlab="x-space", ylab="y-space") |> display
savefig("stimulus_example.png")
heatmap(NeuralModels.kernel(connectivity_p, space), xlab="x-space", ylab="y-space") |> display
savefig("gaussian_connectivity_example.png")

# %%
