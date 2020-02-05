# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.3.1
#     language: julia
#     name: julia-1.3
# ---

# %%
using Revise
using Simulation73, NeuralModels, TravelingWaveSimulations, Plots, 
    LinearAlgebra, Distances, Statistics,
    IterTools, DataFrames, GLM, JuliaDB,
    DifferentialEquations

# %%
mtxconv_example = get_example("replicate_neuman")
exec_mtxconv = execute(mtxconv_example(; n = 1000, stop_time = 100.0, callback=nothing));

# %%
anim = custom_animate(exec_mtxconv)
mp4(anim, "tmp/replicate_neuman_anim.mp4")

# %%
fft_example = get_example("replicate_neuman_fft")
exec_fft = execute(fft_example(; n = 100, stop_time = 100.0, dt=0.1, nonlinearityI = SigmoidNonlinearity(; a=1.2, \theta)));
anim = custom_animate(exec_fft)
mp4(anim, "tmp/replicate_neuman_fft_anim.mp4")

# %%
A = map(1.0:0.1:30.0) do stim_strength
    exec = execute(fft_example(;
            α = (1.1, 1.0),
            β = (1.0, 1.0),
            τ = (10.0, 18.0),
            Aee = 16.0, Aei = 18.2,
            Aie = 27.0, Aii = 4.0,
            See = 25.0, Sei = 27.0,
            Sie = 27.0, Sii = 25.0,
            algorithm = Euler(),
            time_windows = [[(0.0, 55.0)], [(0.0, 55.0)]],
            stim_strength = stim_strength,
            stim_width = 2.81,
            common_baseline=0.1,
            stop_time=28.1,
            callback=nothing,
            x = 100, n = 101
            ))
    return maximum(exec.solution.u[t][:,1][end] for t in 1:length(exec.solution.u))
end |> maximum
