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
#     display_name: Julia 1.3.0
#     language: julia
#     name: julia-1.3
# ---

# %%
using TravelingWaveSimulations, Simulation73, Plots

# %%
example=TravelingWaveSimulations.examples_dict["sigmoid_normal_fft"]

# %%
execution = execute(example(n=512, x=700.0, amplitude=([25.0 -25.2; 35.0 -4.0]), stop_time=30.0));

# %%
anim = TravelingWaveSimulations.custom_animate(execution)
mp4(anim, "tmp/_tmp.mp4")
