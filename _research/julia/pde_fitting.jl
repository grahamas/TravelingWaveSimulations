# ---
# jupyter:
#   jupytext:
#     formats: ipynb,julia//jl
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.0
#   kernelspec:
#     display_name: Julia 1.4.1
#     language: julia
#     name: julia-1.4
# ---

using FindPDE, TravelingWaveSimulations, Simulation73

example_name = "reduced_line_dos_effectively_sigmoid"
line_example = get_example(example_name)
exec = execute(line_example(; other_opts=Dict()));

diff(exec.solution.u)


