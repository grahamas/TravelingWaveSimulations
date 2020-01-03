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
using TravelingWaveSimulations

# %%
function search_for_traveling_behavior(example, mod_bounds)
    # Need to map from arry to kwargs
    # Should properly be done at the DEProblem interface...
    minfn(p) = 1 - tw_score(execute(example(; p...)))
    result = optimize(minfn, starting_mods, mod_bounds, SAMIN)
end
