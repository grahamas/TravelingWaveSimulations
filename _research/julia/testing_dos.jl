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
dos_example = get_example("dos_normal_fft")
dos_exec = execute(dos_example(;  n=256, x=700.0, See=25.0, Sii=25.0, Sie=27.0, Sei=27.0,
                                Aee=250.0, Aei=75.0, Aie=50.0, Aii=10.0, strengthE=10.0, widthE=50.0,
                                algorithm=Tsit5()));

# %%
anim = custom_animate(dos_exec)
mp4(anim, "dos_tmp/dos_normal_fft/no_mods/anim.mp4")

# %%
data_root = joinpath(homedir(), "sim_data")
(example, mdb) = TravelingWaveSimulations.load_data(data_root, "dos_normal_fft");
example_name = TravelingWaveSimulations.get_example_name(mdb.fns[1])
sim_name = TravelingWaveSimulations.get_sim_name(mdb.fns[1])

# %%
# Analyse and extract twscore
# mdb_execs = MultiDBExecIter(example, dbs, ())
GC.gc()
mods = TravelingWaveSimulations.get_mods(mdb)
mod_names = keys(mods) |> collect
mod_values = values(mods) |> collect
A_tws = Array{Float64}(undef, length.(values(mods))...)
A_velocity= Array{Union{Float64,Missing}}(undef, length.(values(mods))...)
A_velocity_errors= Array{Union{Float64,Missing}}(undef, length.(values(mods))...)
for db in mdb
    for (this_mod, exec) in DBExecIter(example, db, ())
        this_mod_key = keys(this_mod)
        this_mod_val = values(this_mod)
        A_idx = TravelingWaveSimulations.mod_idx(this_mod_key, this_mod_val, mod_names, mod_values)
        (results, score, tw_df) = tw_metrics(SolitaryWaveformMetrics, exec);
        if tws === nothing
            A_is_traveling[A_idx] = false
            A_velocity[A_idx] = missing
            A_velocity_errors[A_idx] = missing
        else
            A_is_traveling[A_idx] = true
            A_velocity[A_idx], A_velocity_errors[A_idx] = velocity_results(results)
        end
    end
end
@show sum(ismissing.(A_velocity))
@show prod(size(A_velocity))

# %%
exec.solution.u |> size

# %%
