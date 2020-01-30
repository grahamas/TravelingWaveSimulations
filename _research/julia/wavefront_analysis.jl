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
using Simulation73, NeuralModels, TravelingWaveSimulations, Plots, Optim, LinearAlgebra, Distances, Statistics,
    IterTools, Combinatorics, DataFrames, GLM, JuliaDB

equals_str(key,val) = "$key=$val"
equals_strs(mods) = [equals_str(p...) for p in pairs(mods)]
mods_filename(x) = join(equals_strs(x), "_")

# %%
function velocity_results(results)
    (coef(results[:apex_loc])[1], deviance(results[:apex_loc]))
end

# %%
function find_first_satisfying_execution(mdb, example, dict_min=Dict(), dict_max=Dict())
    function filter_fn(row)
        above_mins = [row[key] >= val for (key, val) in pairs(dict_min)]
        below_maxes = [row[key] <= val for (key, val) in pairs(dict_max)]
        return all(above_mins) && all(below_maxes)
    end     
    for db in mdb
        satisfactory_rows = filter(filter_fn, db)
        if length(satisfactory_rows) > 0
            mods = JuliaDB.select(satisfactory_rows, Keys())
            @show "found $(mods[1])"
            sols = JuliaDB.select(satisfactory_rows, JuliaDB.Not(Keys()))
            return Execution(example(;mods[1]...), BareSolution(; pairs(sols[1])...))
        end
    end
    return nothing    
end

# %%
# Load most recent simulation data
data_root = joinpath(homedir(), "sim_data")
(example, mdb) = TravelingWaveSimulations.load_data(data_root, "sigmoid_normal_fft", 5);
example_name = TravelingWaveSimulations.get_example_name(mdb.fns[1])
sim_name = TravelingWaveSimulations.get_sim_name(mdb.fns[1])

# %%
# Analyse and extract twscore
# mdb_execs = MultiDBExecIter(example, dbs, ())
GC.gc()
mods = TravelingWaveSimulations.get_mods(mdb)
mod_names = keys(mods) |> collect
mod_values = values(mods) |> collect
A_is_traveling = Array{Bool}(undef, length.(values(mods))...)
A_velocity= Array{Union{Float64,Missing}}(undef, length.(values(mods))...)
A_velocity_errors= Array{Union{Float64,Missing}}(undef, length.(values(mods))...)
for db in mdb
    for (this_mod, exec) in DBExecIter(example, db, ())
        this_mod_key = keys(this_mod)
        this_mod_val = values(this_mod)
        A_idx = TravelingWaveSimulations.mod_idx(this_mod_key, this_mod_val, mod_names, mod_values)
        (results, score, tw_df) = tw_metrics(SolitaryWaveformMetrics, exec);
        if results === nothing || (:apex_loc ∉ keys(results))
            A_is_traveling[A_idx] = false
            A_velocity[A_idx] = missing
            A_velocity_errors[A_idx] = missing
        else
            A_is_traveling[A_idx] = true
            A_velocity[A_idx], sse = velocity_results(results)
            A_velocity_errors[A_idx] = sse / length(score) #make sum into mean
        end
    end
end
@show sum(ismissing.(A_velocity))
@show prod(size(A_velocity))

# %%
function mean_skip_missing(A; dims)
    missings = ismissing.(A)
    zeroed = copy(A)
    zeroed[missings] .= 0
    nonmissingsum = sum(zeroed; dims=dims)
    nonmissingmean = nonmissingsum ./ sum(.!missings; dims=dims)
    return nonmissingmean
end
all_dims = 1:length(mod_names)
for (x,y) in IterTools.subsets(all_dims, Val{2}())
    collapsed_dims = Tuple(setdiff(all_dims, (x,y)))
    is_traveling = dropdims(mean_skip_missing(A_is_traveling, dims=collapsed_dims), dims=collapsed_dims)
    velocities = dropdims(mean_skip_missing(A_velocity, dims=collapsed_dims), dims=collapsed_dims)
    velocity_errors = dropdims(mean_skip_missing(A_velocity_errors, dims=collapsed_dims), dims=collapsed_dims)
    prop_notmissing = dropdims(mean(.!ismissing.(A_velocity), dims=collapsed_dims), dims=collapsed_dims)
    plot(
        #heatmap(mod_values[x], mod_values[y], is_traveling, xlab=mod_names[x], ylab=mod_names[y], title="\"peakiness\" avgd across other spreads"),
        heatmap(mod_values[y], mod_values[x], velocities, xlab=mod_names[y], ylab=mod_names[x], title="velocity avgd"),
        heatmap(mod_values[y], mod_values[x], velocity_errors, xlab=mod_names[y], ylab=mod_names[x], title="error"),
        heatmap(mod_values[y], mod_values[x], prop_notmissing, xlab=mod_names[y], ylab=mod_names[x], title="prop not missing")
        ) |> display
    path = "wavefront_tmp/$(example_name)/$(sim_name)/$(mod_names[x])_$(mod_names[y])_centerfiterror.png"
    mkpath(dirname(path))
    png(path)
end

# %%
dict_max = Dict()
dict_min = Dict(:See => 40.0, :Sei => 100.0)
mod_names = keys(TravelingWaveSimulations.get_mods(mdb))
test_exec = find_first_satisfying_execution(mdb, example, dict_min, dict_max);

# %%
anim = custom_animate(test_exec)
mp4(anim, "wavefront_tmp/$(example_name)/$(sim_name)/anim_$(mods_filename(mods)).mp4")
