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
    IterTools, Combinatorics, DataFrames, GLM, JuliaDB, DifferentialEquations

equals_str(key,val) = "$key=$val"
equals_strs(mods) = [equals_str(p...) for p in pairs(mods)]
mods_filename(x) = join(equals_strs(x), "_")

# %%
function velocity_results(results)
    (coef(results[:apex_loc])[1], deviance(results[:apex_loc]))
end

function mean_skip_missing(A; dims)
    missings = ismissing.(A)
    zeroed = copy(A)
    zeroed[missings] .= 0
    nonmissingsum = sum(zeroed; dims=dims)
    nonmissingmean = nonmissingsum ./ sum(.!missings; dims=dims)
    return nonmissingmean
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
            return (mods[1], Execution(example(;mods[1]...), BareSolution(; pairs(sols[1])...)))
        end
    end
    return nothing    
end

# %%
function is_traveling(l_frame_fronts::Array{<:Array,1}, t, min_motion)
    rightmost_locations = [frame_fronts[end].slope.loc for frame_fronts in l_frame_fronts]
    rightward_motion = diff(rightmost_locations) ./ diff(t)
    return all(rightward_motion[(length(rightward_motion) ÷ 2):end-5] .> min_motion)
end

function persistent_activation(l_frames, t, min_activation)
    # Test if the left-most value is ever high and non-decreasing
    leftmost_value = [frame[1] for frame in l_frames[(length(t) ÷ 2):end]]
    d_leftmost_value = diff(leftmost_value) ./ diff(t[(length(t) ÷ 2):end])
    # test leftmost is not ONLY decreasing
    low_enough = leftmost_value[2:end] .< min_activation
    decreasing = d_leftmost_value .< 0.0
    return !all(low_enough .| decreasing)
end



function is_epileptic_radial_slice(exec, min_motion=1e-2, min_activation=5e-2)
    # We'll call it epileptic if:
    #  1. There is a traveling front
    #  2. There is persistently elevated activity following the front.
    # How to deal with oscillatory activity? 
    # Want persistent oscillation to be considered epileptic, even if sometimes zero
    # In practice we'll ignore this for now --- FIXME
    l_frames = exec.solution.u
    t = timepoints(exec)
    x = [x[1] for x in space(exec).arr]
    l_frame_fronts = TravelingWaveSimulations.substantial_fronts.(l_frames, Ref(x))
    persistent_activation(l_frames, t, min_activation) && (mean(population(l_frames[end],1)) > 0.2 || is_traveling(l_frame_fronts, t, min_motion)) 
end

function is_traveling_solitary_radial_slice(exec, min_motion=1e-2, max_activation=5e-2)
    # We'll call it traveling solitary if:
    #  1. There is a traveling front
    #  2. No elevated activity trails the front
    l_frames = exec.solution.u
    t = timepoints(exec)
    x = [x[1] for x in space(exec).arr]
    l_frame_fronts = TravelingWaveSimulations.substantial_fronts.(l_frames, Ref(x))
    !persistent_activation(l_frames, t, max_activation) && is_traveling(l_frame_fronts, t, min_motion)
end

# %%
# Load most recent simulation data
data_root = joinpath(homedir(), "sim_data")
(example, mdb) = TravelingWaveSimulations.load_data(data_root, "line_dos_effectively_sigmoid");
example_name = TravelingWaveSimulations.get_example_name(mdb.fns[1])
sim_name = TravelingWaveSimulations.get_sim_name(mdb.fns[1])

# %%
# Analyse and extract twscore
# mdb_execs = MultiDBExecIter(example, dbs, ())
GC.gc()
mods = TravelingWaveSimulations.get_mods(mdb)

all_mod_names = keys(mods) |> collect
all_mod_values = values(mods) |> collect
varied_mods = length.(all_mod_values) .> 1

mod_names = all_mod_names[varied_mods]
mod_values = all_mod_values[varied_mods]

# A_velocity= Array{Union{Float64,Missing}}(undef, length.(values(mods))...)
# A_velocity_errors= Array{Union{Float64,Missing}}(undef, length.(values(mods))...)
# A_mean_N_fronts = Array{Union{Float64,Missing}}(undef, length.(values(mods))...)
A_is_epileptic = Array{Union{Bool,Missing}}(undef, length.(mod_values)...)
A_is_traveling_solitary = Array{Union{Bool,Missing}}(undef, length.(mod_values)...)

# for db in mdb
#     for (this_mod, exec) in DBExecIter(example, db, ())
for (this_mod, exec) in MultiDBExecIter(example, mdb, ())
    names_and_vals = map(mod_names) do name
        val = this_mod[name]
        return (name, val)
    end
    this_mod_key, this_mod_val = zip(names_and_vals...)
    A_idx = TravelingWaveSimulations.mod_idx(this_mod_key, this_mod_val, mod_names, mod_values)
    if exec isa Execution{Missing}
#         A_velocity[A_idx] = missing
#         A_velocity_errors[A_idx] = missing
#         A_mean_N_fronts[A_idx] = missing
        A_is_epileptic[A_idx] = missing
        A_is_traveling_solitary[A_idx] = missing
    else
        A_is_epileptic[A_idx] = is_epileptic_radial_slice(exec)
        A_is_traveling_solitary[A_idx] = is_traveling_solitary_radial_slice(exec)
    end
end


# %%
@show sum(ismissing.(A_is_epileptic))
@show prod(size(A_is_epileptic))
@show sum(A_is_epileptic) / prod(size(A_is_epileptic))
@show sum(A_is_traveling_solitary) / prod(size(A_is_traveling_solitary))

# %%

all_dims = 1:length(mod_names)
for (x,y) in IterTools.subsets(all_dims, Val{2}())
    collapsed_dims = Tuple(setdiff(all_dims, (x,y)))
    mean_is_traveling_solitary = dropdims(mean_skip_missing(A_is_traveling_solitary, dims=collapsed_dims), dims=collapsed_dims)
    mean_is_epileptic = dropdims(mean_skip_missing(A_is_epileptic, dims=collapsed_dims), dims=collapsed_dims)
    prop_notmissing = dropdims(mean(.!ismissing.(A_is_epileptic), dims=collapsed_dims), dims=collapsed_dims)
    plot(
        #heatmap(mod_values[x], mod_values[y], mean_traveling, xlab=mod_names[x], ylab=mod_names[y], title="\"peakiness\" avgd across other spreads"),
#         heatmap(mod_values[y], mod_values[x], velocities, xlab=mod_names[y], ylab=mod_names[x], title="velocity avgd"),
#         heatmap(mod_values[y], mod_values[x], velocity_errors, xlab=mod_names[y], ylab=mod_names[x], title="error"),
        heatmap(mod_values[y], mod_values[x], prop_notmissing, xlab=mod_names[y], ylab=mod_names[x], title="prop not missing"),
        heatmap(mod_values[y], mod_values[x], mean_is_epileptic, xlab=mod_names[y], ylab=mod_names[x], title="prop epileptic"),
        heatmap(mod_values[y], mod_values[x], mean_is_traveling_solitary, xlab=mod_names[y], ylab=mod_names[x], title="prop traveling solitary")
        #layout = (1,3)
        ) |> display
    path = "wavefront_tmp/$(example_name)/$(sim_name)/$(mod_names[x])_$(mod_names[y])_centerfiterror.png"
    mkpath(dirname(path))
    png(path)
end

# %%
dict_max = Dict()#:blocking_θI => 14.0, :blocking_θE => 14.0)
dict_min = Dict(:See => 34.)#, :blocking_θI => 10.0)
mod_names = keys(TravelingWaveSimulations.get_mods(mdb))
test_mods, test_exec = find_first_satisfying_execution(mdb, example, dict_min, dict_max);

# %%
@show is_epileptic_radial_slice(test_exec)
@show is_traveling_solitary_radial_slice(test_exec)
anim = custom_animate(test_exec)
mp4(anim, "wavefront_tmp/$(example_name)/$(sim_name)/anim_$(mods_filename(mods)).mp4")

# %%
test_exec.solution.u[end][end,1]

# %%
@show is_traveling_solitary_radial_slice(test_exec)
rerun_anim = custom_animate(test_exec)
mp4(rerun_anim, "wavefront_tmp/$(example_name)/$(sim_name)/rerun_anim_$(mods_filename(mods)).mp4")

# %%
function add_noise(sim, noise_coef, stochastic_algorithm=SRIW1())
    model = sim.model
    noisy_model = NoisyInputModel(model, WeinerProcess(0.0, 0.5))
    Simulation(noisy_model; space=sim.space, tspan=sim.tspan, initial_value=sim.initial_value,
        dt=sim.dt, algorithm=stochastic_algorithm, sim.solver_options...)
end

noisy_example = add_noise(example(), 1.0)
noisy_exec = execute(noisy_example);

# %%
@show is_traveling_solitary_radial_slice(noisy_exec)
noisy_anim = custom_animate(noisy_exec)
mp4(noisy_anim, "wavefront_tmp/$(example_name)/$(sim_name)/noisy_anim_$(mods_filename(mods)).mp4")

# %% jupyter={"source_hidden": true}
dict_max = Dict(:blocking_aI => 1.0, :blocking_aE => 1.0)
dict_min = Dict(:blocking_aI => 1.0, :blocking_aE => 1.0)#, :blocking_θE => 10.0)#, :blocking_θI => 10.0)
test_mods, test_exec = find_first_satisfying_execution(mdb, example, dict_min, dict_max);
