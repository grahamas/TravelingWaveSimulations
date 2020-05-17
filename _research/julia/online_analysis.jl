# -*- coding: utf-8 -*-
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

using AxisIndices

]update

# +
# %%
using Revise, DrWatson
 
using Simulation73, NeuralModels, TravelingWaveSimulations, LinearAlgebra, Distances, Statistics, AxisIndices,
    IterTools, Combinatorics, DataFrames, GLM, JuliaDB, DifferentialEquations
using Makie, MakieLayout
# using Plots, PlotThemes
# theme(:ggplot2)

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
# -

using Distributed
using Revise
using IterTools, Statistics, Plots, JuliaDB
using Simulation73, TravelingWaveSimulations
addprocs(5)
@everywhere using Revise
@everywhere using DrWatson
@everywhere using AxisIndices
@everywhere quickactivate(@__DIR__, "TravelingWaveSimulations")
@everywhere using TravelingWaveSimulations


function find_and_reconstruct_first_satisfying_execution(mdb, example, dict_bounds=Dict(); exec_opts...)
    @show pairs(dict_bounds)
    function filter_fn(row)
        within_bounds = [row[key] >= val[1] && row[key] <= val[2] for (key, val) in pairs(dict_bounds)]
        return all(within_bounds)
    end     
    for db in mdb
        satisfactory_rows = filter(filter_fn, db)
        if length(satisfactory_rows) > 0
            mods = JuliaDB.select(satisfactory_rows, Keys())
            @show "found $(mods[1])"
            sols = JuliaDB.select(satisfactory_rows, JuliaDB.Not(Keys()))
            @show sols[1]
            return (mods[1], execute(example(;mods[1]..., exec_opts...)))
        end
    end
    @warn "Could not find sim matching $dict_bounds"
    return nothing    
end

function solve_as_fronts(simulation::Simulation; solver_opts...)
    sv = SavedValues(Float64,Array{TravelingWaveSimulations.Wavefront{Float64,Float64,TravelingWaveSimulations.Value{Float64,Float64}},1})
    function save_func(u, t, integrator)
        # Need to get x from integrator (or simulation)
        sub_idx = integrator.opts.save_idxs
        sub_u = u[sub_idx]
        sub_x = [x[1] for x in simulation.space.arr[population(sub_idx,1)]]
        fronts = TravelingWaveSimulations.substantial_fronts(sub_u, sub_x)
        return fronts
    end                
    cb = SavingCallback(save_func, sv)
    sol = solve(simulation; callback=cb, solver_opts...)
    return (sv, sol)
end

example = get_example("line_dos_effectively_sigmoid")

# Load most recent simulation data
data_root = joinpath(homedir(), "data")#datadir())
(example, mdb) = TravelingWaveSimulations.load_data(data_root, 
    "reduced_line_dos_effectively_sigmoid", -6);
example_name = TravelingWaveSimulations.get_example_name(mdb.fns[1])
sim_name = TravelingWaveSimulations.get_sim_name(mdb.fns[1])

# +
#### Analyse and extract twscore
# mdb_execs = MultiDBExecIter(example, dbs, ())
GC.gc()
mods = TravelingWaveSimulations.get_mods(mdb)

all_mod_names = keys(mods) |> collect
all_mod_values = values(mods) |> collect
varied_mods = length.(all_mod_values) .> 1

mod_names = all_mod_names[varied_mods] |> sort
mod_values = [mods[name] for name in mod_names]
mod_dict = Dict(name => val for (name, val) in zip(mod_names, mod_values))

mod_array(T) = AxisIndicesArray(Array{Union{Bool,Missing}}(undef, length.(mod_values)...), Tuple(mod_values))
A_is_epileptic = mod_array(Bool)
A_is_traveling_solitary = mod_array(Bool)
A_is_decaying = mod_array(Bool)

@show A_is_epileptic.axes

# for db in mdb
#     for (this_mod, exec) in DBExecIter(example, db, ())
for (this_mod, this_result) in MultiDBRowIter(mdb)
#     exec = execute(example(; this_mod...))
#     wave_properties = TravelingWaveSimulations.get_wave_properties(exec)
    wave_properties = this_result[:wave_properties]
    A_idx = this_mod[mod_names]
    if wave_properties === missing
        A_is_epileptic[A_idx...] = missing
        A_is_traveling_solitary[A_idx...] = missing
        A_is_decaying[A_idx...] = missing
    else
        A_is_epileptic[A_idx...] = wave_properties.epileptic
        A_is_traveling_solitary[A_idx...] = wave_properties.traveling_solitary
        A_is_decaying[A_idx...] = wave_properties.decaying
    end
end

# -

mods

@show sum(ismissing.(A_is_epileptic))
@show prod(size(A_is_epileptic))
@show sum(skipmissing(A_is_epileptic)) 
@show sum(skipmissing(A_is_traveling_solitary)) 

# + jupyter={"outputs_hidden": true}
# Separate plots
plot_size = (600,300)
plot_color = :magma
all_dims = 1:length(mod_names)
scene, layout = GridLay
for (x,y) in IterTools.subsets(all_dims, Val{2}())
    collapsed_dims = Tuple(setdiff(all_dims, (x,y)))
    mean_is_traveling_solitary = dropdims(mean_skip_missing(A_is_traveling_solitary, dims=collapsed_dims), dims=collapsed_dims)
    mean_is_epileptic = dropdims(mean_skip_missing(A_is_epileptic, dims=collapsed_dims), dims=collapsed_dims)
    #prop_notmissing = dropdims(mean(.!ismissing.(A_is_epileptic), dims=collapsed_dims), dims=collapsed_dims)
    plot(
        #heatmap(mod_values[x], mod_values[y], mean_traveling, xlab=mod_names[x], ylab=mod_names[y], title="\"peakiness\" avgd across other spreads"),
#         heatmap(mod_values[y], mod_values[x], velocities, xlab=mod_names[y], ylab=mod_names[x], title="velocity avgd"),
#         heatmap(mod_values[y], mod_values[x], velocity_errors, xlab=mod_names[y], ylab=mod_names[x], title="error"),
        #heatmap(mod_values[y], mod_values[x], prop_notmissing, xlab=mod_names[y], ylab=mod_names[x], title="prop not missing"),
        heatmap(mod_values[y], mod_values[x], mean_is_epileptic, xlab=mod_names[y], ylab=mod_names[x], size=plot_size, color=plot_color, title="prop epileptic"),
        heatmap(mod_values[y], mod_values[x], mean_is_traveling_solitary, xlab=mod_names[y], ylab=mod_names[x], size=plot_size, color=plot_color, title="prop traveling solitary")
        #layout = (1,3)
        ) |> display
    path = "wavefront_tmp/$(example_name)/$(sim_name)/$(mod_names[x])_$(mod_names[y])_slice.png"
    mkpath(dirname(path))
    png(path)
end
# -

typeof(string.(mod_names))

# +
# One plot
using Makie, MakieLayout, DrWatson

let mod_names = string.(mod_names)
    plot_color = :magma
    all_dims = 1:length(mod_names)
    slices_2d = IterTools.subsets(all_dims, Val{2}())
    plot_side_size = 300 * (length(all_dims) - 1)
    plot_size = (plot_side_size, plot_side_size)
    epileptic_scene, epileptic_layout = layoutscene(resolution=plot_size)
    tw_scene, tw_layout = layoutscene(resolution=plot_size, title="Traveling Solitary Waves")

    heatmap_pairs = map(slices_2d) do (x,y)
        (x,y) = x < y ? (x,y) : (y,x)
        collapsed_dims = Tuple(setdiff(all_dims, (x,y)))
        mean_is_traveling_solitary = dropdims(mean_skip_missing(A_is_traveling_solitary, dims=collapsed_dims), dims=collapsed_dims) |> collect
        mean_is_epileptic = dropdims(mean_skip_missing(A_is_epileptic, dims=collapsed_dims), dims=collapsed_dims) |> collect
        my = mod_values[y] |> collect
        mx = mod_values[x] |> collect
        epileptic_layout[x,y] = epileptic_ax = LAxis(epileptic_scene); 
        tw_layout[x,y] = tw_ax = LAxis(tw_scene)
        (Makie.heatmap!(epileptic_ax, mx, my, mean_is_epileptic, colorrange=(01,1)),
            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop epileptic"),
        Makie.heatmap!(tw_ax, mx, my, mean_is_traveling_solitary))
            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop traveling solitary")
    end
    epileptic_layout[:,1] = LText.(epileptic_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
    epileptic_layout[end+1,2:end] = LText.(epileptic_scene, mod_names[2:end], tellwidth=false)
    epileptic_layout[0, :] = LText(epileptic_scene, "Traveling fronts", textsize = 30)
    ep_cbar = epileptic_layout[2:end-1, end+1] = LColorbar(scene, heatmap_pairs[1][1], label = "Activity Level")
    ep_cbar.width = 25
    ep_path = plotsdir("$(example_name)/$(sim_name)/epileptic_slices.png")
    mkpath(ep_path |> dirname)
    Makie.save(ep_path, epileptic_scene)
    
    tw_layout[:,1] = LText.(tw_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
    tw_layout[end+1,2:end] = LText.(tw_scene, mod_names[2:end], tellwidth=false)
    tw_layout[0, :] = LText(tw_scene, "Traveling waves", textsize = 30)
    tw_path = plotsdir("$(example_name)/$(sim_name)/tw_slices.png")
    mkpath(tw_path |> dirname)
    Makie.save(tw_path, tw_scene)
end
@show mod_names

# +
# %%
dict_bounds = Dict(:firing_θI => (8.0, 8.), :blocking_θI => (10.0, 11.0), 
                   :firing_θE => (-Inf, 7.9), :blocking_θE => (17.0, 22.0))
mod_names = keys(TravelingWaveSimulations.get_mods(mdb))
test_mods, test_exec = find_and_reconstruct_first_satisfying_execution(mdb, example, dict_bounds; other_opts=Dict());
# rerun = execute(test_exec.simulation)

# %%
wp = TravelingWaveSimulations.get_wave_properties(test_exec)
@show wp.epileptic
@show wp.traveling_solitary
anim = custom_animate(test_exec)
mp4(anim, "wavefront_tmp/$(example_name)/$(sim_name)/anim_$(mods_filename(mods)).mp4")
# -

example_name = "reduced_line_dos_effectively_sigmoid"
line_example = get_example(example_name)
mods = (Aee=200.0,Aei=200.0, Aie=200.0, blocking_θE=25.0,blocking_θI=10.0,firing_θE=6.0,firing_θI=7.0)
(mod_name, exec) = TravelingWaveSimulations.execute_single_modification(line_example, mods)
#exec = execute(line_example(;mods..., other_opts=Dict()))
wp = TravelingWaveSimulations.get_wave_properties(exec)
@show wp.epileptic
@show wp.traveling_solitary
anim = custom_animate(exec)
mp4(anim, "wavefront_tmp/$(example_name)/line_anim_$(mods_filename(mods)).mp4")

these_data = TravelingWaveSimulations.extract_data_namedtuple(exec)

pfronts = TravelingWaveSimulations.persistent_fronts(TravelingWaveSimulations.all_fronts(exec), exec.solution.t)
anim = custom_animate(exec, pfronts)
mp4(anim, "wavefront_tmp/$(example_name)/fronts_line_anim_$(mods_filename(mods)).mp4")

TravelingWaveSimulations.get_velocities(pfronts[5])

TravelingWaveSimulations.based_on_example_NO_PARALLEL(; example_name="reduced_line_dos_effectively_sigmoid", 
    modifications=["blocking_θI=6.0:30.0", "blocking_θE=6.0:30.0"], 
    data_root=projectdir("tmp"))

based_on_example(; example_name="reduced_line_dos_effectively_sigmoid", modifications=["stim_strength=0.4:0.1:2.0", "stim_width=20.0:5.0:50.0", "Aee=22.0", "Aei=20.2"], max_batch_size=5000, data_root=projectdir("tmp"), max_sims_in_mem=30001, progress=true)




using FindPDE


