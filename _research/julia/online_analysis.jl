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

]update

# +
# %%
using Revise, DrWatson
 
using Simulation73, NeuralModels, TravelingWaveSimulations, LinearAlgebra, Distances, Statistics, AxisIndices,
    IterTools, Combinatorics, DataFrames, GLM, JuliaDB, DifferentialEquations
using Makie, MakieLayout, AbstractPlotting
AbstractPlotting.inline!(true)
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

# +
#### FROM plots.ipynb
### DO NOT MODIFY HERE ###

function exec_heatmap_slices(exec::AbstractExecution, n_slices=5, resolution=(1600,1200))
    scene, layout = layoutscene(resolution=resolution)
    
    # adding timepoint slices
    soln = exec.solution
    t = soln.t
    xs = coordinate_axes(Simulation73.reduced_space(exec))[1] |> collect
    pop_names = exec.simulation.model.pop_names
    pop_idxs = 1:length(pop_names)
    
    n_x, n_p, n_t = size(soln)
    step = (length(soln.t)) ÷ n_slices
    t_idxs = 2:step:length(soln.t)
    
    hm_axes = [LAxis(scene, title = "$pop_name activity", aspect=1.0) for pop_name in pop_names]
    heatmaps = map(1:length(pop_names)) do idx_pop
        ax = hm_axes[idx_pop]
        pop_activity = cat(population.(soln.u, idx_pop)..., dims=2)
        heatmap!(ax, t, xs, pop_activity')
    end
    tightlimits!.(hm_axes)
    linkaxes!(hm_axes...)
    hideydecorations!.(hm_axes[2:end])
    
    layout[2,pop_idxs] = map(pop_idxs) do pop_idx
        slices_layout = GridLayout(rowsizes=[Auto()], alignmode=Outside(10), tellheight=true)
        slice_axes = slices_layout[:h] = map(t_idxs) do t_idx
            ax = LAxis(scene, aspect=1.0, tellheight=true)
            lines!(ax, xs, soln[:,pop_idx,t_idx])
            tightlimits!(ax)
            hideydecorations!(ax)
            hidexdecorations!(ax)
            ax
        end
        linkaxes!(slice_axes...)
        slices_layout[2,1:length(t_idxs)] = [LText(scene, "t=$(round(time, digits=1))", textsize=16, tellwidth=false) for time in t[t_idxs]]
        trim!(slices_layout)
        slices_layout        
    end
    layout[1,pop_idxs] = hm_axes
    cbar = layout[1, end+1] = LColorbar(scene, heatmaps[1], label = "Activity Level")
    cbar.width = 25

    ylabel = layout[1,0] = LText(scene, "space (μm)", rotation=pi/2, tellheight=false)
    xlabel = layout[end+1,2:3] = [ LText(scene, "time (ms)", tellwidth=false), LText(scene, "time (ms)", tellwidth=false)]
    return (scene, layout)
end

# + jupyter={"source_hidden": true}
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


# + jupyter={"source_hidden": true}
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

# + jupyter={"source_hidden": true}
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
# -

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
fixed_mods = length.(all_mod_values) .== 1
fixed_mods_dict = Dict(name => mods[name] for name in all_mod_names[fixed_mods])

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

# +
# Separate plots
all_dims = 1:length(mod_names)
let mod_names = string.(mod_names)
    for (x,y) in IterTools.subsets(all_dims, Val{2}())
        @show (x,y)
        scene, layout = layoutscene(resolution=(600, 600))
        collapsed_dims = Tuple(setdiff(all_dims, (x,y)))
        mean_is_traveling_solitary = dropdims(mean_skip_missing(A_is_traveling_solitary, dims=collapsed_dims), dims=collapsed_dims) |> collect
        mean_is_epileptic = dropdims(mean_skip_missing(A_is_epileptic, dims=collapsed_dims), dims=collapsed_dims) |> collect
        
        @assert size(mean_is_epileptic) == length.((mod_values[x], mod_values[y]))
        
        ep_layout = GridLayout()
        ep_layout[1,2] = epileptic_ax = LAxis(scene); 
        ep_heatmap = Makie.heatmap!(epileptic_ax, mod_values[x], mod_values[y], mean_is_epileptic')#, ylabel=mod_names[y], xlabel=mod_names[x])#, title="prop epileptic"),
        ep_layout[1,3] = ep_cbar = LColorbar(scene, ep_heatmap)
        ep_layout[1,1] = LText(scene, mod_names[y], tellheight=false, rotation=pi/2)
        ep_layout[2,2] = LText(scene, mod_names[x], tellwidth=false)
        ep_cbar.width = 25
        #ep_layout[0,:] = LText(scene, "fronts")
        
#         tw_layout = GridLayout()
#         tw_layout[1,1] = LText(scene, mod_names[y], tellheight=false, rotation=pi/2)
#         tw_layout[1,2] = tw_ax = LAxis(scene); 
#         tw_layout[2,2] = LText(scene, mod_names[x], tellwidth=false)
#         tw_heatmap = Makie.heatmap!(tw_ax, mod_values[x], mod_values[y], mean_is_traveling_solitary')#, ylabel=mod_names[y], xlabel=mod_names[x])#, title="prop epileptic"),
#         tw_layout[1,3] = tw_cbar = LColorbar(scene, tw_heatmap)
#         tw_cbar.width = 25
#         tw_layout[0,:] = LText(scene, "solitary waves")
        
        tightlimits!.([epileptic_ax])#, tw_ax])
        
        layout[1,1] = [ep_layout]#, tw_layout]
        
        layout[0,:] = LText(scene, "Proportion of simulations exhibiting traveling fronts")
        
        path = plotsdir("$(example_name)/$(sim_name)/single_slices/$(mod_names[x])_$(mod_names[y])_front_slice.png")
        mkpath(dirname(path))
        Makie.save(path, scene)
    end
end

# +
# One plot
using Makie, MakieLayout, DrWatson

let mod_names = string.(mod_names)
    plot_color = :magma
    all_dims = 1:length(mod_names)
    slices_2d = IterTools.subsets(all_dims, Val{2}())
    plot_side_size = 350 * (length(all_dims) - 1)
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
        
        @assert size(mean_is_epileptic) == length.((mod_values[x], mod_values[y]))
        
        epileptic_layout[x,y] = epileptic_ax = LAxis(epileptic_scene); 
        tw_layout[x,y] = tw_ax = LAxis(tw_scene)
        tightlimits!.([epileptic_ax, tw_ax])
        (Makie.heatmap!(epileptic_ax, my, mx, mean_is_epileptic', colorrange=(0,1)),
            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop epileptic"),
        Makie.heatmap!(tw_ax, my, mx, mean_is_traveling_solitary', colorrange=(0,0.100001)))
            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop traveling solitary")
    end
    ep_heatmaps, tw_heatmaps = zip(heatmap_pairs...)
    epileptic_layout[:,1] = LText.(epileptic_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
    epileptic_layout[end+1,2:end] = LText.(epileptic_scene, mod_names[2:end], tellwidth=false)
    epileptic_layout[0, :] = LText(epileptic_scene, "Traveling fronts", textsize = 30)
    ep_cbar = epileptic_layout[2:end-1, end+1] = LColorbar(epileptic_scene, ep_heatmaps[1], label = "Proportion traveling fronts")
    ep_cbar.width = 25
    ep_path = plotsdir("$(example_name)/$(sim_name)/epileptic_slices.png")
    mkpath(ep_path |> dirname)
    Makie.save(ep_path, epileptic_scene)
    
    tw_layout[:,1] = LText.(tw_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
    tw_layout[end+1,2:end] = LText.(tw_scene, mod_names[2:end], tellwidth=false)
    tw_layout[0, :] = LText(tw_scene, "Traveling waves", textsize = 30)
    tw_cbar = tw_layout[2:end-1, end+1] = LColorbar(tw_scene, tw_heatmaps[1], label = "Proportion traveling waves")
    tw_cbar.width = 25
    tw_path = plotsdir("$(example_name)/$(sim_name)/tw_slices.png")
    mkpath(tw_path |> dirname)
    Makie.save(tw_path, tw_scene)
end
@show mod_names
# -

# ## Search for modification 

# + jupyter={"source_hidden": true}
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

# ## Run specific modification

# + jupyter={"outputs_hidden": true, "source_hidden": true}
function animate_execution(filename, execution::AbstractFullExecution{T,<:Simulation{T}}; kwargs...) where T
    solution = execution.solution
    pop_names = execution.simulation.model.pop_names
    x = coordinate_axes(Simulation73.reduced_space(execution))[1]
    t = timepoints(execution)
    max_val = maximum(solution)
	min_val = minimum(solution)
    
    scene = Scene();
    time_idx_node = Node(1)
    single_pop = lift(idx -> population_timepoint(solution, 1, idx), time_idx_node)
    lines!(scene, x, single_pop)
    ylims!(scene, (min_val, max_val))
    
    @show t
    record(scene, filename, 1:length(t); framerate=20) do time_idx # TODO @views
        time_idx_node[] = time_idx
    end
end

# + jupyter={"source_hidden": true}
example_name = "reduced_line_dos_effectively_sigmoid"
line_example = get_example(example_name)
these_mods = (Aee=40.0,Aei=200.0, Aie=73.0, blocking_θE=25.0,blocking_θI=25.0,firing_θE=6.0,firing_θI=7.0, other_opts=Dict())
(mod_name, exec) = TravelingWaveSimulations.execute_single_modification(line_example, these_mods)
#exec = execute(line_example(;mods..., other_opts=Dict()))
wp = TravelingWaveSimulations.get_wave_properties(exec)
@show wp.epileptic
@show wp.traveling_solitary
path = plotsdir("$(example_name)/line_anim_$(mods_filename(mods)).mp4")
animate_execution(path, exec);

# +
# One plot with targets (no traveling wave)

let target_mods_nt = (Aee=137.0,Aei=115.0, Aie=132.0, Aii=30.0), 
        mod_names = string.(mod_names),
        line_example = get_example("reduced_line_dos_effectively_sigmoid")
    target_mods = [target_mods_nt[Symbol(name)] for name in mod_names]
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
        
        @assert size(mean_is_epileptic) == length.((mod_values[x], mod_values[y]))
        
        epileptic_layout[x,y] = epileptic_ax = LAxis(epileptic_scene); 
        tw_layout[x,y] = tw_ax = LAxis(tw_scene)
        tightlimits!.([epileptic_ax, tw_ax])
        
        heatmaps = (Makie.heatmap!(epileptic_ax, my, mx, mean_is_epileptic', colorrange=(0,1)),
            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop epileptic"),
        Makie.heatmap!(tw_ax, my, mx, mean_is_traveling_solitary', colorrange=(0,0.100001)))
            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop traveling solitary")
        
        Makie.scatter!(epileptic_ax, [target_mods[y]], [target_mods[x]], color=:red, markersize=5)
        Makie.scatter!(tw_ax, [target_mods[y]], [target_mods[x]], color=:red, markersize=5)
        
        heatmaps
    end
    
    savedir(bn) = plotsdir("$(example_name)/$(sim_name)/$(mods_filename(target_mods_nt))", bn)
    
    ep_heatmaps, tw_heatmaps = zip(heatmap_pairs...)
    epileptic_layout[:,1] = LText.(epileptic_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
    epileptic_layout[end+1,2:end] = LText.(epileptic_scene, mod_names[2:end], tellwidth=false)
    epileptic_layout[0, :] = LText(epileptic_scene, "Traveling fronts", textsize = 30)
    ep_cbar = epileptic_layout[2:end-1, end+1] = LColorbar(epileptic_scene, ep_heatmaps[1], label = "Proportion traveling fronts")
    ep_cbar.width = 25
    ep_path = savedir("epileptic_slices.png")
    mkpath(ep_path |> dirname)
    Makie.save(ep_path, epileptic_scene)
    
    tw_layout[:,1] = LText.(tw_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
    tw_layout[end+1,2:end] = LText.(tw_scene, mod_names[2:end], tellwidth=false)
    tw_layout[0, :] = LText(tw_scene, "Traveling waves", textsize = 30)
    tw_cbar = tw_layout[2:end-1, end+1] = LColorbar(tw_scene, tw_heatmaps[1], label = "Proportion traveling waves")
    tw_cbar.width = 25
    tw_path = savedir("tw_slices.png")
    mkpath(tw_path |> dirname)
    Makie.save(tw_path, tw_scene)
    
    these_mods =  (fixed_mods_dict..., target_mods_nt..., save_idxs=nothing, other_opts=Dict())
    @show these_mods
    (mod_name, exec) = TravelingWaveSimulations.execute_single_modification(line_example, these_mods)
    #exec = execute(line_example(;mods..., other_opts=Dict()))
    wp = TravelingWaveSimulations.get_wave_properties(exec)
    @show wp.epileptic
    @show wp.traveling_solitary
    anim_path = savedir("sim_animation.mp4")
    animate_execution(anim_path, exec);
    
    scene, layout = exec_heatmap_slices(exec, 5, (1400,1000))
    layout[0,:] = LText(scene, "Simulation: $(target_mods_nt)")#"\n Inhibition blocking threshold: $(these_mods[:blocking_θI])") 
    sim_heatmap_path = savedir("sim_heatmap.png")
    Makie.save( sim_heatmap_path, scene)
    
end

# +
# One plot with targets (traveling front)

let target_mods_nt = (Aee=200.0,Aei=75.0, Aie=65.0, Aii=4.0), 
        mod_names = string.(mod_names),
        line_example = get_example("reduced_line_dos_effectively_sigmoid")
    target_mods = [target_mods_nt[Symbol(name)] for name in mod_names]
    plot_color = :magma
    all_dims = 1:length(mod_names)
    slices_2d = IterTools.subsets(all_dims, Val{2}())
    plot_side_size = 600 * (length(all_dims) - 1)
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
        
        @assert size(mean_is_epileptic) == length.((mod_values[x], mod_values[y]))
        
        epileptic_layout[x,y] = epileptic_ax = LAxis(epileptic_scene); 
        tw_layout[x,y] = tw_ax = LAxis(tw_scene)
        tightlimits!.([epileptic_ax, tw_ax])
        
        heatmaps = (Makie.heatmap!(epileptic_ax, my, mx, mean_is_epileptic', colorrange=(0,1)),
            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop epileptic"),
        Makie.heatmap!(tw_ax, my, mx, mean_is_traveling_solitary', colorrange=(0,0.100001)))
            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop traveling solitary")
        
        Makie.scatter!(epileptic_ax, [target_mods[y]], [target_mods[x]], color=:red, markersize=5)
        Makie.scatter!(tw_ax, [target_mods[y]], [target_mods[x]], color=:red, markersize=5)
        
        heatmaps
    end
    
    savedir(bn) = plotsdir("$(example_name)/$(sim_name)/$(mods_filename(target_mods_nt))", bn)
    
    ep_heatmaps, tw_heatmaps = zip(heatmap_pairs...)
    epileptic_layout[:,1] = LText.(epileptic_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
    epileptic_layout[end+1,2:end] = LText.(epileptic_scene, mod_names[2:end], tellwidth=false)
    epileptic_layout[0, :] = LText(epileptic_scene, "Traveling fronts", textsize = 30)
    ep_cbar = epileptic_layout[2:end-1, end+1] = LColorbar(epileptic_scene, ep_heatmaps[1], label = "Proportion traveling fronts")
    ep_cbar.width = 25
    ep_path = savedir("epileptic_slices.png")
    mkpath(ep_path |> dirname)
    Makie.save(ep_path, epileptic_scene)
    
    tw_layout[:,1] = LText.(tw_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
    tw_layout[end+1,2:end] = LText.(tw_scene, mod_names[2:end], tellwidth=false)
    tw_layout[0, :] = LText(tw_scene, "Traveling waves", textsize = 30)
    tw_cbar = tw_layout[2:end-1, end+1] = LColorbar(tw_scene, tw_heatmaps[1], label = "Proportion traveling waves")
    tw_cbar.width = 25
    tw_path = savedir("tw_slices.png")
    mkpath(tw_path |> dirname)
    Makie.save(tw_path, tw_scene)
    
    these_mods =  (fixed_mods_dict..., target_mods_nt..., save_idxs=nothing, other_opts=Dict())
    @show these_mods
    (mod_name, exec) = TravelingWaveSimulations.execute_single_modification(line_example, these_mods)
    #exec = execute(line_example(;mods..., other_opts=Dict()))
    wp = TravelingWaveSimulations.get_wave_properties(exec)
    @show wp.epileptic
    @show wp.traveling_solitary
    anim_path = savedir("sim_animation.mp4")
    animate_execution(anim_path, exec);
    
    scene, layout = exec_heatmap_slices(exec)
    sim_heatmap_path = savedir("sim_heatmap.png")
    Makie.save( sim_heatmap_path, scene)
    
end

# +
# One plot with targets (traveling wave)

let target_mods_nt = (Aee=167.0,Aei=95.0, Aie=132.0, Aii=30.0), 
        mod_names = string.(mod_names),
        line_example = get_example("reduced_line_dos_effectively_sigmoid")
    target_mods = [target_mods_nt[Symbol(name)] for name in mod_names]
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
        
        @assert size(mean_is_epileptic) == length.((mod_values[x], mod_values[y]))
        
        epileptic_layout[x,y] = epileptic_ax = LAxis(epileptic_scene); 
        tw_layout[x,y] = tw_ax = LAxis(tw_scene)
        tightlimits!.([epileptic_ax, tw_ax])
        
        heatmaps = (Makie.heatmap!(epileptic_ax, my, mx, mean_is_epileptic', colorrange=(0,1)),
            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop epileptic"),
        Makie.heatmap!(tw_ax, my, mx, mean_is_traveling_solitary', colorrange=(0,0.100001)))
            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop traveling solitary")
        
        Makie.scatter!(epileptic_ax, [target_mods[y]], [target_mods[x]], color=:red, markersize=5)
        Makie.scatter!(tw_ax, [target_mods[y]], [target_mods[x]], color=:red, markersize=5)
        
        heatmaps
    end
    
    savedir(bn) = plotsdir("$(example_name)/$(sim_name)/$(mods_filename(target_mods_nt))", bn)
    
    ep_heatmaps, tw_heatmaps = zip(heatmap_pairs...)
    epileptic_layout[:,1] = LText.(epileptic_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
    epileptic_layout[end+1,2:end] = LText.(epileptic_scene, mod_names[2:end], tellwidth=false)
    epileptic_layout[0, :] = LText(epileptic_scene, "Traveling fronts", textsize = 30)
    ep_cbar = epileptic_layout[2:end-1, end+1] = LColorbar(epileptic_scene, ep_heatmaps[1], label = "Proportion traveling fronts")
    ep_cbar.width = 25
    ep_path = savedir("epileptic_slices.png")
    mkpath(ep_path |> dirname)
    Makie.save(ep_path, epileptic_scene)
    
    tw_layout[:,1] = LText.(tw_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
    tw_layout[end+1,2:end] = LText.(tw_scene, mod_names[2:end], tellwidth=false)
    tw_layout[0, :] = LText(tw_scene, "Traveling waves", textsize = 30)
    tw_cbar = tw_layout[2:end-1, end+1] = LColorbar(tw_scene, tw_heatmaps[1], label = "Proportion traveling waves")
    tw_cbar.width = 25
    tw_path = savedir("tw_slices.png")
    mkpath(tw_path |> dirname)
    Makie.save(tw_path, tw_scene)
    
    these_mods =  (fixed_mods_dict..., target_mods_nt..., save_idxs=nothing, other_opts=Dict())
    @show these_mods
    (mod_name, exec) = TravelingWaveSimulations.execute_single_modification(line_example, these_mods)
    #exec = execute(line_example(;mods..., other_opts=Dict()))
    wp = TravelingWaveSimulations.get_wave_properties(exec)
    @show wp.epileptic
    @show wp.traveling_solitary
    anim_path = savedir("sim_animation.mp4")
    animate_execution(anim_path, exec);
    
    scene, layout = exec_heatmap_slices(exec, 5, (1400,1000))
    layout[0,:] = LText(scene, "Simulation: $(target_mods_nt)") 
    sim_heatmap_path = savedir("sim_heatmap.png")
    Makie.save( sim_heatmap_path, scene)
    
end

# +
# One plot with targets (S)
using Makie, MakieLayout, DrWatson

let target_mods_nt = (See=60.0,Sei=110.0, Sie=90.0, Sii=100.0), 
        mod_names = string.(mod_names),
        line_example = get_example("reduced_line_dos_effectively_sigmoid")
    target_mods = [target_mods_nt[Symbol(name)] for name in mod_names]
    plot_color = :magma
    all_dims = 1:length(mod_names)
    slices_2d = IterTools.subsets(all_dims, Val{2}())
    plot_side_size = 600 * (length(all_dims) - 1)
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
        
        @assert size(mean_is_epileptic) == length.((mod_values[x], mod_values[y]))
        
        epileptic_layout[x,y] = epileptic_ax = LAxis(epileptic_scene); 
        tw_layout[x,y] = tw_ax = LAxis(tw_scene)
        tightlimits!.([epileptic_ax, tw_ax])
        
        heatmaps = (Makie.heatmap!(epileptic_ax, my, mx, mean_is_epileptic', colorrange=(0,1)),
            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop epileptic"),
        Makie.heatmap!(tw_ax, my, mx, mean_is_traveling_solitary', colorrange=(0,0.100001)))
            #xlab=mod_names[y], ylab=mod_names[x], color=plot_color, title="prop traveling solitary")
        
        Makie.scatter!(epileptic_ax, [target_mods[y]], [target_mods[x]], color=:red, markersize=5)
        Makie.scatter!(tw_ax, [target_mods[y]], [target_mods[x]], color=:red, markersize=5)
        
        heatmaps
    end
    
    savedir(bn) = plotsdir("$(example_name)/$(sim_name)/$(mods_filename(target_mods_nt))", bn)
    
    ep_heatmaps, tw_heatmaps = zip(heatmap_pairs...)
    epileptic_layout[:,1] = LText.(epileptic_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
    epileptic_layout[end+1,2:end] = LText.(epileptic_scene, mod_names[2:end], tellwidth=false)
    epileptic_layout[0, :] = LText(epileptic_scene, "Traveling fronts", textsize = 30)
    ep_cbar = epileptic_layout[2:end-1, end+1] = LColorbar(epileptic_scene, ep_heatmaps[1], label = "Proportion traveling fronts")
    ep_cbar.width = 25
    ep_path = savedir("epileptic_slices.png")
    mkpath(ep_path |> dirname)
    Makie.save(ep_path, epileptic_scene)
    
    tw_layout[:,1] = LText.(tw_scene, mod_names[1:end-1], tellheight=false, rotation=pi/2)
    tw_layout[end+1,2:end] = LText.(tw_scene, mod_names[2:end], tellwidth=false)
    tw_layout[0, :] = LText(tw_scene, "Traveling waves", textsize = 30)
    tw_cbar = tw_layout[2:end-1, end+1] = LColorbar(tw_scene, tw_heatmaps[1], label = "Proportion traveling waves")
    tw_cbar.width = 25
    tw_path = savedir("tw_slices.png")
    mkpath(tw_path |> dirname)
    Makie.save(tw_path, tw_scene)
    
    these_mods =  (fixed_mods_dict..., target_mods_nt..., save_idxs=nothing, other_opts=Dict())
    @show these_mods
    (mod_name, exec) = TravelingWaveSimulations.execute_single_modification(line_example, these_mods)
    #exec = execute(line_example(;mods..., other_opts=Dict()))
    wp = TravelingWaveSimulations.get_wave_properties(exec)
    @show wp.epileptic
    @show wp.traveling_solitary
    anim_path = savedir("sim_animation.mp4")
    animate_execution(anim_path, exec);
    
    scene, layout = exec_heatmap_slices(exec)
    sim_heatmap_path = savedir("sim_heatmap.png")
    Makie.save( sim_heatmap_path, scene)
    
end


# + jupyter={"outputs_hidden": true}
these_data = TravelingWaveSimulations.extract_data_namedtuple(exec)
# -

pfronts = TravelingWaveSimulations.persistent_fronts(TravelingWaveSimulations.all_fronts(exec), exec.solution.t)
anim = custom_animate(exec, pfronts)
mp4(anim, "wavefront_tmp/$(example_name)/fronts_line_anim_$(mods_filename(mods)).mp4")

TravelingWaveSimulations.get_velocities(pfronts[5])

TravelingWaveSimulations.based_on_example_NO_PARALLEL(; example_name="reduced_line_dos_effectively_sigmoid", 
    modifications=["blocking_θI=6.0:30.0", "blocking_θE=6.0:30.0"], 
    data_root=projectdir("tmp"))

based_on_example(; example_name="reduced_line_dos_effectively_sigmoid", modifications=["stim_strength=0.4:0.1:2.0", "stim_width=20.0:5.0:50.0", "Aee=22.0", "Aei=20.2"], max_batch_size=5000, data_root=projectdir("tmp"), max_sims_in_mem=30001, progress=true)




using FindPDE

example = get_example("reduced_line_dos_effectively_sigmoid")
