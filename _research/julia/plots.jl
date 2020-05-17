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

]up

# +
using DrWatson, Revise

using TravelingWaveSimulations, NeuralModels, Simulation73
using Makie, AbstractPlotting, MakieLayout, Colors

function autoscale!(sc::Scene)
    raw_scales = 1 ./ scene_limits(sc).widths
    scaling = raw_scales ./ minimum(raw_scales)
    scale!(sc, scaling)
end

AbstractPlotting.inline!(true)
# -

# # Plot nonlinearity

# +

TOL = 0.001

function _bounds_fns(nonl::SigmoidNonlinearity, tol=TOL)
    left = (xs) -> findfirst(x .> tol)
    right = (xs) -> findfirst(xs .> (1.0 - tol))
    return (left, right)
end

function _bounds_fns(nonl::DifferenceOfSigmoids, tol=TOL)
    left = (xs) -> findfirst(xs .> tol)
    right = (xs) -> findfirst(xs[left(xs):end] .< tol)
    return (left, right)
end

function _auto_range_idx(nonl::AbstractNonlinearity, output)
    (left_fn, right_fn) = _bounds_fns(nonl)
    left_idx, right_idx = (left_fn(output), right_fn(output))
    # TODO error messages
    left_idx = left_idx isa Int ? left_idx : 0
    right_idx = right_idx isa Int ? right_idx : length(output)
    return (left_idx:right_idx)
end

function _auto_range_idx(ps::AbstractArray, outputs::AbstractArray{<:AbstractArray})
    all_idxs = _auto_range_idx.(ps, outputs)
    firsts = [idx[1] for idx in all_idxs]
    lasts = [idx[end] for idx in all_idxs]
    return (minimum(firsts):maximum(lasts))
end

function _auto_range(nonl::AbstractNonlinearity, potential_range=collect(0.0:0.1:100.0))
    output = copy(potential_range)
    nonl(output, nothing, nothing)
    idxs = _auto_range_idx(nonl, output)
    (potential_range[idxs], output[idxs])
end
    
function _auto_range(nonls::AbstractArray{<:AbstractNonlinearity})
    max_range = collect(0.0:0.1:100.0)
    max_outputs = [nonl(copy(max_range), nothing, nothing) for nonl in nonls]
    idxs = _auto_range_idx(nonls, max_outputs)
    @show idxs
    ([max_range[idxs] for _ in max_outputs], [output[idxs] for output in max_outputs])
end    

function convert_arguments(nonl::AbstractNonlinearity)
    _auto_range(nonl)
end

function AbstractPlotting.convert_arguments(P::Type{<:AbstractPlot}, pops::Simulation73.AbstractPOneD{NPOPS,<:AbstractNonlinearity}) where NPOPS
    return _auto_range(Simulation73.array(pops))
end

# +
scene_sig = let line_example = get_example("reduced_line_dos_effectively_sigmoid")
    exec_sig = execute(line_example(; blocking_θ=[25.0,25.0], firing_θ = [6.0, 7.0], other_opts=Dict()));
    
    pop_names = exec_sig.simulation.model.pop_names
    nonl_xs, nonl_ys = _auto_range(exec_sig.simulation.model.nonlinearity |> Simulation73.array)
    sc, layout = layoutscene()
    ax = layout[1,1] = LAxis(sc, xlabel="input (a.u.)", ylabel="pop. activity (proportion)", title="Population activation functions")
    @show parse(Colorant, ax.attributes[:backgroundcolor][])
    colors = distinguishable_colors(length(pop_names), parse(Colorant, ax.attributes[:backgroundcolor][]), dropseed=true)
    plots = [plot!(ax, xs, ys, linestyle=:dash, linewidth=3, color=color) 
        for (xs, ys, color) in zip(nonl_xs, nonl_ys, colors)]
    legend_names = ["$(pop), eSig" for pop in pop_names]
    leg = LLegend(sc, plots, legend_names)
    layout[1,2] = leg

    exec_dos = execute(line_example(; blocking_θ=[25.0,10.0], firing_θ = [6.0, 7.0], other_opts=Dict()));

    nonl_xs, nonl_ys = _auto_range(exec_dos.simulation.model.nonlinearity |> Simulation73.array)
    append!(plots, [plot!(ax, xs, ys, color=color) for (xs, ys, color) in zip(nonl_xs, nonl_ys, colors)])
    append!(legend_names, ["$(pop), eDoS" for pop in pop_names])
    leg = LLegend(sc, plots, legend_names)
    layout[1,2] = leg
    
    sc
end


# -

Makie.save(plotsdir("sig_and_dos_nonlinearities.png"), scene_sig)

# # Connectivity plot

# +
using LaTeXStrings

function plot_connectivity_lims(simulation::Simulation73.Simulation, xlimits=nothing)
    scene, layout = layoutscene()
    ax = layout[1,1] = LAxis(scene, xlabel="distance (μm)", ylabel="connectivity strength (a.u.)",
        title="Connectivity kernels")
    
    pop_names = simulation.model.pop_names
    n_pops = length(pop_names)
    colors = distinguishable_colors(n_pops, parse(Colorant, ax.attributes[:backgroundcolor][]), dropseed=true)
    
    # Connectivity must be calculated before subsampling
    raw_space = simulation.space
    midpoint_idx = NeuralModels.fft_center_idx(raw_space)
    midpoint = coordinates(raw_space)[midpoint_idx]
    dists = differences(raw_space, midpoint) .|> x -> x[1]
    dists[CartesianIndex(1):midpoint_idx] .*= -1
    
    x_plot_idxs = if xlimits !== nothing
        findfirst(dists .>= xlimits[1]):findlast(dists .<= xlimits[2])
    else
        Colon()
    end
    onto_idx_idx = 1
    connectivity_array = Simulation73.array(simulation.model.connectivity)
    plots = map(pairs(IndexCartesian(), connectivity_array) |> collect) do (idx, conn)
        kern = NeuralModels.kernel(conn, raw_space)
        plot!(ax, dists[x_plot_idxs], kern[x_plot_idxs], color=colors[Tuple(idx)[onto_idx_idx]])
    end
    ontos = [Tuple(idx)[onto_idx_idx] for idx in CartesianIndices(connectivity_array)]
    legend_names = ["onto $pop" for pop in pop_names]
    leg = LLegend(scene, plots[1:n_pops], legend_names)
    layout[1,2] = leg
    
    scene
end
# -

scene_conn = let example_name = "reduced_line_dos_effectively_sigmoid"
    line_example = get_example(example_name)
    exec = execute(line_example(; blocking_θ=[25.0,25.0], firing_θ = [6.0, 7.0], other_opts=Dict()));

    sc = plot_connectivity_lims(exec.simulation, (-100., 100.))
end

Makie.save(plotsdir("gaussian_connectivity_example.png"), scene_conn)

# # Animation

function exec_animate(exec::AbstractExecution, output_file="test_anim.mp4"; framerate=25)
    @show size(exec.solution)
    scene = Scene()
    time_idx = Node(1)
    solution = exec.solution
    @show Simulation73.reduced_space(exec) |> size
    x = Simulation73.reduced_space(exec) |> coordinate_axes |> x -> x[1]
    @show size(x)
    @warn "only plotting one pop"
    activity = lift(tdx -> population(solution[tdx], 1), time_idx)
    time = lift(tdx -> solution.t[tdx], time_idx)
    scene = plot!(scene, x, activity)
    title(scene, "time: $time")
    xlabel!(scene, "space (μm)")
    ylabel!(scene, "activity (a.u.)")
    record(scene, output_file, 1:length(solution.t); framerate=framerate) do tdx
        time_idx[] = tdx
    end
end

scene_anim = let line_example = get_example("reduced_line_dos_effectively_sigmoid")
    exec_dos = execute(line_example(; blocking_θ=[25.0,10.0], firing_θ = [6.0, 7.0], other_opts=Dict()));
    sc = exec_animate(exec_dos)
end

# # Heatmap plot

# +
# function exec_heatmap(exec::AbstractExecution)
#     scene, layout = layoutscene(resolution=(1200, 1200))
#     soln = exec.solution
#     t = soln.t
#     xs = coordinate_axes(Simulation73.reduced_space(exec))[1] |> collect
#     pop_names = exec.simulation.model.pop_names

#     hm_axes = layout[1,1:length(pop_names)] = [LAxis(scene, title = "$pop_name activity") for pop_name in pop_names]
#     heatmaps = map(1:length(pop_names)) do idx_pop
#         ax = hm_axes[idx_pop]
#         pop_activity = cat(population.(soln.u, idx_pop)..., dims=2)
#         heatmap!(ax, t, xs, pop_activity')
#     end
#     tightlimits!.(hm_axes)
#     linkaxes!(hm_axes...)
#     hideydecorations!.(hm_axes[2:end])
#     cbar = layout[:, length(pop_names) + 1] = LColorbar(scene, heatmaps[1], label = "Activity Level")
#     cbar.width = 25
    
#     ylabel = layout[:,0] = LText(scene, "space (μm)", rotation=pi/2, tellheight=false)
#     xlabel = layout[end+1,2:3] = LText(scene, "time (ms)")
#     return (scene, layout)
# end 

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
        slices_layout[2,1:length(t_idxs)] = [LText(scene, "t=$(round(time, digits=1))", textsize=14, tellwidth=false) for time in t[t_idxs]]
        trim!(slices_layout)
        slices_layout        
    end
    layout[1,pop_idxs] = hm_axes
    cbar = layout[1, end+1] = LColorbar(scene, heatmaps[1], label = "Activity Level")
    cbar.width = 25

    ylabel = layout[:,0] = LText(scene, "space (μm)", rotation=pi/2, tellheight=false)
    xlabel = layout[end+1,2:3] = LText(scene, "time (ms)")
    return (scene, layout)
end

# -

scene_heatmap = let line_example = get_example("reduced_line_dos_effectively_sigmoid")
    exec_dos = execute(line_example(; blocking_θ=[25.0,10.0], firing_θ = [6.0, 7.0],save_idxs=nothing, other_opts=Dict()));
    (sc, layout) = exec_heatmap_slices(exec_dos)
    sc
end
scene_heatmap

save(plotsdir("heatmap_test.png"), scene_heatmap)

# # Waterfall plot

# +
# function default_theme(scene::SceneLike, ::Type{<: AbstractPlotting.Plot(AbstractExecution)})
#     Theme(
#         n_excerpts = 5
#     )
# end

# function AbstractPlotting.plot!(p::AbstractPlotting.Plot(AbstractExecution), n_excerpts=5)    @warn "Only plotting one pop"
#     exec = to_value(p[1])
#     n_excerpts = to_value(p[:n_excerpts])
#     soln = exec.solution
#     n_x, n_p, n_t = size(soln)
#     step = (length(soln.t) + 1) ÷ n_excerpts
#     ts = Array{Float64,1}(undef, n_x)
    
#     xs = coordinate_axes(space(exec))[1] |> collect
#     for idx in 1:step:length(soln.t)
#         ts .= soln.t[idx] |> collect
#         single_pop = soln[:,1,idx]
#         lines!(p, xs, ts, single_pop)
#     end
# end
# -

scene_hello = let line_example = get_example("reduced_line_dos_effectively_sigmoid")
    exec_dos = execute(line_example(; blocking_θ=[25.0,10.0], firing_θ = [6.0, 7.0], other_opts=Dict()));
    plot(exec_dos; n_excerpts=6)
    xlabel!(sc, "hello")
    autoscale!(sc)
    Makie.save("hello.png", sc)
    sc
end
