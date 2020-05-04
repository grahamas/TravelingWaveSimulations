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

# +
using DrWatson

using TravelingWaveSimulations, NeuralModels, Simulation73, Makie, AbstractPlotting
using MakieLayout, Colors

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
example_name = "reduced_line_dos_effectively_sigmoid"
line_example = get_example(example_name)
exec_sig = execute(line_example(; blocking_θ=[25.0,25.0], firing_θ = [6.0, 7.0], other_opts=Dict()));

pop_names = this_exec.simulation.model.pop_names
this_exec = exec_sig
nonl_xs, nonl_ys = _auto_range(this_exec.simulation.model.nonlinearity |> Simulation73.array)
sc, layout = layoutscene()
ax = layout[1,1] = LAxis(sc, xlabel="input (a.u.)", ylabel="pop. activity (proportion)", title="Population activation functions")
@show parse(Colorant, ax.attributes[:backgroundcolor][])
colors = distinguishable_colors(length(pop_names), parse(Colorant, ax.attributes[:backgroundcolor][]), dropseed=true)
plots = [plot!(ax, xs, ys, linestyle=:dash, linewidth=3, color=color) 
    for (xs, ys, color) in zip(nonl_xs, nonl_ys, colors)]
legend_names = ["$(pop), eSig" for pop in pop_names]
leg = LLegend(sc, plots, legend_names)
layout[1,2] = leg

example_name = "reduced_line_dos_effectively_sigmoid"
line_example = get_example(example_name)
exec_dos = execute(line_example(; blocking_θ=[25.0,10.0], firing_θ = [6.0, 7.0], other_opts=Dict()));

nonl_xs, nonl_ys = _auto_range(exec_dos.simulation.model.nonlinearity |> Simulation73.array)
append!(plots, [plot!(ax, xs, ys, color=color) for (xs, ys, color) in zip(nonl_xs, nonl_ys, colors)])
append!(legend_names, ["$(pop), eDoS" for pop in pop_names])
leg = LLegend(sc, plots, legend_names)
layout[1,2] = leg

sc
# -

Makie.save(plotsdir("sig_and_dos_nonlinearities.png"), sc)

# # Connectivity plot

# +
using LaTeXStrings

example_name = "reduced_line_dos_effectively_sigmoid"
line_example = get_example(example_name)
exec = execute(line_example(; blocking_θ=[25.0,25.0], firing_θ = [6.0, 7.0], other_opts=Dict()));

function plot_connectivity(simulation::Simulation73.Simulation)
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
    onto_idx_idx = 1
    connectivity_array = Simulation73.array(simulation.model.connectivity)
    plots = map(pairs(IndexCartesian(), connectivity_array) |> collect) do (idx, conn)
        kern = NeuralModels.kernel(conn, raw_space)
        plot!(ax, dists, kern, color=colors[Tuple(idx)[onto_idx_idx]])
    end
    ontos = [Tuple(idx)[onto_idx_idx] for idx in CartesianIndices(connectivity_array)]
    legend_names = ["onto $pop" for pop in pop_names]
    leg = LLegend(scene, plots[1:n_pops], legend_names)
    layout[1,2] = leg
    
    scene
end


conn_scene = plot_connectivity(exec.simulation)
# -

Makie.save(plotsdir("gaussian_connectivity_example.png"), conn_scene)

# # Waterfall plot

# +
function default_theme(scene::SceneLike, ::Type{<: AbstractPlotting.Plot(AbstractExecution)})
    Theme(
        n_excerpts = 5
    )
end

function AbstractPlotting.plot!(p::AbstractPlotting.Plot(AbstractExecution), n_excerpts=5)    @warn "Only plotting one pop"
    exec = to_value(p[1])
    n_excerpts = to_value(p[:n_excerpts])
    soln = exec.solution
    n_x, n_p, n_t = size(soln)
    step = (length(soln.t) + 1) ÷ n_excerpts
    ts = Array{Float64,1}(undef, n_x)
    
    xs = coordinate_axes(space(exec))[1] |> collect
    for idx in 1:step:length(soln.t)
        ts .= soln.t[idx] |> collect
        single_pop = soln[:,1,idx]
        lines!(p, xs, ts, single_pop)
    end
end
# -

sc = plot(exec_dos; n_excerpts=6)
xlabel!(sc, "hello")
auto_scale!(sc)
Makie.save("hello.png", sc)
sc

function auto_scale!(sc::Scene)
    raw_scales = 1 ./ scene_limits(sc).widths
    scaling = raw_scales ./ minimum(raw_scales)
    scale!(sc, scaling)
end
auto_scale!(sc)
sc

sc = surface(1.0:0.01:2.0, 100.0:200.0, rand(10,10))
auto_scale!(sc)
sc
