
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

function nonlinearity_plot!(scene::Scene, simulation::Simulation; save_dir=nothing)
    layout = GridLayout()

    pop_names = simulation.model.pop_names
    nonl_xs, nonl_ys = _auto_range(simulation.model.nonlinearity |> Simulation73.array)
    ax = layout[1,1] = LAxis(scene, xlabel="input (a.u.)", ylabel="pop. activity (proportion)", title="Population activation functions")
    @show parse(Colorant, ax.attributes[:backgroundcolor][])
    colors = distinguishable_colors(length(pop_names), parse(Colorant, ax.attributes[:backgroundcolor][]), dropseed=true)
    plots = [plot!(ax, xs, ys, linestyle=:dash, linewidth=3, color=color) 
        for (xs, ys, color) in zip(nonl_xs, nonl_ys, colors)]
    legend_names = ["$(pop)" for pop in pop_names]
    leg = LLegend(scene, plots, legend_names,
                  tellheight=false, tellwidth=false,
                  halign=:left, valign=:top, orientation=:vertical)
    layout[1,1] = leg

    layout
end


