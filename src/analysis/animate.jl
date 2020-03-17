export Animate

struct Animate <: AbstractPlotSpecification
    fps::Int
    output_name::String
    kwargs::Dict
end
Animate(; fps=20, output_name="animation.mp4", kwargs...) = Animate(fps, output_name, kwargs)
function analyse(plot_spec::Animate, execution::Execution, output_dir::AbstractString, prefix="")
    path = full_name(output_name(plot_spec); path=output_dir, prefix=prefix)
    recursively_clear_path(path)
    anim = custom_animate(execution; plot_spec.kwargs...)
    @show mp4(anim, path; fps=plot_spec.fps)
end
function custom_animate(execution::AbstractFullExecution{T,<:Simulation{T}}; kwargs...) where T
    solution = execution.solution
    pop_names = execution.simulation.model.pop_names
    x = space(execution)
    t = timepoints(execution)
    max_val = maximum(solution)
	min_val = minimum(solution)
    @animate for time_dx in 1:length(t) # TODO @views
        plot([
                plot(
                    x, population_timepoint(solution, i_pop, time_dx); label=pop_names[i_pop],
                    val_lim=(min_val,max_val), title="t = $(round(t[time_dx], digits=4))",
                    xlab = "Space (a.u. approx. um)",kwargs...
                    )
                for i_pop in 1:length(pop_names)
            ]...)
    end
end
