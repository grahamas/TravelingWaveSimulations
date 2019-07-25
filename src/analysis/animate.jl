export Animate

struct Animate <: AbstractPlotSpecification
    fps::Int
    output_name::String
    kwargs::Dict
end
Animate(; fps=20, output_name="animation.mp4", kwargs...) = Animate(fps, output_name, kwargs)
function plot_and_save(plot_spec::Animate, execution::Execution, output_dir::AbstractString, prefix="")
    path = full_name(output_name(plot_spec); path=output_dir, prefix=prefix)
    recursively_clear_path(path)
    anim = RecipesBase.animate(execution; plot_spec.kwargs...)
    mp4(anim, path; fps=plot_spec.fps)
end
function RecipesBase.animate(execution::Execution{T,<:Simulation{T,M}}; kwargs...) where {T,M<:WCMSpatial{T}}
    solution = execution.solution
    simulation = execution.simulation
    pop_names = simulation.model.pop_names
    space = simulation.model.space
    @warn "not subsampling"
    t = saved_time_arr(simulation)
    max_val = maximum(solution)
    @animate for time_dx in 1:length(t) # TODO @views
        plot([
                plot(
                    space, pop_frame(solution, 1, time_dx); label=pop_names[1],
                    val_lim=(0,max_val), title="t = $(round(t[time_dx], digits=4))",
                    xlab = "Space (a.u. approx. um)", legend=:none, kwargs...
                    )
                for i_pop in 1:length(pop_names)
            ]...; size=(1600,800))
    end
end
