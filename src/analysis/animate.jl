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
function RecipesBase.animate(execution::Execution{T,<:Simulation{T,M}}; kwargs...) where {T,M<:WCMSpatial{T,1}}
    solution = execution.solution
    simulation = execution.simulation
    pop_names = simulation.model.pop_names
    x = saved_space_arr(simulation)
    t = saved_time_arr(simulation)
    max_val = maximum(solution)
    @animate for time_dx in 1:length(t) # TODO @views
        plot(x, pop_frame(solution, 1, time_dx); label=pop_names[1],
            ylim=(0,max_val), title="t = $(round(t[time_dx], digits=4))", kwargs...)
        for i_pop in 2:length(pop_names)
            plot!(x, pop_frame(solution, i_pop, time_dx); label=pop_names[i_pop], kwargs...)
        end
    end
end
function RecipesBase.animate(execution::Execution{T,<:Simulation{T,M}}; kwargs...) where {T,M<:WCMSpatial{T,2}}
    solution = execution.solution
    simulation = execution.simulation
    pop_names = simulation.model.pop_names
    x = saved_space_arr(simulation)
    xs = [coords[1] for coords in x[1,:]]
    ys = [coords[2] for coords in x[:,1]]
    t = saved_time_arr(simulation)
    max_val = maximum(solution)
    @animate for time_dx in 1:length(t) # TODO @views
        plot(
            [heatmap( pop_frame(solution, i_pop, time_dx);
                zlim=(0,max_val), clim=(0,max_val),
                title="t = $(round(t[time_dx], digits=4)), $(pop_names[i_pop])", kwargs...)
                for i_pop in 1:length(pop_names)]...
                    )
    end
end
