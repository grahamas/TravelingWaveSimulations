export radial_slice

struct AnimateRadialSlice <: AbstractPlotSpecification
    fps::Int
    output_name::String
    kwargs::Dict
end
AnimateRadialSlice(; fps=20, output_name="radial_slice_animation.mp4", kwargs...) = AnimateRadialSlice(fps, output_name, kwargs)
function analyse(plot_spec::AnimateRadialSlice, execution::Execution, output_dir, prefix="")
    path = full_name(output_name(plot_spec); path=output_dir, prefix=prefix)
    recursively_clear_path(path)
    anim = animate_radial_slice(execution; plot_spec.kwargs...)
    mp4(anim, path; fps=plot_spec.fps)
end

function animate_radial_slice(execution::Execution; kwargs...)
    o_dx = origin_dx(execution)
    # Slice from middle out to right
    slice_dxs = [o_dx[1], SliceToEnd(1,start=o_dx[2])]
    slice_coordinates = coordinates(execution)[slice_dxs...]
    t = time_points(execution)
    min_val, max_val = extrema(execution.solution)
    min_val = min(min_val, 0.0)
    pop_names = simulation.model.pop_names
    solution = execution.solution
    @animate for time_dx in 1:length(t)
        plot(x, pop_frame(solution, 1, time_dx)[slice_dxs...];
            label=pop_names[1], ylim=(min_val, max_val),
            title="t = $(round(t[time_dx], digits=4))",
            xlab = "Space (a.u. approx. um)", ylab = "Prop. of cells active",
            kwargs...)
        for i_pop in 2:length(pop_names)
            plot!(x, pop_frame(solution, i_pop, time_dx)[slice_dxs...]; label=pop_names[i_pop], kwargs...)
        end
    end
end
