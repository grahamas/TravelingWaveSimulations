export radial_slice
export AnimateRadialSlice

struct AnimateRadialSlice <: AbstractPlotSpecification
    fps::Int
    output_name::String
    kwargs::Dict
end
AnimateRadialSlice(; fps=20, output_name="radial_slice_animation.mp4", kwargs...) = AnimateRadialSlice(fps, output_name, kwargs)
function analyse(plot_spec::AnimateRadialSlice, execution::Execution, output_dir, prefix="")
    @warn "Slicing radially... and animating!..."
    path = full_name(output_name(plot_spec); path=output_dir, prefix=prefix)
    recursively_clear_path(path)
    @warn "saving to $path"
    anim = animate_radial_slice(execution; plot_spec.kwargs...)
    @warn "... done animating! Saving..."
    mp4(anim, path; fps=plot_spec.fps)
    @warn "... done saving!"
end

function animate_radial_slice(execution::Execution; kwargs...)
    o_dx = origin_idx(execution)
    # Slice from middle out to right
    slice_dxs = [StrideToEnd(1, o_dx[1]), o_dx[2]] # Straight down from center
    @show slice_dxs
    @show coordinates(execution) |> size
    slice_coordinates = [coord[1] for coord in coordinates(execution)[slice_dxs...]]
    t = timepoints(execution)
    min_val = minimum(execution.solution)
    max_val = maximum(execution.solution)
    min_val = min(min_val, 0.0)
    pop_names = execution.simulation.model.pop_names
    solution = execution.solution
    @show slice_coordinates |> size
    @show pop_frame(solution, 1, 1)[slice_dxs...] |> size
    @show slice_coordinates
    @animate for time_dx in 1:length(t)
        plot(slice_coordinates, pop_frame(solution, 1, time_dx)[slice_dxs...];
            label=pop_names[1], ylim=(min_val, max_val),
            title="t = $(round(t[time_dx], digits=4))",
            xlab = "Space (a.u. approx. um)", ylab = "Prop. of cells active",
            kwargs...)
        for i_pop in 2:length(pop_names)
            plot!(slice_coordinates, pop_frame(solution, i_pop, time_dx)[slice_dxs...]; label=pop_names[i_pop], kwargs...)
        end
    end
end
