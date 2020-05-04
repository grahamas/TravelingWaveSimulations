# export Animate

# struct Animate <: AbstractPlotSpecification
#     fps::Int
#     output_name::String
#     kwargs::Dict
# end
# Animate(; fps=20, output_name="animation.mp4", kwargs...) = Animate(fps, output_name, kwargs)
# function analyse(plot_spec::Animate, execution::Execution, output_dir::AbstractString, prefix="")
#     path = full_name(output_name(plot_spec); path=output_dir, prefix=prefix)
#     recursively_clear_path(path)
#     anim = custom_animate(execution; plot_spec.kwargs...)
#     @show mp4(anim, path; fps=plot_spec.fps)
# end
# function custom_animate(execution::AbstractFullExecution{T,<:Simulation{T}}; kwargs...) where T
#     solution = execution.solution
#     pop_names = execution.simulation.model.pop_names
#     x = space(execution)
#     t = timepoints(execution)
#     max_val = maximum(solution)
# 	min_val = minimum(solution)
#     @animate for time_dx in 1:length(t) # TODO @views
#         plot([
#                 plot(
#                     x, population_timepoint(solution, i_pop, time_dx); label=pop_names[i_pop],
#                     val_lim=(min_val,max_val), title="t = $(round(t[time_dx], digits=4))",
#                     xlab = "Space (a.u. approx. um)",kwargs...
#                     )
#                 for i_pop in 1:length(pop_names)
#             ]...)
#     end
# end


# function TravelingWaveSimulations.custom_animate(execution::Execution{T,<:Simulation{T}}, fronts::AbstractArray{<:AbstractArray{<:TravelingWaveSimulations.AbstractWaveform}}; kwargs...) where T
#     solution = execution.solution
#     pop_names = execution.simulation.model.pop_names
#     x = space(execution)
#     x1 = [xcoord[1] for xcoord in x.arr]
#     t = timepoints(execution)
#     max_val = maximum(solution)
# 	min_val = minimum(solution)
#     i_pop = 1
#     @animate for time_dx in 1:length(t) # TODO @views
#         frame = population_timepoint(solution, i_pop, time_dx)
#         vs = ValuedSpace(frame, x1)
#         plot(
#             x, frame; label=pop_names[i_pop],
#             val_lim=(min_val,max_val), title="t = $(round(t[time_dx], digits=4))",
#             xlab = "Space (a.u. approx. um)",kwargs...
#             )
#         wf_arr = fronts[time_dx]
#         scatter!([wf.slope.loc for wf in wf_arr], [vs[wf.slope.loc].val for wf in wf_arr])
#     end
# end

# using Colors, AxisIndices
# function TravelingWaveSimulations.custom_animate(execution::Execution{T,<:Simulation{T}}, p_fronts::AbstractArray{<:Persistent{WF}}; kwargs...) where {T,WF}
#     ts = timepoints(execution)
#     if ts[end-1] == ts[end]
#         ts = ts[1:end-1]
#     end
#     n_p_fronts = length(p_fronts)
#     cmap = distinguishable_colors(n_p_fronts+1, RGB(1,1,1), dropseed=true)
#     fronts_by_time = AxisIndicesArray(Array{ID{WF},1}[ID{WF}[] for _ in ts], (ts,))
#     for (i_p_front, p_front) in enumerate(p_fronts)
#          for (wf, t) in zip(p_front.waveforms, p_front.t)
#             push!(fronts_by_time[t], ID{WF}(wf, i_p_front))
#         end
#     end
#     solution = execution.solution
#     pop_names = execution.simulation.model.pop_names
#     xs = space(execution)
#     x1 = [xcoord[1] for xcoord in xs.arr]
#     max_val = maximum(solution)
# 	min_val = minimum(solution)
#     i_pop = 1
    
#     @animate for time_dx in 1:length(ts) # TODO @views
#         frame = population_timepoint(solution, i_pop, time_dx)
#         vs = ValuedSpace(frame, x1)
#         plot(
#             xs, frame; label=pop_names[i_pop],
#             val_lim=(min_val,max_val), title="t = $(round(ts[time_dx], digits=4))",
#             xlab = "Space (a.u. approx. um)",kwargs...
#             )
#         wf_arr = fronts_by_time[time_dx]
#         if length(wf_arr) > 0
#             @show get_id.(wf_arr)
#             scatter!([wf.obj.slope.loc for wf in wf_arr], [vs[wf.obj.slope.loc].val for wf in wf_arr], color=[cmap[get_id(wf)] for wf in wf_arr])
#         end
#     end

# end