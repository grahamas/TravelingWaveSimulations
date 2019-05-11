struct ConnectivityPlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
ConnectivityPlot(; output_name = "connectivity.png", kwargs...) = ConnectivityPlot(output_name, kwargs)
# @recipe function f(plot_spec::ConnectivityPlot, execution::E) where {T,P,M<:WCMSpatial{T,1,P}, S<:Simulation{T,M}, E<:Execution{T,S}}
#     simulation = execution.simulation
#
#     pop_names = simulation.model.pop_names
#
#     connectivity = calculate.(simulation.model.connectivity, Ref(Calculated(simulation.model.space)))
#     map(Iterators.product(1:P, 1:P)) do (dst_pop, src_pop)
#         @series begin
#             lab --> "$(pop_names[src_pop]) → $(pop_names[dst_pop])"
#             seriestype := :heatmap
#             subplot := (dst_pop, src_pop)
#             connectivity[:,:,dst_pop,src_pop]
#         end
#     end
# end
@recipe function f(plot_spec::ConnectivityPlot, execution::E) where {T,P,M<:WCMSpatial{T,2,P}, S<:Simulation{T,M}, E<:Execution{T,S}}
    simulation = execution.simulation

    pop_names = simulation.model.pop_names

    connectivity = tensor(simulation.model.connectivity, simulation.model.space)
    layout := (2,2)
    for (dst_pop, src_pop) in Iterators.product(1:P, 1:P)
        @series begin
            lab --> "$(pop_names[src_pop])[1,1] → $(pop_names[dst_pop])"
            seriestype := :heatmap
            subplot := dst_pop + (src_pop-1) * P
            connectivity[:,:,1,1,dst_pop,src_pop]
        end
    end
end
