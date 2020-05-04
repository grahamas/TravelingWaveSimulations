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

@recipe function f(connectivity::AbstractConnectivityParameter{T,N_CDT},
                   lattice::AbstractLattice{T,N_ARR,N_CDT};
                   source_location = Tuple(zero(T) for _ in 1:N_CDT)) where {T,N_ARR,N_CDT}
    weights = directed_weights(connectivity, lattice, source_location)
    (lattice, weights)
end

@recipe function f(plot_spec::ConnectivityPlot, execution::Execution) where {T}
    model = execution.simulation.model
    pop_names = model.pop_names
    connectivity = model.connectivity
    x = space(execution)

    (plot_spec, connectivity, x, pop_names)
end

