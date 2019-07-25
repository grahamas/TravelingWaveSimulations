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

@recipe function f(connectivity::AbstractConnectivity{T,N_CDT},
                   lattice::AbstractLattice{T,N_ARR,N_CDT};
                   source_location = Tuple(zero(T) for _ in 1:N_CDT)) where {T,N_ARR,N_CDT}
    weights = directed_weights(connectivity, lattice, source_location)
    (lattice, weights)
end

@recipe function f(plot_spec::ConnectivityPlot, execution::E) where {T,
                            M<:WCMSpatial,
                            S<:Simulation{T,M},
                            E<:Execution{T,S}}
    model = execution.simulation.model
    pop_names = model.pop_names
    space = model.space
    connectivity = model.connectivity

    (plot_spec, connectivity, space, pop_names)
end

using Plots: GridLayout
@recipe function f(plot_spec::ConnectivityPlot, connectivity::C, space::S,
        pop_names::SVector{P,<:AbstractString}) where {
            T, P,
            C<:SMatrix{P,P},
            S<:AbstractSpace{T}
        }
    layout := GridLayout(P,P)
    for (dst_pop, src_pop) in Iterators.product(1:P, 1:P)
        @series begin
            lab --> "$(pop_names[src_pop])[1,1] → $(pop_names[dst_pop])"
            subplot := (dst_pop, src_pop)
            (connectivity[dst_pop, src_pop], space)
        end
    end
end

@recipe function f(plot_spec::ConnectivityPlot, connectivity::C, space::RandomlyEmbeddedLattice{T,N_ARR,N_CDT},
        pop_names::SVector{P,<:AbstractString}; source_location = Tuple(zero(T) for _ in 1:N_CDT)) where {
            T, N_ARR, N_CDT, P,
            C<:SMatrix{P,P}
        }
    @assert P == 2
    layout := @layout [[a; b] [c; d]; [e; f] [g; h]]
    size := (1600,1600)
    subplot = 1
    for (dst_pop, src_pop) in Iterators.product(1:P, 1:P)
        weights = directed_weights(connectivity[dst_pop, src_pop], space, source_location)
        @series begin
            title := "$(pop_names[src_pop]) → $(pop_names[dst_pop])"
            subplot := subplot
            (space.lattice, weights)
        end
        @series begin
            title := "$(pop_names[src_pop]) → $(pop_names[dst_pop])"
            subplot := subplot+1
            (space.embedded_lattice, unembed_values(space, weights))
        end
        subplot += 2
    end
end
