struct SpaceTimePlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
SpaceTimePlot(; output_name = "spacetime.png", kwargs...) = SpaceTimePlot(output_name, kwargs)
@recipe function f(plot_spec::SpaceTimePlot, execution::Execution{T,<:Simulation{T,M}}) where {T,M<:WCMSpatial}
    simulation = execution.simulation
    solution = execution.solution
    v_space = saved_space_arr(simulation)
    v_time = saved_time_arr(simulation)
    clims := (minimum(solution), maximum(solution))
    grid := false
    layout := (2,1)
    for i_pop in 1:2 # TODO!!
        @series begin
            seriestype --> :heatmap
            subplot := i_pop
            x := v_time
            y := v_space
            solution[:,i_pop,:]
        end
    end
end
