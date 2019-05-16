export NonlinearityPlot

struct NonlinearityPlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
NonlinearityPlot(; output_name = "nonlinearity.png", kwargs...) = NonlinearityPlot(output_name, kwargs)
@recipe function f(plot_spec::NonlinearityPlot, execution::Execution{T,<:Simulation{T,M}}; resolution=100, fn_bounds=(-1.0,15.0)) where {T,M<:WCMSpatial}
    simulation = execution.simulation

    pop_names = simulation.model.pop_names
    n_pops = length(pop_names)

    nonlinearity_mutator! = make_mutator(simulation.model.nonlinearity)

    xlab := "Input current"
    ylab := "Proportion pop. reaching at least threshold"

    one_pop_x = range(fn_bounds[1], stop=fn_bounds[2], length=resolution)
    output = repeat(one_pop_x, outer=(1,n_pops))
    nonlinearity_mutator!(output, output, 0.0)

    for i_pop in 1:length(pop_names)
        @series begin
            lab --> pop_names[i_pop]
            seriestype := :line
            x := one_pop_x
            y := output[:,i_pop]
            ()
        end
    end
end
