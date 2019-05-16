export SubsampledPlot

mutable struct SubsampledPlot <: AbstractPlotSpecification
    plot_type::Type{<:AbstractSpaceTimePlotSpecification}
    time_subsampler::Subsampler
    space_subsampler::Subsampler
    output_name::String
    kwargs::Iterators.Pairs
end
SubsampledPlot(; plot_type=nothing, time_subsampler=Subsampler(), space_subsampler=Subsampler(), output_name="", kwargs...) = SubsampledPlot(plot_type, time_subsampler, space_subsampler, output_name, kwargs)
@recipe function f(subsampledplot::SubsampledPlot, execution::Execution{T,<:Simulation{T,M}}) where {T,M<:WCMSpatial}
    simulation = execution.simulation
    t, x, wave = subsample(execution, time_subsampler=subsampledplot.time_subsampler, space_subsampler=subsampledplot.space_subsampler)

    dt = subsampledplot.time_subsampler.Δ == nothing ? save_dt(simulation) : subsampledplot.time_subsampler.Δ

    dx = subsampledplot.space_subsampler.Δ == nothing ? save_dx(simulation) : subsampledplot.space_subsampler.Δ

    plot_spec = subsampledplot.plot_type(;dt=dt, dx=dx, subsampledplot.kwargs...)
    if subsampledplot.output_name == ""
        subsampledplot.output_name = plot_spec.output_name
    end

    (plot_spec, t, x, wave)
end
