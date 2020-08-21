
using Makie, MakieLayout, AbstractPlotting, StatsMakie
using TravelingWaveSimulations.Plots
using Simulation73Plotting

AbstractPlotting.convert_arguments(P::AbstractPlotting.PointBased, nt::NamedTuple{names}) where names = AbstractPlotting.convert_arguments(P, [names...], [values(nt)...])
