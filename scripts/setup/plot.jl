
using Makie, StatsMakie
using TravelingWaveSimulationsPlotting
using Simulation73Plotting

Makie.convert_arguments(P::Makie.PointBased, nt::NamedTuple{names}) where names = Makie.convert_arguments(P, [names...], [values(nt)...])
