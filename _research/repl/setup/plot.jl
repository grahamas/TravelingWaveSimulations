
using Makie, MakieLayout, AbstractPlotting

AbstractPlotting.convert_arguments(P::AbstractPlotting.PointBased, nt::NamedTuple{names}) where names = AbstractPlotting.convert_arguments(P, [names...], [values(nt)...])
