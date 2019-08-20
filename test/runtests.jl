using Pkg
Pkg.activate("docs/")
using Test, Documenter, TravelingWaveSimulations
doctest(TravelingWaveSimulations)
