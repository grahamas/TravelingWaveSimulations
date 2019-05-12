module TravelingWaveSimulations

using Simulation73, WilsonCowanModel
using MacroTools
import DifferentialEquations: Euler

include("examples.jl")
include("analysis.jl")

end #module
