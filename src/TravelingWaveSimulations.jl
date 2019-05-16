module TravelingWaveSimulations

using Simulation73, WilsonCowanModel
using MacroTools
import DifferentialEquations: Euler
using Plots

export get_example
export plot_and_save

include("examples.jl")
include("analysis.jl")

end #module
