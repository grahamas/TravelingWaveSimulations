module TravelingWaveSimulations

using DrWatson, Lazy, Dates, BSON, Logging
using Simulation73, NeuralModels, WilsonCowanModel
using MacroTools
import DifferentialEquations: Euler
using Plots

export get_example
export plot_and_save

include("examples.jl")
include("analysis.jl")

include(joinpath(scriptdir(), "modifications.jl"))
include(joinpath(scriptdir(), "plot_specs.jl"))

include("script_helpers.jl")

end #module
