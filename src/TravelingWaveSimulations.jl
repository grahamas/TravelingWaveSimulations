module TravelingWaveSimulations

using DrWatson, Lazy, Dates, BSON, Logging
using Simulation73, NeuralModels, WilsonCowanModel
using MacroTools
import DifferentialEquations: Euler
using Plots
using StaticArrays
using ArgParse

export get_example
export plot_and_save

include("examples.jl")
include("analysis.jl")

include(joinpath(scriptsdir(), "modifications.jl"))
include(joinpath(scriptsdir(), "plot_specs.jl"))

include("script_helpers.jl")

end #module
