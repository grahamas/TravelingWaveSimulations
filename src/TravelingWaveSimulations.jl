module TravelingWaveSimulations

using DrWatson, Lazy, Dates, BSON, Logging
using Simulation73, NeuralModels, WilsonCowanModel
using MacroTools
import DifferentialEquations: Euler
using Plots
using StaticArrays

export get_example
export plot_and_save
export custom_animate
export based_on_example

include("connectivity.jl")
include("examples.jl")
include("analysis.jl")

include(joinpath(scriptsdir(), "modifications.jl"))
include(joinpath(scriptsdir(), "analyses.jl"))

include("script_helpers.jl")

end #module
