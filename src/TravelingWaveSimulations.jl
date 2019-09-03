module TravelingWaveSimulations
using DrWatson

using Distributed

using Lazy, Dates, BSON, Logging
using Simulation73, NeuralModels, WilsonCowanModel
using MacroTools
import DifferentialEquations: Euler
using Plots
using StaticArrays
using ArgParse
using IterTools

using JuliaDB, CSV

export get_example
export plot_and_save
export custom_animate
export based_on_example

#include("saving.jl")

include("connectivity.jl")
include("examples.jl")
include("analysis.jl")

include(joinpath(scriptsdir(), "modifications.jl"))
include(joinpath(scriptsdir(), "analyses.jl"))

include("script_helpers.jl")

end #module
