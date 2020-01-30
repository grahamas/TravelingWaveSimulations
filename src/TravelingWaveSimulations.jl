module TravelingWaveSimulations
using DrWatson

using Distributed

using Lazy, Dates, BSON, Logging
using Simulation73, NeuralModels, WilsonCowanModel
using MacroTools
using LSODA
using DifferentialEquations
using Plots
using ArgParse
using IterTools
using Statistics, LinearAlgebra, Distances, GLM

using JuliaDB, CSV, DataFrames

export get_example
export plot_and_save
export custom_animate
export based_on_example
export load_directory

#include("saving.jl")

include("connectivity.jl")
include("examples.jl")
include("analysis.jl")

include("loading.jl")

include(joinpath(@__DIR__, "..", "scripts", "modifications.jl"))
include(joinpath(@__DIR__, "..", "scripts", "analyses.jl"))

include("script_helpers.jl")
include(joinpath("post", "loading.jl"))

end #module
