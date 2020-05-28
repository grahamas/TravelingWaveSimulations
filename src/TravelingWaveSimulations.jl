module TravelingWaveSimulations
using DrWatson

using Distributed

using Lazy, Dates, BSON, Logging
using Simulation73, NeuralModels, WilsonCowanModel
using MacroTools
using DifferentialEquations
using ArgParse
using IterTools
using Statistics, LinearAlgebra, Distances, GLM
using RecipesBase

using JuliaDB, CSV, DataFrames

export execute_single_modification
export get_example
export plot_and_save
export custom_animate
export based_on_example
export load_directory
export get_wave_properties

#include("saving.jl")

include("valued_space.jl")

include("connectivity.jl")

include("callbacks.jl")
include("step_reductions.jl")
include("global_reductions.jl")

include("analysis.jl")
include("examples.jl")
include("loading.jl")

include("script_helpers/saving_utils.jl")
include("script_helpers.jl")
include(joinpath("post", "loading.jl"))

end #module
