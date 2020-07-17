module TravelingWaveSimulations
using DrWatson

using Distributed

using Lazy, Dates, BSON, Logging
using Simulation73, NeuralModels, WilsonCowanModel
using DifferentialEquations
using ArgParse
using IterTools
using Statistics, LinearAlgebra, Distances, GLM
using RecipesBase
using AxisIndices

using JuliaDB, CSV, DataFrames

include("util/valued_space.jl")

include("phenomena/fronts.jl")
include("phenomena/classifications.jl")
export ExecutionClassifications

include("callbacks.jl")
export terminate_when_E_fully_propagates
include("step_reductions.jl")
export front_array_type, reduce_to_fronts
include("global_reductions.jl")
export reduce_to_wave_properties

include("prototypes/prototypes.jl")
export get_prototype

include("util/loading.jl")
export DBRowIter, MultiDBRowIter, DBExecIter, MultiDBExecIter, MultiDB, 
       load_data_recent,
       load_ExecutionClassifications_recent
include("util/saving.jl")
include("plot/plotting.jl")
export custom_animate

include("util/parsing.jl")
export parse_modifications_argument, parse_analyses_argument
include("prototypes/iteration.jl")
export iterate_prototype, execute_single_modification
include("prototypes/replications.jl")

end #module
