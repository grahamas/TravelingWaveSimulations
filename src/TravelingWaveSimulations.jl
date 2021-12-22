module TravelingWaveSimulations
using DrWatson

using Distributed

using Lazy 
using Dates
using Logging
using Simulation73, NeuralModels, WilsonCowanModel
using DiffEqBase, DiffEqOperators, OrdinaryDiffEq
using ArgParse
using IterTools
using Statistics, LinearAlgebra
using Interpolations
using FileIO, BSON, TypedTables
# using JuliaDB

include("callbacks.jl")
export terminate_when_E_fully_propagates#, is_propagated

include("prototypes/prototypes.jl")
export get_prototype

include("util/parsing.jl")
export parse_modifications_argument, parse_analyses_argument
include("util/io.jl")
export parse_prototype_name_from_mdb_path, write_modifications!
include("util/saving.jl")
include("util/loading.jl")
export load_simulation_data_recent, load_simulation_data, get_recent_simulation_data_path
       load_classifications, load_classifications_recent

include("prototypes/iteration.jl")
export iterate_prototype, execute_single_modification
include("prototypes/replications.jl")

end #module
