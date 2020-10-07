module TravelingWaveSimulations
using DrWatson

using Distributed

using Lazy 
using Dates
using Logging
using Simulation73, NeuralModels, WilsonCowanModel
using DifferentialEquations, DiffEqOperators
using ArgParse
using IterTools
using Statistics, LinearAlgebra
using AxisIndices
using Interpolations

using JuliaDB

include("util/axisarray.jl")

include("phenomena/waveforms.jl")
include("phenomena/wavefronts.jl")
include("phenomena/persistent_waveforms.jl")
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

include("util/parsing.jl")
export parse_modifications_argument, parse_analyses_argument
include("util/io.jl")
include("util/saving.jl")
include("util/loading.jl")
export DBRowIter, MultiDBRowIter, DBExecIter, MultiDBExecIter, MultiDB, 
       load_simulation_data_recent, load_simulation_data, get_recent_simulation_data_path
       load_ExecutionClassifications_recent,
       load_ExecutionClassifications

include("prototypes/iteration.jl")
export iterate_prototype, execute_single_modification
include("prototypes/replications.jl")

end #module
