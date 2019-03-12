module WilsonCowanModel

#region imports
using Parameters
using CalculatedTypes
import CalculatedTypes: Calculated, update
using RecipesBase
using StaticArrays
using Simulation73
import Simulation73: update_from_p!, make_calculated_function, target_loss
using JLD2
using Plots; pyplot()
using Random
using MacroTools, IterTools, Espresso
using Optim
using TensorOperations
#endregion

export WCMSpatial

export AbstractNonlinearity, SigmoidNonlinearity, Sech2Nonlinearity, GaussianNonlinearity

export AbstractStimulus, SharpBumpStimulus, NoisyStimulus, GaussianNoiseStimulus, NoStimulus

export AbstractConnectivity, ShollConnectivity, MeijerConnectivity

export MatchExample, StretchExample, SpatioTemporalFnTarget, @optim_st_target

include("helpers.jl")
include("nonlinearity.jl")
include("stimulus.jl")
include("connectivity.jl")
include("models.jl")
include("target.jl")
include("analysis.jl")

function run_simulation_example(example_name)
    include(example_name)
    filecopy(output, example_name, "parameters.jl")
    return simulation
end

function run_search_example(example_name)
    include(example_name)
    filecopy(output, example_name, "parameters.jl")
    return p_search
end

export run_simulation_example, run_search_example

end #module
