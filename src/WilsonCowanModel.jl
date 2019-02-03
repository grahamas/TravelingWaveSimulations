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
#endregion

export WCMSpatial1D

export SigmoidNonlinearity, Sech2Nonlinearity, GaussianNonlinearity

export SharpBumpStimulus, NoisyStimulus, GaussianNoiseStimulus, NoStimulus

export ShollConnectivity

export MatchExample, StretchExample

include("nonlinearity.jl")
include("stimulus.jl")
include("connectivity.jl")
include("models.jl")
include("target.jl")
include("analysis.jl")

function run_example(example_name)
    include(example_name)
    filecopy(output, example_name, "parameters.jl")
end

export run_example

end #module
