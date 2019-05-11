module WilsonCowanModel

using Parameters
using RecipesBase
using StaticArrays
using Simulation73
import Simulation73: target_loss
using Random
using MacroTools, IterTools, Espresso
using TensorOperations
using Statistics

export WCMSpatial

export AbstractNonlinearity, SigmoidNonlinearity, Sech2Nonlinearity, GaussianNonlinearity

export AbstractStimulus, SharpBumpStimulus, NoisyStimulus, GaussianNoiseStimulus, NoStimulus

export AbstractConnectivity, ShollConnectivity, MeijerConnectivity, GaussianConnectivity

include("helpers.jl")
include("nonlinearity.jl")
include("stimulus.jl")
include("connectivity.jl")
include("models.jl")

# using Optim
# export MatchExample, StretchExample, SpatioTemporalFnTarget, @optim_st_target
# include("target.jl")

end #module
