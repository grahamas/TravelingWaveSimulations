module Simulation73

using DrWatson
using Markdown # for doc_str
using DifferentialEquations, DiffEqBase#, DiffEqParamEstim
#using BlackBoxOptim, Optim
using JLD2
import DifferentialEquations: DESolution, OrdinaryDiffEqAlgorithm, solve, Euler, ODEProblem
using OrdinaryDiffEq
#using RecipesBase
using Parameters
using Lazy
using AxisIndices
using AbstractPlotting

# ENV["GKSwstype"] = "100" # For headless plotting (on server)
# ENV["MPLBACKEND"]="Agg"
# using Plots

# "variables.jl"
export AbstractVariable, UnboundedVariable, BoundedVariable,
	default_value, bounds, pops, MaybeVariable,
	AbstractParameter, AbstractAction, AbstractSpaceAction


# space.jl
export AbstractSpace, AbstractLattice, AbstractPeriodicLattice, AbstractCompactLattice,
	AbstractEmbeddedLattice

export CompactLattice, PeriodicLattice

export RandomlyEmbeddedLattice, unembed_values

export coordinates, origin_idx, differences, coordinate_axes, timepoints, space,
    extent, abs_difference, abs_difference_periodic, discrete_segment, fft_center_dx

# "action.jl"
export AbstractAction, AbstractInteraction, AbstractSpaceAction, AbstractSpaceInteraction
export NullifyParameter, NullAction

# "population.jl"
export AbstractPopulationActionsParameters, AbstractPopulationInteractionsParameters,
    AbstractPopulationActions, AbstractPopulationInteractions,
    population, population_coordinates, population_timepoint, population_repeat,
    PopulationActionsParameters, PopulationInteractionsParameters,
    PopulationActions2, PopActParam, PopInteractParam, PopAct, PopInteract

# "subsampling.jl" (note: should probably be meshed with meshes)
export scalar_to_idx_window, subsampling_Δidx, subsampling_idxs,
	subsampling_time_idxs, subsampling_space_idxs, AbstractSubsampler,
    IndexSubsampler, ValueSubsampler, ValueWindower, RadialSlice, StrideToEnd

# "stimulus.jl"
export AbstractStimulusParameter, AbstractStimulusAction,
    NoStimulusParameter, NoStimulusAction

# "analysing.jl"
export AbstractPlotSpecification, AbstractSpaceTimePlotSpecification, Analyses,
	output_name, plot_and_save, analyse, subsample, subsampling_idxs

# "targets.jl"
export AbstractTarget, target_loss

export execute

# "simulating.jl"
export AbstractModel, AbstractModelwithDelay, Solver, Simulation, 
    AbstractExecution, Execution, FailedExecution,
	initial_value, history, time_span, saved_dt, saved_dx,
	generate_problem, solve, run_simulation,
	make_mutators, make_system_mutator,
    BareSolution, reduced_space, frame_xs

# # "exploring.jl"
# export Search, SearchExecution, make_problem_generator, search, run_search

abstract type AbstractParameter{T} end
include("helpers.jl")
include("deconstructing.jl")
include("variables.jl")
include("space.jl")
include("parameter.jl")
include("action.jl")
include("population.jl")
include("subsampling.jl") # depends on space.jl
include("stimulus.jl") # depends on variables.jl, population.jl, and space.jl
include("solutions.jl")
#include("solvers.jl")
include("simulating.jl")
include("targets.jl")
# include("exploring.jl")
include("analysing.jl")
include("plotting.jl")
export heatmap_slices_execution, animate_execution,
    exec_heatmap, exec_heatmap!,
    exec_heatmap_slices

end
