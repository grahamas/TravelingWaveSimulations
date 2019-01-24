module WCMTarget

using WCM
using Targets
using Parameters
using Modeling
using Simulating
import DifferentialEquations: ODESolution
using JLD2

abstract type AbstractExampleTarget{T} <: LossFunction{T} end

struct MatchExample{T} <: AbstractExampleTarget{T}
	data::Array{T}
	x::Array{T,1}
	t::Array{T,1}
end
function MatchExample(; file_name::String="")
	@load file_name wave x t
	MatchExample(wave, x, t)
end

function (p::MatchExample{T})(soln::ODESolution{T}, x_dxs, pop_dxs, t_dxs) where {T, AT<:Array{T,2}}
	sum((soln[x_dxs, pop_dxs, t_dxs] .- p.data) .^ 2)
end

function Targets.loss(fn::MatchExample{T}, model::WCMSpatial1D{T}, solver::Solver{T}) where T
	t_target, x_target = fn.t, fn.x
	t_dxs = subsampling_time_idxs(t_target, solver)
	x_dxs = subsampling_space_idxs(x_target, model, solver)
	pop_dxs = 1
	(soln) -> fn(soln, x_dxs, pop_dxs, t_dxs)
end
	
export MatchExample

struct StretchExample{T} <: AbstractExampleTarget{T}
	data::Array{T}
	x::Array{T,1}
	t::Array{T,1}
	stretch_dx::Int
end
function StretchExample(; file_name::String="", stretch_dx=nothing)
	@load file_name wave x t
	StretchExample(wave, x, t, stretch_dx)
end

function (p::StretchExample{T})(soln::ODESolution{T}, x_dxs, pop_dxs, t_dxs) where {T, AT<:Array{T,2}}
	sum((soln[x_dxs, pop_dxs, t_dxs] .- p.data) .^ 2)
end

function Targets.loss(fn::StretchExample{T}, model::WCMSpatial1D{T}, solver::Solver{T}) where T
	t_target, x_target = fn.t, fn.x
	t_dxs = subsampling_time_idxs(t_target, solver)
	x_dxs = subsampling_space_idxs(x_target, model, solver)
	pop_dxs = 1
	(soln) -> fn(soln, x_dxs, pop_dxs, t_dxs)
end
	
export StretchExample

@with_kw struct DecayingTraveling{T} <: LossFunction{T}
	timepoints::AbstractArray{T}
	target_pop::Int
	space_start::T
end

@with_kw struct Traveling{T} <: LossFunction{T}
	timepoints::AbstractArray{T}
	target_pop::Int
	space_start::T
end

function (p::DecayingTraveling{T})(soln::ODESolution, calc_space) where {T}
	time_vec = p.timepoints
	pop = p.target_pop
	timeseries = soln(time_vec)[calc_space[:,pop] .> p.space_start,pop,:]
	peak_vals = Array{T,1}(undef, length(time_vec))
	peak_dxs = Array{Int,1}(undef, length(time_vec))
	for time_dx in 1:length(time_vec)
		peak_vals[time_dx], peak_dxs[time_dx] = findmax(timeseries[:,time_dx])
	end

	return (penalize_decreasing(peak_dxs)
			+ must_travel(peak_vals, peak_dxs)
			+ must_not_increase(peak_vals)
			+ must_decrease(peak_vals)
			+ must_be_pulse(timeseries[1,:]))
end

end