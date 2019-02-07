
abstract type AbstractExampleTarget{T} <: AbstractTarget{T} end
abstract type AbstractFunctionTarget{T} <: AbstractTarget{T} end

#region Sech2Target

struct TravelingSech2Target{T} <: AbstractFunctionTarget{T}
	velocity::T
	width::T
	decay::T
end


#endregion

#region MatchExample

struct MatchExample{T} <: AbstractExampleTarget{T}
	data::Array{T}
	x::Array{T,1}
	t::Array{T,1}
end
function MatchExample(; file_name::String="")
	@load file_name wave x t
	MatchExample(wave, x, t)
end

function (p::MatchExample{T})(soln::AbstractArray{T}) where T#, x_dxs::AbstractArray{Int,1}, pop_dxs::AbstractArray{Int,1}, t_dxs::AbstractArray{Int,1}) where {T}
	res = sum((soln .- p.data) .^ 2)
	return res
end

function target_loss(fn::MatchExample{T}, model::WCMSpatial1D{T}, solver::Solver{T}) where T
	t_target, x_target = fn.t, fn.x
	t_dxs = subsampling_time_idxs(t_target, solver)
	x_dxs = subsampling_space_idxs(x_target, model, solver)
	pop_dxs = 1
	(soln) -> (fn(soln[x_dxs, pop_dxs, t_dxs]))
end

#endregion

#region StretchExample

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

function (p::StretchExample{T})(soln::AbstractArray{T}, x_dxs::AbstractArray{Int,1}, pop_dxs::AbstractArray{Int,1}, t_dxs::AbstractArray{Int,1}) where {T}
	sum((soln[x_dxs, pop_dxs, t_dxs] .- p.data) .^ 2)
end

function target_loss(fn::StretchExample{T}, model::WCMSpatial1D{T}, solver::Solver{T}) where T
	t_target, x_target = fn.t, fn.x
	t_dxs = subsampling_time_idxs(t_target, solver)
	x_dxs = subsampling_space_idxs(x_target, model, solver)
	pop_dxs = [1]
	(soln) -> fn(soln, x_dxs, pop_dxs, t_dxs)
end

#endregion
