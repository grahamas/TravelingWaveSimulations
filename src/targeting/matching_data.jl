
abstract type AbstractDataTarget{T} <: AbstractTarget{T} end

struct MatchData{T} <: AbstractDataTarget{T}
	data::Array{T}
	x::Array{T,1}
	t::Array{T,1}
end
function MatchData(; file_name::String="")
	@load file_name wave x t
	MatchData(wave, x, t)
end

function (p::MatchData{T})(soln::AbstractArray{T}) where T#, x_dxs::AbstractArray{Int,1}, pop_dxs::AbstractArray{Int,1}, t_dxs::AbstractArray{Int,1}) where {T}
	res = sumsqdiff(soln, p.data)
	return res
end

function target_loss(fn::MatchData{T}, model::WCMSpatial{T}, solver::Solver{T}) where T
	t_target, x_target = fn.t, fn.x
	t_dxs = subsampling_time_idxs(solver, t_target)
	x_dxs = subsampling_space_idxs(model, solver, x_target)
	pop_dxs = 1
	loss_fn(soln) = fn(soln[x_dxs, pop_dxs, t_dxs])
	return loss_fn
end
