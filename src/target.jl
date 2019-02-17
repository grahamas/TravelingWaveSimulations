
function sumsqdiff(a::AbstractArray,b::AbstractArray)
	sum((a .- b) .^ 2)
end

include("targeting/optimizing_functions.jl")

@optim_st_target sech2(x,t; amplitude, steepness, velocity) = amplitude * sech((1/width) * (x - velocity * t)) ^ 2
@optim_st_target function sech2_state_change(x,t; bump_amp, bump_wid, bump_vel, state_inc)
	bump = bump_amp * sech((1/bump_wid) * (x - bump_vel * t)) ^ 2
	bump -= state_inc * tanh((1/bump_wid) * (x - bump_vel * t))
	return bump
end

include("targeting/matching_data.jl")
