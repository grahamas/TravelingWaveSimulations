
function sumsqdiff(a::AbstractArray,b::AbstractArray)
	sum((a .- b) .^ 2)
end

include("targeting/optimizing_functions.jl")

@function_target sech2(x,t; amplitude, width, velocity) = amplitude * sech((1/width) * (x - velocity * t)) ^ 2
@function_target function sech2_state_change(x,t; bump_amp, bump_wid, bump_vel, state_inc)
	bump = bump_amp * sech((1/bump_wid) * (x - bump_vel * t)) ^ 2
	bump -= state_inc * tanh((1/bump_wid) * (x - bump_vel * t))
	return bump
end
@function_target function decaying_sech2(x,t; amplitude, velocity, decay, n)
	A = amplitude * exp(-decay * t)
	B = n / sqrt(2 * velocity * (n + 1) * (n + 2))
	A * sech(B * (x - (velocity * t)))^(2/n) 
end

include("targeting/matching_data.jl")
