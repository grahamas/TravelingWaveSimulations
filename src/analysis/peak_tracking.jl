struct PeakTravelingWavePlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
PeakTravelingWavePlot(; output_name="peak_traveling_wave.png", kwargs...) = PeakTravelingWavePlot(output_name, kwargs)
@recipe function f(plot_spec::PeakTravelingWavePlot, t::AbstractArray{T,1}, x::AbstractArray{T,1}, wave::AbstractArray{T,2}, interpolation_n::Int=1) where {T}
    title := "Traveling wave with labeled peaks"
    for time_dx in 1:length(t)
        @series begin
            seriestype := :line
            x := x
            y := wave[:,time_dx] # Subtract inhibitory activity...
            ()
        end
    end
    xs, peaks = track_wave_peak(x, wave, interpolation_n)
    @series begin
        seriestype := :scatter
        x := xs
        y := peaks
        ()
    end
end

function wave_maxima(single_wave_data::SPACE1DTIME) where {T, SPACE1DTIME<:AbstractArray{T,2}}
	# [space, time]
	max_vals, max_ixs = findmax(single_wave_data, dims=1)
	return (max_vals, max_ixs)
end

using LsqFit, OffsetArrays
function interpolate_parabola(space::OffsetArray{T,1}, wave::OffsetArray{T,1}) where T
	@. parabola(x,p) = p[1] + p[2] * ((x - p[3]) ^ 2)
	ub = [Inf, 0.0, space[1]]
	lb = [minimum(wave), -Inf, space[-1]]
	guess = [maximum(wave), -1.0, space[0]]
	fit = curve_fit(parabola, parent(space), parent(wave), guess, upper=ub, lower=lb)
    return (coef(fit)[3], coef(fit)[1])
end

function track_wave_peak(x::SPACE1D, wave::SPACE1DTIME, side::Int) where {T, SPACE1D<:AbstractArray{T,1}, SPACE1DTIME<:AbstractArray{T,2}}
	max_vals, max_ixs = wave_maxima(wave)
	space_ixs = [ix[1] for ix in max_ixs]
	circa_space_ixs = [(ix-side):(ix+side) for ix in space_ixs]
	circa_wave_ixs = [(ix-CartesianIndex(side,0)):(ix+CartesianIndex(side,0)) for ix in max_ixs]
	interpolated = map(zip(circa_space_ixs, circa_wave_ixs)) do (circa_space_ix, circa_wave_ix)
		if any(circa_space_ix .< 1) || any(circa_space_ix .> length(x))
            return (NaN, NaN)
        end
        circa_space = OffsetArray(x[circa_space_ix], -side:side)
		circa_wave = OffsetArray(dropdims(wave[circa_wave_ix], dims=2), -side:side)
		interpolate_parabola(circa_space, circa_wave)
    end
	return zip(interpolated...) .|> collect
end
