export WaveStatsPlot, WaveWidthPlot, WaveVelocityPlot, NeumanTravelingWavePlot

function calculate_width(single_wave_data::AT) where {T, AT<:AbstractArray{T,2}}
    # [space, time]
    space_max = findmax(single_wave_data)[1]
    above_half = single_wave_data .> (space_max / 2)
    frame_first(frame)::Union{Int64,Missing} = any(frame) ? findfirst(frame)[1] : 0
    frame_last(frame)::Union{Int64,Missing} = any(frame) ? findlast(frame)[1] : 0
    half_start = mapslices(frame_first, above_half, dims=1)
    half_end = mapslices(frame_last, above_half, dims=1)
    return half_end - half_start
end

function calculate_wave_velocity(x::SPACE1D, wave::SPACE1DTIME, dt::T=one(T), interpolation_n::Int=1) where {T, SPACE1D<:AbstractArray{T,1}, SPACE1DTIME<:AbstractArray{T,2}}
    # [space, time]
	xs, vals = track_wave_peak(x, wave, interpolation_n)
	diff(xs) ./ dt
end

function smooth(t::AbstractArray, arr::AbstractArray, window)
    window = round(Int, window)
    @assert window < length(arr) / 2
    output_len::Int = length(arr) - window
    output_arr = [mean(arr[dx:dx+window]) for dx in 1:output_len]
    return t[1:output_len], output_arr
end


struct WaveVelocityPlot{T} <: AbstractPlotSpecification
    output_name::String
    dt::T
    dx::T
    smoothing::Union{T,Nothing}
	interpolation_n::Int
    kwargs::Dict
end
WaveVelocityPlot(; output_name="wave_velocity.png", dt::Union{Nothing,T}=nothing, dx::Union{Nothing,T}=nothing, smoothing::Union{Nothing,T}=nothing, interpolation_n=1, kwargs...) where {T<:Float64} = WaveVelocityPlot{T}(output_name, dt, dx, smoothing, interpolation_n, kwargs)
@recipe function f(plot_spec::WaveVelocityPlot{T}, t::AbstractArray{T,1}, x::AbstractArray{T,1}, wave::AbstractArray{T,2}, transform::Function=identity, naive=false) where {T}
    velocity = calculate_wave_velocity(x, wave, plot_spec.dt, plot_spec.interpolation_n)
    if plot_spec.smoothing != nothing
        t, velocity = smooth(t, velocity, plot_spec.smoothing/plot_spec.dt)
    end
	xlab := "Time (a.u. approx. ms)"
	ylab := "Velocity (a.u. approx. um/ms)"
    @series begin
        title --> "Wave velocity over time"
        seriestype := :scatter
        ylims := (0.0, maximum(velocity))
        x := t
        y := velocity
        ()
    end
end


struct WaveWidthPlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
WaveWidthPlot(; output_name="wave_width.png", kwargs...) = WaveWidthPlot(output_name, kwargs)
@recipe function f(plot_spec::WaveWidthPlot, t::AbstractArray{T,1}, wave::AbstractArray{T,2}, transform=identity) where {T}
	title := "Wave half-width over time"
	xlab := "Time (a.u. approx. ms)"
	ylab := "Width (a.u. approx. um) "
    @series begin
        seriestype := :scatter
        x := t
        y := transform.(calculate_width(wave)')
        ()
    end
end

struct WaveStatsPlot{T} <: AbstractSpaceTimePlotSpecification
    output_name::String
    dt::T
    dx::T
    #smoothing::Union{T,Nothing}
    kwargs::Dict
end
WaveStatsPlot(; output_name="wave_stats.png", dt::T=nothing, dx::T=nothing, kwargs...) where T = WaveStatsPlot{T}(output_name, dt, dx, kwargs)
@recipe function f(plot_spec::WaveStatsPlot{T}, t::AbstractArray{T,1}, x::AbstractArray{T,1}, wave::AbstractArray{T,2})  where {T}
    layout := (2,2)
    legend := false
    delete!.(Ref(plotattributes), [:velocity_smoothing, :peak_interpolation_n])
    velocity_smoothing = pop!(plot_spec.kwargs, :velocity_smoothing, nothing)
    peak_interpolation_n = pop!(plot_spec.kwargs, :peak_interpolation_n, 1)

    begin
    @series begin
        subplot := 1
        (WaveWidthPlot(), t, wave)
    end
    @series begin
        subplot := 3
        (PeakTravelingWavePlot(dt=plot_spec.dt, dx=plot_spec.dx), t, x, wave, peak_interpolation_n)
    end
    # @series begin
    #     subplot := 2
    #     (WaveVelocityPlot(dt=plot_spec.dt), t, wave)
    # end
    @series begin
        subplot := 4
        # title := "Wave log velocity"
        # (WaveVelocityPlot(dt=plot_spec.dt, dx=plot_spec.dx), t, wave, (x) -> sign(x) * log10(abs(x)))
        (WaveVelocityPlot(;dt=plot_spec.dt, dx=plot_spec.dx,
				smoothing=velocity_smoothing,
				interpolation_n=peak_interpolation_n),
			t, x, wave)
    end
    end
    # plot(WaveWidthPlot(), wave, t, subplot=1)
    # plot(PeakTravelingWavePlot(dt=dt), simulation, subplot=2)
    # plot(WaveVelocityPlot(dt=dt), wave, t, subplot=3)
end


struct NeumanTravelingWavePlot{T} <: AbstractPlotSpecification
    output_name::String
    dt::T
    kwargs::Dict
end
NeumanTravelingWavePlot(; output_name="traveling_wave.png", dt::Union{Nothing,T}=nothing, kwargs...) where {T<:Float64} = NeumanTravelingWavePlot{T}(output_name, dt, kwargs)
@recipe function f(plot_spec::NeumanTravelingWavePlot{T}, execution::Execution{T,<:Simulation{T,M}}) where {T,M<:WCMSpatial}
    solution = execution.solution
    simulation = execution.simulation
    t = timepoints(simulation)
    space_origin::Int = origin_idx(simulation) # TODO: Remove 1D return assumption
    di = max(1, round(Int, simulation.solver.simulated_dt / plot_spec.dt))
    x = coordinates(simulation)[space_origin:di:end] # TODO: Remove 1D return assumption
	for time_dx in 1:length(t)
        @series begin
            seriestype := :line
            x := x
            y := solution[space_origin:di:end,:,time_dx] * [1.0, -1.0] # Subtract inhibitory activity...
            ()
        end
    end
end
