module WCMAnalysis

#region imports
using Modeling
using Analysis
using Records
import DifferentialEquations: DESolution
using Meshes
using RecipesBase
using Plots
using Parameters
using WCM
using CalculatedParameters
using Simulating
using WCMNonlinearity
using Subsampling
#endregion

#region Animate
struct Animate <: AbstractPlotSpecification
    fps::Int
    output_name::String
    kwargs::Dict
end
Animate(; fps=20, output_name="animation.mp4", kwargs...) = Animate(fps, output_name, kwargs)
function Analysis.plot_and_save(plot_spec::Animate, simulation::Simulation)
    save_fn(name, anim) = mp4(anim, name; fps=plot_spec.fps)
    simulation.output(save_fn, output_name(plot_spec), animate(simulation; plot_spec.kwargs...))
end
export Animate
function RecipesBase.animate(simulation::Simulation{T,M}; kwargs...) where {T,M<:WCMSpatial1D}
    solution = simulation.solution
    pop_names = simulation.model.pop_names
    x = space_arr(simulation)
    t = time_arr(simulation)
    max_val = maximum(simulation)
    @animate for time_dx in 1:length(t) # TODO @views
        plot(x, solution[:, 1, time_dx]; label=pop_names[1],
            ylim=(0,max_val), title="t = $(round(t[time_dx], digits=4))", kwargs...)
        for i_pop in 2:length(pop_names)
            plot!(x, solution[:, i_pop, time_dx]; label=pop_names[i_pop], kwargs...)
        end

    end
end
#endregion

#region SpaceTimePlot
struct SpaceTimePlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
SpaceTimePlot(; output_name = "spacetime.png", kwargs...) = SpaceTimePlot(output_name, kwargs)
@recipe function f(plot_spec::SpaceTimePlot, simulation::Simulation{T,M}) where {T,M<:WCMSpatial1D}
    v_space = space_arr(simulation)
    v_time = time_arr(simulation)
    clims := (minimum(simulation), maximum(simulation))
    grid := false
    layout := (2,1)
    for i_pop in 1:2 # TODO!!
        @series begin
            seriestype --> :heatmap
            subplot := i_pop
            x := v_time
            y := v_space
            simulation.solution[:,i_pop,:]
        end
    end
end
export SpaceTimePlot
#endregion

#region NonlinearityPlot
struct NonlinearityPlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
NonlinearityPlot(; output_name = "nonlinearity.png", kwargs...) = NonlinearityPlot(output_name, kwargs)
@recipe function f(plot_spec::NonlinearityPlot, simulation::Simulation{T,M}; resolution=100, fn_bounds=(-1.0,15.0)) where {T,M<:WCMSpatial1D}
    pop_names = simulation.model.pop_names
    n_pops = length(pop_names)

    nonlinearity_fns = get_value.(Calculated(simulation.model).nonlinearity)

    one_pop_x = 
    #delete!.(Ref(plotattributes),[:resolution,:fn_bounds])

    xlab := "Input current"
    ylab := "Proportion pop. reaching at least threshold"

    one_pop_x = range(fn_bounds[1], stop=fn_bounds[2], length=resolution)

    for i_pop in 1:length(pop_names)
        @series begin
            lab --> pop_names[i_pop]
            seriestype := :line
            x := one_pop_x
            y := nonlinearity(nonlinearity_fns[i_pop], collect(one_pop_x))
            ()
        end
    end
end
export NonlinearityPlot
#endregion

#region NeumanTravelingWavePlot
struct NeumanTravelingWavePlot{T} <: AbstractPlotSpecification
    output_name::String
    dt::T
    kwargs::Dict
end
NeumanTravelingWavePlot(; output_name="traveling_wave.png", dt::Union{Nothing,T}=nothing, kwargs...) where {T<:Float64} = NeumanTravelingWavePlot{T}(output_name, dt, kwargs)
@recipe function f(plot_spec::NeumanTravelingWavePlot{T}, simulation::Simulation{T,M}) where {T,M<:WCMSpatial1D}
    t = time_arr(simulation)
    space_origin::Int = get_origin(simulation) # TODO: Remove 1D return assumption
    di = max(1, round(Int, simulation.solver.simulated_dt / plot_spec.dt))
    x = space_arr(simulation)[space_origin:di:end] # TODO: Remove 1D return assumption
    for time_dx in 1:length(t)
        @series begin
            seriestype := :line
            x := x
            y := simulation.solution[space_origin:di:end,:,time_dx] * [1.0, -1.0] # Subtract inhibitory activity...
            ()
        end
    end
end
export NeumanTravelingWavePlot #, TravelingWavePlot
#endregion

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

"""
    During sparse movement (velocity as difference of location usually 0),
    calculate velocity as the amount of movement divided by the length of
    the window in which no movement but that was recorded. The window is 
    determined by bisecting each run of no-movement.

    Discard first and last movements.
"""
# function calculate_velocity(single_wave_data::AT, dt::T=one(T)) where {T, AT<:AbstractArray{T,2}}
#     # [space, time]
#     naive_velocity = calculate_naive_velocity(single_wave_data, 1.0)'
#     @show naive_velocity
#     nonzero_vel = naive_velocity .!= zero(T)
#     n_moves = sum(nonzero_vel)
#     if n_moves < 3
#         return (T[], T[])
#     end
    
#     first_movement = findfirst(nonzero_vel)[1]
#     second_movement = findfirst(nonzero_vel[first_movement+1:end])[1]
#     prev_leading_window = dt * (second_movement - first_movement) / 2.0
#     prev_movement = second_movement

#     t = Array{T,1}(undef, n_moves-2) # Doesn't store the first or last move
#     velocity = similar(t)
#     i_velocity = 1
#     for i_naive in (second_movement+1):length(naive_velocity)
#         vel = naive_velocity[i_naive]
#         if vel != zero(T)
#             this_movement = i_naive
#             leading_window = dt * (this_movement - prev_movement) / 2.0
#             prev_total_window = prev_leading_window + leading_window
#             @show i_velocity, i_naive, n_moves, length(naive_velocity)
#             velocity[i_velocity] = naive_velocity[i_naive] / prev_total_window
#             t[i_velocity] = prev_movement * dt

#             i_velocity += 1
#             prev_movement = this_movement
#         end
#     end
#     return (t, velocity)
# end

function calculate_naive_velocity(single_wave_data::AT, dt::T=one(T)) where {T, AT<:AbstractArray{T,2}}
    # [space, time]
    peak_through_time = track_wave_peak(single_wave_data)[2]
    peak_through_time = [cart[1] for cart in peak_through_time]
    peak_vel = diff(peak_through_time, dims=2) ./ dt
end

function track_wave_peak(single_wave_data::AT) where {T, AT<:AbstractArray{T,2}}
    # [space, time]
    #space_diff = diff(single_wave_data, dims=1)
    space_max = findmax(single_wave_data, dims=1) # TODO: Make this more robust
end


mutable struct SubsampledPlot <: AbstractPlotSpecification
    plot_type::Type{<:AbstractSpaceTimePlotSpecification}
    time_subsampling::Dict
    space_subsampling::Dict
    output_name::String
    kwargs::Iterators.Pairs
end
SubsampledPlot(; plot_type=nothing, time_subsampling=Dict(), space_subsampling=Dict(), output_name="", kwargs...) = SubsampledPlot(plot_type, time_subsampling, space_subsampling, output_name, kwargs)
@recipe function f(subsampledplot::SubsampledPlot, simulation::Simulation{T,M}) where {T,M<:WCMSpatial1D}
    space_origin::Int = get_origin(simulation)
    t = time_arr(simulation)
    x = space_arr(simulation)

    t_dxs = subsampling_idxs(save_dt(simulation), length(t); origin_idx=1, subsampledplot.time_subsampling...)
    x_dxs = subsampling_idxs(save_dx(simulation), length(x); origin_idx=space_origin, subsampledplot.space_subsampling...)
    pop_dxs = 1

    dt = get(subsampledplot.time_subsampling, :Δsubsampled) do 
        save_dt(simulation)
    end
    dx = get(subsampledplot.space_subsampling, :Δsubsampled) do 
        save_dx(simulation)
    end

    t = t[t_dxs]
    x = x[x_dxs] # TODO: Remove 1D return assumption
    wave = simulation.solution[x_dxs,pop_dxs,t_dxs]

    plot_spec = subsampledplot.plot_type(dt=dt, dx=dx, subsampledplot.kwargs...)
    if subsampledplot.output_name == ""
        subsampledplot.output_name = plot_spec.output_name
    end

    (plot_spec, t, x, wave)
end
export SubsampledPlot

#region PeakTravelingWavePlot
struct PeakTravelingWavePlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
PeakTravelingWavePlot(; output_name="peak_traveling_wave.png", kwargs...) = PeakTravelingWavePlot(output_name, kwargs)
@recipe function f(plot_spec::PeakTravelingWavePlot, t::AT1, x::AT1, wave::AT2) where {T,AT1<: AbstractArray{T,1},AT2<:AbstractArray{T,2}}
    title := "Traveling wave with labeled peaks"
    for time_dx in 1:length(t)
        @series begin
            seriestype := :line
            x := x
            y := wave[:,time_dx] # Subtract inhibitory activity...
            ()
        end
    end
    peaks, peak_dxs = track_wave_peak(wave)
    x_peak_dxs = [dx[1] for dx in peak_dxs]
    @series begin
        seriestype := :scatter
        x := x[x_peak_dxs]'
        y := peaks'
        ()
    end
end
export PeakTravelingWavePlot #, TravelingWavePlot
#endregion

#region WaveStats
struct WaveVelocityPlot{T} <: AbstractPlotSpecification
    output_name::String
    dt::T
    kwargs::Dict
end
WaveVelocityPlot(; output_name="wave_velocity.png", dt::Union{Nothing,T}=nothing, kwargs...) where {T<:Float64} = WaveVelocityPlot{T}(output_name, dt, kwargs)
@recipe function f(plot_spec::WaveVelocityPlot{T}, t::AT1, wave::AT2, transform::Function=identity, naive=false) where {T,AT1<: AbstractArray{T,1},AT2<:AbstractArray{T,2}}
    @series begin
        title --> "Wave velocity over time"
        seriestype := :scatter
        t, velocity = calculate_velocity(wave, plot_spec.dt)
        x := t
        y := velocity
        ()
    end
end
export WaveVelocityPlot

struct WaveNaiveVelocityPlot{T} <: AbstractPlotSpecification
    output_name::String
    dt::T
    kwargs::Dict
end
WaveNaiveVelocityPlot(; output_name="wave_velocity.png", dt::Union{Nothing,T}=nothing, kwargs...) where {T<:Float64} = WaveNaiveVelocityPlot{T}(output_name, dt, kwargs)
@recipe function f(plot_spec::WaveNaiveVelocityPlot{T}, t::AT1, wave::AT2, transform::Function=identity) where {T,AT1<: AbstractArray{T,1},AT2<:AbstractArray{T,2}}
    @series begin
        title --> "Wave velocity over time"
        seriestype := :scatter
        velocity = calculate_naive_velocity(wave, plot_spec.dt)
        x := t
        y := velocity
        ()
    end
end
export WaveNaiveVelocityPlot

struct WaveWidthPlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
WaveWidthPlot(; output_name="wave_width.png", kwargs...) = WaveWidthPlot(output_name, kwargs)
@recipe function f(plot_spec::WaveWidthPlot, t::AT1, wave::AT2, transform=identity) where {T,AT1<: AbstractArray{T,1},AT2<:AbstractArray{T,2}}
    @series begin
        title --> "Wave width over time"
        seriestype := :scatter
        x := t
        y := transform.(calculate_width(wave)')
        ()
    end
end
export WaveWidthPlot

struct WaveStatsPlot{T} <: AbstractSpaceTimePlotSpecification
    output_name::String
    dt::T
    kwargs::Dict
end
WaveStatsPlot(; output_name="wave_stats.png", dt::Union{Nothing,T}=nothing, kwargs...) where {T<:Float64} = WaveStatsPlot{T}(output_name, dt, kwargs)
@recipe function f(plot_spec::WaveStatsPlot{T}, t::AT1, x::AT1, wave::AT2)  where {T,AT1<: AbstractArray{T,1},AT2<:AbstractArray{T,2}}
    layout := (2,2)
    legend := false

    @series begin
        subplot := 1
        (WaveWidthPlot(), t, wave)
    end
    @series begin
        subplot := 3
        (PeakTravelingWavePlot(dt=plot_spec.dt), t, x, wave)
    end
    # @series begin
    #     subplot := 2
    #     (WaveVelocityPlot(dt=plot_spec.dt), t, wave)
    # end
    @series begin
        subplot := 4
        # title := "Wave log velocity"
        # (WaveVelocityPlot(dt=plot_spec.dt), t, wave, (x) -> sign(x) * log10(abs(x)))
        title := "Wave naive velocity"
        (WaveNaiveVelocityPlot(dt=plot_spec.dt), t, wave)
    end
    # plot(WaveWidthPlot(), wave, t, subplot=1)
    # plot(PeakTravelingWavePlot(dt=dt), simulation, subplot=2)
    # plot(WaveVelocityPlot(dt=dt), wave, t, subplot=3)
end
export WaveStatsPlot
#endregion

end