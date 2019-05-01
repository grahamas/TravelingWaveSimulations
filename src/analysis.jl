
export SubsampledPlot

export NonlinearityPlot
export Animate
export SpaceTimePlot

export NeumanTravelingWavePlot #, TravelingWavePlot
export PeakTravelingWavePlot #, TravelingWavePlot

export WaveVelocityPlot

export WaveWidthPlot

export WaveStatsPlot

export plot_and_save

function full_name(name; path="", prefix="", sep="_")
	if prefix != ""
		name = join([prefix, name], sep)
	end
	return joinpath(path, name)
end

function plot_and_save(plot_spec::AbstractPlotSpecification, execution::Execution, output_dir::AbstractString, prefix="")
	plot_obj = RecipesBase.plot(plot_spec, execution; plot_spec.kwargs...)
	path = full_name(output_name(plot_spec); path=output_dir, prefix=prefix) # Must be after plotting in case name changes
	recursively_clear_path(path)
	savefig(plot_obj, path)
end

#region Animate
struct Animate <: AbstractPlotSpecification
    fps::Int
    output_name::String
    kwargs::Dict
end
Animate(; fps=20, output_name="animation.mp4", kwargs...) = Animate(fps, output_name, kwargs)
function plot_and_save(plot_spec::Animate, execution::Execution, output_dir::AbstractString, prefix="")
    path = full_name(output_name(plot_spec); path=output_dir, prefix=prefix)
    recursively_clear_path(path)
    anim = RecipesBase.animate(execution; plot_spec.kwargs...)
    mp4(anim, path; fps=plot_spec.fps)
end
function RecipesBase.animate(execution::Execution{T,<:Simulation{T,M}}; kwargs...) where {T,M<:WCMSpatial{T,1}}
    solution = execution.solution
    simulation = execution.simulation
    pop_names = simulation.model.pop_names
    x = saved_space_arr(simulation)
    t = saved_time_arr(simulation)
    max_val = maximum(solution)
    @animate for time_dx in 1:length(t) # TODO @views
        plot(x, pop_frame(solution, 1, time_dx); label=pop_names[1],
            ylim=(0,max_val), title="t = $(round(t[time_dx], digits=4))", kwargs...)
        for i_pop in 2:length(pop_names)
            plot!(x, pop_frame(solution, i_pop, time_dx); label=pop_names[i_pop], kwargs...)
        end
    end
end
function RecipesBase.animate(execution::Execution{T,<:Simulation{T,M}}; kwargs...) where {T,M<:WCMSpatial{T,2}}
    solution = execution.solution
    simulation = execution.simulation
    pop_names = simulation.model.pop_names
    x = saved_space_arr(simulation)
    xs = [coords[1] for coords in x[1,:]]
    ys = [coords[2] for coords in x[:,1]]
    t = saved_time_arr(simulation)
    max_val = maximum(solution)
    @animate for time_dx in 1:length(t) # TODO @views
        plot(
            [heatmap( pop_frame(solution, i_pop, time_dx);
                zlim=(0,max_val), clim=(0,max_val),
                title="t = $(round(t[time_dx], digits=4)), $(pop_names[i_pop])", kwargs...)
                for i_pop in 1:length(pop_names)]...
                    )
    end
end
#endregion

#region SpaceTimePlot
struct SpaceTimePlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
SpaceTimePlot(; output_name = "spacetime.png", kwargs...) = SpaceTimePlot(output_name, kwargs)
@recipe function f(plot_spec::SpaceTimePlot, execution::Execution{T,<:Simulation{T,M}}) where {T,M<:WCMSpatial}
    simulation = execution.simulation
    solution = execution.solution
    v_space = saved_space_arr(simulation)
    v_time = saved_time_arr(simulation)
    clims := (minimum(solution), maximum(solution))
    grid := false
    layout := (2,1)
    for i_pop in 1:2 # TODO!!
        @series begin
            seriestype --> :heatmap
            subplot := i_pop
            x := v_time
            y := v_space
            solution[:,i_pop,:]
        end
    end
end
#endregion

#region NonlinearityPlot
struct NonlinearityPlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
NonlinearityPlot(; output_name = "nonlinearity.png", kwargs...) = NonlinearityPlot(output_name, kwargs)
@recipe function f(plot_spec::NonlinearityPlot, execution::Execution{T,<:Simulation{T,M}}; resolution=100, fn_bounds=(-1.0,15.0)) where {T,M<:WCMSpatial}
    simulation = execution.simulation

    pop_names = simulation.model.pop_names
    n_pops = length(pop_names)

    nonlinearity_mutator! = make_mutator(simulation.model.nonlinearity)

    xlab := "Input current"
    ylab := "Proportion pop. reaching at least threshold"

    one_pop_x = range(fn_bounds[1], stop=fn_bounds[2], length=resolution)
    output = repeat(one_pop_x, outer=(1,n_pops))
    nonlinearity_mutator!(output, output, 0.0)

    for i_pop in 1:length(pop_names)
        @series begin
            lab --> pop_names[i_pop]
            seriestype := :line
            x := one_pop_x
            y := output[:,i_pop]
            ()
        end
    end
end

struct ConnectivityPlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
ConnectivityPlot(; output_name = "connectivity.png", kwargs...) = ConnectivityPlot(output_name, kwargs)
# @recipe function f(plot_spec::ConnectivityPlot, execution::E) where {T,P,M<:WCMSpatial{T,1,P}, S<:Simulation{T,M}, E<:Execution{T,S}}
#     simulation = execution.simulation
#
#     pop_names = simulation.model.pop_names
#
#     connectivity = calculate.(simulation.model.connectivity, Ref(Calculated(simulation.model.space)))
#     map(Iterators.product(1:P, 1:P)) do (dst_pop, src_pop)
#         @series begin
#             lab --> "$(pop_names[src_pop]) → $(pop_names[dst_pop])"
#             seriestype := :heatmap
#             subplot := (dst_pop, src_pop)
#             connectivity[:,:,dst_pop,src_pop]
#         end
#     end
# end
@recipe function f(plot_spec::ConnectivityPlot, execution::E) where {T,P,M<:WCMSpatial{T,2,P}, S<:Simulation{T,M}, E<:Execution{T,S}}
    simulation = execution.simulation

    pop_names = simulation.model.pop_names

    connectivity = tensor(simulation.model.connectivity, simulation.model.space)
    layout := (2,2)
    for (dst_pop, src_pop) in Iterators.product(1:P, 1:P)
        @series begin
            lab --> "$(pop_names[src_pop])[1,1] → $(pop_names[dst_pop])"
            seriestype := :heatmap
            subplot := dst_pop + (src_pop-1) * P
            connectivity[:,:,1,1,dst_pop,src_pop]
        end
    end
end

#endregion

#region NeumanTravelingWavePlot
struct NeumanTravelingWavePlot{T} <: AbstractPlotSpecification
    output_name::String
    dt::T
    kwargs::Dict
end
NeumanTravelingWavePlot(; output_name="traveling_wave.png", dt::Union{Nothing,T}=nothing, kwargs...) where {T<:Float64} = NeumanTravelingWavePlot{T}(output_name, dt, kwargs)
@recipe function f(plot_spec::NeumanTravelingWavePlot{T}, execution::Execution{T,<:Simulation{T,M}}) where {T,M<:WCMSpatial}
    solution = execution.solution
    simulation = execution.simulation
    t = saved_time_arr(simulation)
    space_origin::Int = get_origin(simulation) # TODO: Remove 1D return assumption
    di = max(1, round(Int, simulation.solver.simulated_dt / plot_spec.dt))
    x = saved_space_arr(simulation)[space_origin:di:end] # TODO: Remove 1D return assumption
    for time_dx in 1:length(t)
        @series begin
            seriestype := :line
            x := x
            y := solution[space_origin:di:end,:,time_dx] * [1.0, -1.0] # Subtract inhibitory activity...
            ()
        end
    end
end
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




function wave_maxima(single_wave_data::SPACE1DTIME) where {T, SPACE1DTIME<:AbstractArray{T,2}}
	# [space, time]
	max_vals, max_ixs = findmax(single_wave_data, dims=1)
	return (max_vals, max_ixs)
end

function interpolate_true_max_ix(max_ix, arr::SPACE1D) where {T, SPACE1D<:AbstractArray{T,1}}
	circa_max_ixs = (max_ix - CartesianIndex(1,0))
end

function interpolate_parabola(space, wave)
	@. parabola(x,p) = p[1] + p[2] * ((x - p[3]) ^ 2)
	ub = [Inf, 0.0, space[end]]
	lb = [minimum(wave), -Inf, space[1]]
	fit = curve_fit(parabola, [maximum(wave), 0.0, 0.0])
	return (coef(fit)[3], coef(fit)[1])

function track_wave_peak(x::SPACE1D, wave::SPACE1DTIME) where {T, SPACE1D<:AbstractArray{T,1}, SPACE1DTIME<:AbstractArray{T,2}}
	max_vals, max_ixs = wave_maxima(wave)
	space_ixs = [ix[1] for ix in max_ixs]
	side = 1
	circa_space_ixs = [(ix-side):(ix+side) for ix in space_ixs]
	circa_wave_ixs = [(ix-CartesianIndex(side,0)):(ix+CartesianIndex(side,0)) for ix in max_ixs]
	interpolated_xs, interpolated_vals = map(zip(circa_space_ixs, circa_wave_ixs)) do (circa_space_ix, circa_wave_ix)
		interpolate_parabola(x[circa_space_ix], wave[circa_wave_ix])
	end
	return interpolated_xs, interpolated_vals
end

function calculate_wave_velocity(x::SPACE1D, wave::SPACE1DTIME, dt::T=one(T)) where {T, SPACE1D<:AbstractArray{T,1}, SPACE1DTIME<:AbstractArray{T,2}}
    # [space, time]
	xs, vals = track_wave_peak(x, wave)
	diff(xs) ./ dt
end


#region SubsampledPlot
mutable struct SubsampledPlot <: AbstractPlotSpecification
    plot_type::Type{<:AbstractSpaceTimePlotSpecification}
    time_subsampler::Subsampler
    space_subsampler::Subsampler
    output_name::String
    kwargs::Iterators.Pairs
end
SubsampledPlot(; plot_type=nothing, time_subsampler=Subsampler(), space_subsampler=Subsampler(), output_name="", kwargs...) = SubsampledPlot(plot_type, time_subsampler, space_subsampler, output_name, kwargs)
@recipe function f(subsampledplot::SubsampledPlot, execution::Execution{T,<:Simulation{T,M}}) where {T,M<:WCMSpatial}
    simulation = execution.simulation
    t, x, wave = subsample(execution, time_subsampler=subsampledplot.time_subsampler, space_subsampler=subsampledplot.space_subsampler)

    dt = subsampledplot.time_subsampler.Δ == nothing ? save_dt(simulation) : subsampledplot.time_subsampler.Δ

    dx = subsampledplot.space_subsampler.Δ == nothing ? save_dx(simulation) : subsampledplot.space_subsampler.Δ

    plot_spec = subsampledplot.plot_type(;dt=dt, dx=dx, subsampledplot.kwargs...)
    if subsampledplot.output_name == ""
        subsampledplot.output_name = plot_spec.output_name
    end

    (plot_spec, t, x, wave)
end
#endregion

#region PeakTravelingWavePlot
struct PeakTravelingWavePlot <: AbstractPlotSpecification
    output_name::String
    kwargs::Dict
end
PeakTravelingWavePlot(; output_name="peak_traveling_wave.png", kwargs...) = PeakTravelingWavePlot(output_name, kwargs)
@recipe function f(plot_spec::PeakTravelingWavePlot, t::AbstractArray{T,1}, x::AbstractArray{T,1}, wave::AbstractArray{T,2}) where {T}
    title := "Traveling wave with labeled peaks"
    for time_dx in 1:length(t)
        @series begin
            seriestype := :line
            x := x
            y := wave[:,time_dx] # Subtract inhibitory activity...
            ()
        end
    end
    xs, peaks = track_wave_peak(x, wave)
    @series begin
        seriestype := :scatter
        x := xs
        y := peaks
        ()
    end
end
#endregion

#region WaveStats

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
    kwargs::Dict
end
WaveVelocityPlot(; output_name="wave_velocity.png", dt::Union{Nothing,T}=nothing, dx::Union{Nothing,T}=nothing, smoothing::Union{Nothing,T}=nothing, kwargs...) where {T<:Float64} = WaveVelocityPlot{T}(output_name, dt, dx, smoothing, kwargs)
@recipe function f(plot_spec::WaveVelocityPlot{T}, t::AbstractArray{T,1}, x::AbstractArray{T,1}, wave::AbstractArray{T,2}, transform::Function=identity, naive=false) where {T}
    velocity = calculate_wave_velocity(x, wave, plot_spec.dt, plot_spec.dx)
    if plot_spec.smoothing != nothing
        t, velocity = smooth(t, velocity, plot_spec.smoothing/plot_spec.dt)
    end
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
    @series begin
        title --> "Wave width over time"
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
    smoothing::Union{T,Nothing}
    kwargs::Dict
end
WaveStatsPlot(; output_name="wave_stats.png", dt::T=nothing, dx::T=nothing, smoothing::Union{T,Nothing}=nothing, kwargs...) where T = WaveStatsPlot{T}(output_name, dt, dx, smoothing, kwargs)
@recipe function f(plot_spec::WaveStatsPlot{T}, t::AbstractArray{T,1}, x::AbstractArray{T,1}, wave::AbstractArray{T,2})  where {T}
    layout := (2,2)
    legend := false
    @warn "got here"
    delete!(plotattributes, :smoothing)

    begin
    @series begin
        subplot := 1
        (WaveWidthPlot(), t, wave)
    end
    @series begin
        subplot := 3
        (PeakTravelingWavePlot(dt=plot_spec.dt, dx=plot_spec.dx), t, x, wave)
    end
    # @series begin
    #     subplot := 2
    #     (WaveVelocityPlot(dt=plot_spec.dt), t, wave)
    # end
    @series begin
        subplot := 4
        # title := "Wave log velocity"
        # (WaveVelocityPlot(dt=plot_spec.dt, dx=plot_spec.dx), t, wave, (x) -> sign(x) * log10(abs(x)))
        (WaveVelocityPlot(;dt=plot_spec.dt, dx=plot_spec.dx, smoothing=plot_spec.smoothing), t, x, wave)
    end
    end
    # plot(WaveWidthPlot(), wave, t, subplot=1)
    # plot(PeakTravelingWavePlot(dt=dt), simulation, subplot=2)
    # plot(WaveVelocityPlot(dt=dt), wave, t, subplot=3)
end
#endregion
