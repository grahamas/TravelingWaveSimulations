module WCMAnalysis
# Separated into own module to prevent reloading full module during analysis stage
# of research.

using Simulation73, WilsonCowanModel,
	Plots

export SubsampledPlot, NonlinearityPlot, SpaceTimePlot,
	NeumanTravelingWavePlot,PeakTravelingWavePlot,
	WaveStatsPlot, WaveWidthPlot, WaveVelocityPlot,
	Animate,
	plot_and_save

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


include("animation.jl")
include("analysis/animate.jl")
include("analysis/spacetime.jl")
include("analysis/nonlinearity.jl")
include("analysis/connectivity.jl")
include("analysis/traveling_wave_statistics.jl")
include("analysis/peak_tracking.jl")
include("analysis/subsampled_plot.jl")

end
