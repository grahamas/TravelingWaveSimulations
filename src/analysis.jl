export plot_and_save

function full_name(name; path="", prefix="", sep="_")
	if prefix != ""
		name = join([prefix, name], sep)
	end
	full_path = joinpath(path, name)
	@assert length(full_path) < 4096
	return full_path
end

function plot_and_save(plot_spec::AbstractPlotSpecification, execution::Execution, output_dir::AbstractString, prefix="")
	plot_obj = RecipesBase.plot(plot_spec, execution; plot_spec.kwargs...)
	path = full_name(output_name(plot_spec); path=output_dir, prefix=prefix) # Must be after plotting in case name changes
	recursively_clear_path(path)
	savefig(plot_obj, path)
end

analysis_path = "analysis"
include(joinpath(analysis_path,"animate.jl"))
include(joinpath(analysis_path,"spacetime.jl"))
include(joinpath(analysis_path,"nonlinearity.jl"))
include(joinpath(analysis_path,"connectivity.jl"))
include(joinpath(analysis_path,"traveling_wave_statistics.jl"))
include(joinpath(analysis_path,"peak_tracking.jl"))
include(joinpath(analysis_path,"subsampled_plot.jl"))
