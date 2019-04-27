using DrWatson
quickactivate(@__DIR__)

using WilsonCowanModel
using Simulation73
using Dates

ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
using Plots
pyplot()

data_root = ARGS[1]
example_name = ARGS[2]

@eval simulation = Examples.$(Symbol(example_name))()
@eval plots = Examples.$(Symbol(example_name,:_plots))()

sim_output_path = joinpath(data_root, "sim", example_name)
plots_path = joinpath(plotsdir(), example_name)
mkpath.([sim_output_path, plots_path])

execution = execute(simulation)
execution_dict = @dict execution
@tagsave(joinpath(sim_output_path, "$(Dates.now())_$(current_commit()).bson"), execution_dict, true)

plot_and_save.(plots, Ref(execution), plots_path)
