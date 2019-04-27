using DrWatson
quickactivate(@__DIR__)

using WilsonCowanModel

data_root = ARGS[1]
example_name = Symbol(ARGS[2])

@eval simulation = Examples.$(example_name)()
@eval plots = Examples.$(Symbol(example_name,:_plots))()

sim_output_path = joinpath(data_root, "sim", example_name)
plots_path = joinpath(plotsdir(), example_name)

execution = execute(simulation)
execution_dict = @dict execution
@tagsave(joinpath(sim_output_path, "$(Dates.now())_$(current_commit()).bson"), execution_dict, true)

plot_and_save(plots, Ref(execution), plots_path)
