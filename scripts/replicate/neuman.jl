# ENV["GKSwstype"] = "100" # For headless plotting (on server)
# ENV["MPLBACKEND"]="Agg"
using DrWatson
quickactivate(@__DIR__, "WilsonCowanModel")
using BSON
using Simulation73
using WilsonCowanModel
using DifferentialEquations: Euler
using Random
using Dates

simulation = WilsonCowanModel.Examples.replicate_neuman()

execution = execute(simulation);

current_time = string(Dates.now())
generic_replication_dir = joinpath("replicate", "neuman")
replication_directory = joinpath(datadir(), "sim", generic_replication_dir, current_time)
this_commit_filename = current_commit()*".bson"
mkpath(replication_directory)
tagsave(joinpath(replication_directory, this_commit_filename), @dict execution; safe=true)




plots = [
Animate(;
  fps = 20
  ),
NonlinearityPlot(;
  fn_bounds = (-1,15)
  ),
# SpaceTimePlot(),
SubsampledPlot(
  plot_type=WaveStatsPlot,
  time_subsampler=Subsampler(
    Δ = 0.01,
    window = (1.2, 1.8)
  ),
  space_subsampler=Subsampler(
      window = (5.0,Inf)
    )
  )
]

plot_path = joinpath(plotsdir(), generic_replication_dir)
mkpath(plot_path)
plot_and_save.(plots, Ref(execution), plot_path)
