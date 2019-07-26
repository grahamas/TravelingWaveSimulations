using DrWatson
quickactivate(@__DIR__, "TravelingWaveSimulations")

using TravelingWaveSimulations
ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
using Plots
pyplot()

based_on_example(ARGS)

