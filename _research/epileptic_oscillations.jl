using DrWatson; quickactivate(@__DIR__, "TravelingWaveSimulations")
ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
using Plots; pyplot()

using TravelingWaveSimulations


