using DrWatson
quickactivate(@__DIR__, "TravelingWaveSimulations")

using Simulation73, TravelingWaveSimulations, NeuralModels
using BSON
using Dates
using ArgParse
using Lazy
ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
using Plots
pyplot()

arg_settings = ArgParseSettings()
@add_arg_table arg_settings begin
    "--data-root", "-d"
        help = "Location in which to store output data"
        default = datadir()
    "--example-name"
        help = "Name of example defined examples.jl"
    "--modifications", "--mod"
        nargs = '*'
        help = "Name of file specifying dict of modifications"
    "--plot-specs", "--plot"
        nargs = '*'
        help = "Name of file specifying plots"
    "--no-save-raw"
        help = "Don't save raw simulation"
        action = :store_true
end

args = parse_args(ARGS, arg_settings)
@show args
based_on_example(; args...)

#main()
