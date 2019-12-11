if @isdefined nworkers
    @everywhere using DrWatson
    @everywhere quickactivate(@__DIR__, "TravelingWaveSimulations")
    @everywhere using TravelingWaveSimulations
    @everywhere ENV["GKSwstype"] = "100" # For headless plotting (on server)
    @everywhere ENV["MPLBACKEND"]="Agg"
else
    using DrWatson
    quickactivate(@__DIR__, "TravelingWaveSimulations")
    using TravelingWaveSimulations
    ENV["GKSwstype"] = "100" # For headless plotting (on server)
    ENV["MPLBACKEND"]="Agg"
end

using Plots
#pyplot()

using ArgParse
arg_settings = ArgParseSettings(; autofix_names = true)
@add_arg_table arg_settings begin
    "--data-root", "-d"
        help = "Location in which to store output data"
        default = datadir()
    "--example-name"
        help = "Name of example defined examples.jl"
    "--modifications", "--mod"
        nargs = '*'
        help = "Name of file specifying dict of modifications"
    "--analyses"
        nargs = '*'
        help = "Name of file specifying analyses"
    "--no-save-raw"
        help = "Don't save raw simulation"
        action = :store_true
    "--batch"
        default = 1000
        arg_type = Int
    "--max-sims-in-mem"
        arg_type = Int
end

@warn "$ARGS"
args = parse_args(ARGS, arg_settings; as_symbols=true)
@warn "$args"
based_on_example(; args...)

#main()
