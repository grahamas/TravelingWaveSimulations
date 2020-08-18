using InteractiveUtils
using DrWatson
quickactivate(@__DIR__, "TravelingWaveSimulations")
using TravelingWaveSimulations

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
    "--no-save-raw"
        help = "Don't save raw simulation"
        action = :store_true
    "--backup-paths"
        nargs = '*'
        help = "Locations to copy the data to after completion (scp)"
end

@show (versioninfo())
@warn "$ARGS"
args = parse_args(ARGS, arg_settings; as_symbols=true)
@warn "$args"
based_on_example_serial(; args...)

#main()
