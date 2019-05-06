using DrWatson
quickactivate(@__DIR__)

using WilsonCowanModel
using Simulation73
using Dates
using ArgParse
using Lazy

arg_settings = ArgParseSettings()
@add_arg_table arg_settings begin
    "--data-root", "-d"
        help = "Location in which to store output data"
        default = datadir()
    "--example-name"
        help = "Name of example defined in examples.jl"
    "--modifications-case"
        help = "Name of file specifying dict of modifications"
    "--plotspec-case"
        help = "Name of file specifying plots"
    "--no-save-raw"
        help = "Don't save raw simulation"
        action = :store_true
end
args = parse_args(ARGS, arg_settings)
@show args
data_root = args["data-root"]
example_name = args["example-name"]
modifications_case = args["modifications-case"]
plotspec_case = args["plotspec-case"]

ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
using Plots
pyplot()

function read_modification_file(filename::String)
    include(joinpath(scriptdir(), "modifications", "$(filename).jl"))
    return modifications
end

function parse_modification(str::String)
    if ("=" in str)
        name_str, value_str = split(str, "=")
        return Dict(Symbol(name_str) => parse(Float64, value_str))
    else
        return read_modification_file(str)
    end
end

must_be_list(x::AbstractArray) = x
must_be_list(x) = [x]

function parse_modifications_array(array_str::String)
    @assert array_str[1] == "[" && array_str[end] == "]"
    parsed_modifications = @> array_str[2:end-1] begin
        split(",")
        strip.()
        parse_modification.()
        must_be_list.()
    modification_cases = Iterators.product(parsed_modifications...) .|> (case) -> merge(case...)
    return modification_cases
end

if modifications_case != nothing
    if modifications_case[1] == "["
        parse_modifications_array(modifications_case)
    else
        parse_modification(modifications_case) |> must_be_list
    end
    modifications_prefix = "$(modifications_case)_"
else
    modifications = [Dict()]
    modifications_prefix = ""
end

if plotspec_case != nothing
    plotspec_path = joinpath(scriptdir(), "plotspecs", "$(plotspec_case).jl")
    include(joinpath(scriptdir(), "plotspecs", plotspec_path)) # defines plotspecs
    plots_path = joinpath(plotsdir(), example_name)
    modification_prefix = ""
    plots_path = joinpath(plots_path, "$(modifications_prefix)$(plotspec_case)_$(Dates.now())_$(current_commit())")
    mkpath(plots_path)
else
    plotspecs = []
end

if !args["no-save-raw"]
    sim_output_path = joinpath(data_root, "sim", example_name, "$(modifications_prefix)$(Dates.now())_$(current_commit())")
    mkpath(sim_output_path)
end

for modification in modifications
    simulation = Examples.get_example(example_name)(; modification...)
    execution = execute(simulation)
    mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray))
    if !args["no-save-raw"]
        execution_dict = @dict execution
        @tagsave("$(joinpath(sim_output_path, mod_name)).bson", execution_dict, true)
    end
    plot_and_save.(plotspecs, Ref(execution), plots_path, mod_name)
end
