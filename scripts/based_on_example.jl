using DrWatson
quickactivate(@__DIR__)

using WilsonCowanModel, WCMExamples, WCMAnalysis
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
        help = "Name of example defined in WCMExamples"
    "--modifications-cases"
        nargs = '*'
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
modifications_cases = args["modifications-cases"]
plotspec_case = args["plotspec-case"]

ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
using Plots
pyplot()

function read_modification_file(filename::AbstractString)
    include(joinpath(scriptdir(), "modifications", "$(filename).jl"))
    return modifications
end

parse_range(start, stop) = parse(Float64, start):parse(Float64, stop)
parse_range(start, step, stop) = parse(Float64, start):parse(Float64, step):parse(Float64, stop)

function parse_modification(str::AbstractString)
    if occursin("=", str)
        name_str, value_str = split(str, "=")
        if occursin(":", value_str)
            return [Dict(Symbol(name_str) => val) for val in parse_range(split(value_str,":")...)]
        else
            return Dict(Symbol(name_str) => parse(Float64, value_str))
        end
    else
        return read_modification_file(str)
    end
end

must_be_list(x::AbstractArray) = x
must_be_list(x) = [x]

function parse_modifications_array(modification_strs::AbstractArray)
    parsed_modifications = @> modification_strs begin
        parse_modification.()
        must_be_list.()
    end
    modification_cases = Iterators.product(parsed_modifications...)
    modification_cases = map((case) -> merge(case...), modification_cases)
    return modification_cases
end

if modifications_cases != nothing
    modifications = parse_modifications_array(modifications_cases)
    modifications_prefix = """$(join(sort(modifications_cases), ";"))_"""
else
    modifications = [Dict()]
    modifications_prefix = ""
end

if plotspec_case != nothing
    plotspec_path = joinpath(scriptdir(), "plotspecs", "$(plotspec_case).jl")
    include(joinpath(scriptdir(), "plotspecs", plotspec_path)) # defines plotspecs
    plots_path = joinpath(plotsdir(), example_name)
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
    mod_name = savename(modification; allowedtypes=(Real,String,Symbol,AbstractArray), connector=";")
    if !args["no-save-raw"]
        execution_dict = @dict execution
        @tagsave("$(joinpath(sim_output_path, mod_name)).bson", execution_dict, true)
    end
    plot_and_save.(plotspecs, Ref(execution), plots_path, mod_name)
end
