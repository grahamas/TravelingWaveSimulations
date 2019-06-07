using DrWatson
quickactivate(@__DIR__, "TravelingWaveSimulations")
using ArgParse


function parse_commandline(args)
    arg_settings = ArgParseSettings(
                                    autofix_names = false,
                                    )
    @add_arg_table arg_settings begin
        "--ntasks"
            help = "Total number of tasks (cpus)"
            arg_type = Int
            default = 1
        "--mem-per-cpu"
            help = "RAM per CPU (MB)"
            arg_type = Int
            default = 1000
        "--time"
            help = "Job duration (seconds)"
            default = "30:00"
        "--partition"
            help = "Name of partition to run job on"
            default = "debug"
        "--project-root"
            help = "Name of code directory"
            default = projectdir()
        "--data-root"
            help = "Name of data output directory"
            default = datadir()
        "--mail-user"
            help = "Email target for SLURM status updates"
        "--mail-type"
            help = "Types of mail to send user"
        "--base-example", "-b"
            help = "Name of predefined example to use as base"
        "--script-name", "-s"
            help = "Name of script to run (path in script dir)"
            default = "based_on_example.jl"
        "--plot"
            help = "Case specifying plots"
        "--mod"
            nargs = '*'
            help = "Case specifying base model modifications"
        "--no-save-raw"
            help = "Don't save the raw simulation"
            action = :store_true
    end
    return parse_args(args, arg_settings; as_symbols = false)
end

function pop_args!(args::AbstractDict, names::AbstractArray)
    Dict((map(names) do name
        val = pop!(args, name, nothing)
        if val != nothing
            return name => val
        end
    end)...)
end

arg_str(str::AbstractString) = str
arg_str(arr::AbstractArray) = "$(arr...)"
function arg_str(args::AbstractDict)
    arg_str_list = ["--$name $(arg_str(val))" for (name, val) in args]
    return join(arg_str_list, " ")
end

function sbatch_script(ARGS)
    args = parse_commandline(ARGS)

    project_root = pop!(args, "project-root")
    base_example = pop!(args, "base-example")
    script_name = pop!(args, "script-name")
    script_output_dir = joinpath(project_root, "_slurm", "output", base_example)
    mkpath(script_output_dir)
    sbatch_arg_names = ["ntasks", "mem-per-cpu", "time", "partition", "workdir", "mail-user",  "mail-type"]
    sbatch_args = pop_args!(args, sbatch_arg_names)
    sbatch_args["output"] = joinpath(script_output_dir, "%j.stdout")
    sbatch_args["error"] = joinpath(script_output_dir, "%j.stderr")
    sbatch_args["job-name"] = base_example

    script_arg_names = ["mod", "plot", "data-root"]
    script_args = pop_args!(args, script_arg_names)
    script_args["example-name"] = args["base-example"]
    if args["no-save-raw"]
        script_args["no-save-raw"] = args["no-save-raw"]
    end

    script_path = joinpath(scriptdir(), script_name)

    sbatch_script = """#!/bin/bash
    cd $(project_root)
    julia-nightly $(script_path) $(arg_str(script_args))
    /usr/bin/mail -s \${SLURM_JOB_NAME} $(args["mail-user"]) < $(joinpath(script_output_dir, "\${SLURM_JOB_ID}.stderr"))
    """

    run(pipeline(`echo $sbatch_script`, `sbatch $(arg_str(sbatch_args))`))
end

sbatch_script(ARGS)
