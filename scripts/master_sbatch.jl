using DrWatson
quickactivate(@__DIR__, "TravelingWaveSimulations")
using ArgParse


function parse_commandline(args)
    arg_settings = ArgParseSettings(
                                    autofix_names = false,
                                    )
    @add_arg_table! arg_settings begin
        "--ntasks"
            help = "Total number of tasks (cpus)"
            arg_type = Int
            default = 1
        "--mem-per-cpu"
            help = "RAM per CPU"
            default = "1G"
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
        "--mod"
            nargs = '*'
            help = "Case specifying base model modifications"
        "--no-save-raw"
            help = "Don't save the raw simulation"
            action = :store_true
        "--max-batch-size"
            help = "Parallel batch size"
            arg_type = Int
        "--nodes", "-N"
            help="minnodes[-maxnodes]"
        "--max-sims-in-mem"
            help="Limit the number of bare simulations held in memory by a single task"
        "--julia-exe"
            help="Julia to call"
            default="julia"
        "--backup-paths"
            nargs = '*'
            help = "Locations to copy saved results via scp"
    end
    return parse_args(args, arg_settings; as_symbols = false)
end

function pop_args!(args::AbstractDict, names::AbstractArray)
    Dict(filter(!isnothing, (map(names) do name
        val = pop!(args, name, nothing)
        if val != nothing
            return name => val
        end
    end))...)
end

arg_str(val::Number) = string(val)
arg_str(bool::Bool) = ""
arg_str(str::AbstractString) = str
arg_str(arr::AbstractArray) = "$(join(arr, " "))"
function arg_str_list(args::AbstractDict)
    hcat([["--$name", arg_str(val)] for (name, val) in args]...)
end
function single_arg_str(args::AbstractDict)
    join(arg_str_list(args), " ")
end

function sbatch_script(ARGS)
    args = parse_commandline(ARGS)

    julia_exe = pop!(args, "julia-exe")
    ntasks = args["ntasks"] # Needed below
    
    project_root = pop!(args, "project-root")
    base_example = pop!(args, "base-example")
    script_name = pop!(args, "script-name")
    script_output_dir = joinpath(project_root, "_slurm", "output", base_example)
    mkpath(script_output_dir)
    sbatch_arg_names = ["nodes", "ntasks", "mem-per-cpu", "time", "partition", "workdir", "mail-user",  "mail-type"]
    sbatch_args = pop_args!(args, sbatch_arg_names)
    sbatch_args["output"] = joinpath(script_output_dir, "%j.stdout")
    sbatch_args["error"] = joinpath(script_output_dir, "%j.stderr")
    sbatch_args["job-name"] = base_example

    script_arg_names = ["mod", "data-root", "max-batch-size", "max-sims-in-mem", "backup-paths"]
    script_args = pop_args!(args, script_arg_names)
    remaining_args = args
    @show remaining_args
    script_args["example-name"] = base_example
    if args["no-save-raw"]
        script_args["no-save-raw"] = args["no-save-raw"]
    end
    if ntasks > 1
        p_args = "-p $ntasks"
    else
        p_args = ""
    end

    script_path = joinpath(scriptsdir(), script_name)

    sbatch_script = """#!/bin/bash
    cd $(project_root)
    $(julia_exe) $p_args $(script_path) $(single_arg_str(script_args))
    cat $(joinpath(script_output_dir, "\${SLURM_JOB_ID}.stderr")) $(joinpath(script_output_dir, "\${SLURM_JOB_ID}.stdout")) | /usr/bin/mail -s \${SLURM_JOB_NAME} $(sbatch_args["mail-user"])
    """
    @show sbatch_script
    run(`echo "Starting job, using julia exec: $julia_exe"`)
    run(pipeline(`echo $sbatch_script`, `sbatch $(arg_str_list(sbatch_args))`))
end

sbatch_script(ARGS)
