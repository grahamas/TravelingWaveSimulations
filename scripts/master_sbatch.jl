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
        "--job-name"
            help = "Name of slurm job"
            default = "julia-tws"
        "--output-subdir"
            help = "Subdir of _slurm/output for logs"
            default = ""
        "--script-name", "-s"
            help = "Name of script to run (path in script dir)"
        "--no-save-raw"
            help = "Don't save the raw simulation"
            action = :store_true
        "--nodes", "-N"
            help="minnodes[-maxnodes]"
        "--julia-exe"
            help="Julia to call"
            default="julia"
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
    output_subdir = pop!(args, "output-subdir")
    job_name = pop!(args, "job-name")
    script_name = pop!(args, "script-name")
    script_output_dir = joinpath(project_root, "_slurm", "output", output_subdir)
    mkpath(script_output_dir)
    sbatch_arg_names = ["nodes", "ntasks", "mem-per-cpu", "time", "partition", "workdir", "mail-user",  "mail-type"]
    sbatch_args = pop_args!(args, sbatch_arg_names)
    sbatch_args["output"] = joinpath(script_output_dir, "%j.stdout")
    sbatch_args["error"] = joinpath(script_output_dir, "%j.stderr")
    sbatch_args["job-name"] = job_name 

    script_args = args
    script_args["n-tasks"] = ntasks

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
