using DrWatson
quickactivate(@__DIR__, "WilsonCowanModel")
using ArgParse

function parse_commandline()
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
        "--plotspec-case"
            help = "Case specifying plots"
        "--modifications-case"
            help = "Case specifying base model modifications"
        "--no-save-raw"
            help = "Don't save the raw simulation"
            action = :store_true
    end
    return parse_args(arg_settings; as_symbols = false)
end

function main()
    args = parse_commandline()
    project_root = args["project-root"]
    base_example = args["base-example"]
    script_name = args["script-name"]
    sbatch_arg_names = ["ntasks", "mem-per-cpu", "time", "partition", "workdir", "mail-user",  "mail-type"]
    script_output_dir = joinpath(project_root, "_slurm", "output", base_example)
    mkpath(script_output_dir)
    sbatch_args = ["--$name=$val" for (name, val) in args if (name in sbatch_arg_names && val != nothing)]
    sbatch_args = [sbatch_args..., """--output=$(joinpath(script_output_dir, "%j.stdout"))"""]
    sbatch_args = [sbatch_args..., """--error=$(joinpath(script_output_dir, "%j.stderr"))"""]
    sbatch_args = [sbatch_args..., "--job-name=$(base_example)"]

    script_arg_names = ["modifications-case", "plotspec-case", "data-root"]
    script_args = ["--$name=$val" for (name, val) in args if (name in script_arg_names && val != nothing)]
    script_args = [script_args..., """--example-name=$(args["base-example"])"""]
    if args["no-save-raw"]
        script_args = [script_args..., "--no-save-raw"]
    end

    script_path = joinpath(scriptdir(), script_name)

    sbatch_script = """#!/bin/bash
    cd $(project_root)
    julia $(script_path) $(join(script_args, " "))
    /usr/bin/mail -s \${SLURM_JOB_NAME} $(args["mail-user"]) < $(joinpath(script_output_dir, "\${SLURM_JOB_ID}.stderr"))
    """

    run(pipeline(`echo $sbatch_script`, `sbatch $sbatch_args`))
end

main()
