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
        "--workdir"
            help = "Name of working directory"
            default = "."
        "--mail-user"
            help = "Email target for SLURM status updates"
        "--mail-type"
            help = "Types of mail to send user"
        "--base-example", "-b"
            help = "Name of predefined example to use as base"
    end
    return parse_args(arg_settings; as_symbols = false)
end

function main()
    args = parse_commandline()
    @show args
    sbatch_arg_names = ["ntasks", "mem-per-cpu", "time", "partition", "workdir", "mail-user",  "mail-type"]
    script_output_dir = joinpath(args["workdir"], "_slurm", "output", args["base-example"])
    mkpath(script_output_dir)
    sbatch_args = ["--$name=$val" for (name, val) in args if (name in sbatch_arg_names && val != nothing)]
    sbatch_args = [sbatch_args..., """--output=$(joinpath(script_output_dir, "%j.%N.stdout"))"""]
    sbatch_args = [sbatch_args..., """--error=$(joinpath(script_output_dir, "%j.%N.stderr"))"""]
    sbatch_args = [sbatch_args..., "--job-name=$(args["base-example"])"]

    example_script = """#!/bin/bash
    julia scripts/based_on_example.jl $(args["base-example"])
    """

    run(pipeline(`echo $example_script`, `sbatch $sbatch_args`))
end

main()

    
