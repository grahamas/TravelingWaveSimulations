particular_name, particular_exec = execute_single_modification(get_example(particular_example_name), (particular_mods..., other_opts=Dict(), step_reduction=nothing))
animate_execution(plotsdir("particular_test.mp4"), particular_exec)
particular_classifications = ExecutionClassifications(particular_exec)
