particular_name, particular_exec = execute_single_modification(get_example(particular_example_name), (particular_mods..., other_opts=Dict(), step_reduction=nothing))
(particular_heatmap_slices_scene, _layout) = heatmap_slices_exec()
Makie.save(plotsdir("particular_heatmap_slices.png"), particular_heatmap_slices_scene)
particular_classifications = ExecutionClassifications(particular_exec)
