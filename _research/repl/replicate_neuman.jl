include_here(str) = include(joinpath(@__DIR__, str))

include_here("setup.jl")
neuman_exec, neuman_classifications = let neuman_example_name = "FAULTY_neuman_fft"
    neuman_exec = execute(get_example(neuman_example_name)())
    (heatmap_slices_scene, _layout) = heatmap_slices_execution(neuman_exec)
    Makie.save(plotsdir("replicated_neuman_heatmap_slices.png"), heatmap_slices_scene)
    neuman_exec, ExecutionClassifications(neuman_exec)
end
@warn "Set neuman_exec and neuman_classifications"
