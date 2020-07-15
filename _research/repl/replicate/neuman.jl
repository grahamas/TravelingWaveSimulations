include_here(str...) = include(joinpath(@__DIR__, str...))

include_here("..","setup","basic.jl")
include_here("..","setup","plot.jl")
neuman_exec, neuman_classifications = let neuman_replication_name = "faulty_neuman_fft"
    neuman_exec = execute(TravelingWaveSimulations.replications_dict[neuman_replication_name]())
    (heatmap_slices_scene, _layout) = heatmap_slices_execution(neuman_exec)
    savepath = plotsdir("replicated_neuman_heatmap_slices.png")
    @warn "saving $savepath"
    Makie.save(savepath, heatmap_slices_scene)
    neuman_exec, ExecutionClassifications(neuman_exec)
end
@warn "Set neuman_exec and neuman_classifications"
