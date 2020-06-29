include_here(str) = include(joinpath(@__DIR__, str))

include_here("setup.jl")
let neuman_example_name = "oscillating_pulse"
    mkpath(plotsdir("tau_sweep"))
    for tau_i = [1.0, 5.0, 8.0, 10.0, 12.0, 18.0, 20.0, 30.0, 50.0, 100.0]
        GC.gc()
        neuman_exec = execute(get_example(neuman_example_name)(; τ=(10.0,tau_i)))
        (heatmap_slices_scene, _layout) = heatmap_slices_execution(neuman_exec)
        Makie.save(plotsdir("tau_sweep/tau_$(tau_i)_neuman_heatmap_slices.png"), heatmap_slices_scene)
    end
end
