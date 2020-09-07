using Pkg

rel_path(dirs...) = joinpath(@__DIR__, dirs...)
git_pull(path) = run(`git -C $(path) pull`)

s73_path = rel_path("..", "Simulation73")
nm_path = rel_path("..", "NeuralModels")
wcm_path = rel_path("..", "WilsonCowanModel")
pde_path = ""#rel_path("..", "FindPDE")
s73_plot_path = rel_path("..", "Simulation73Plotting")
tws_plot_path = rel_path("..", "TravelingWaveSimulationsPlotting")
tws_path = rel_path(".")

all_paths = filter(!=(""), [s73_path, nm_path, wcm_path, pde_path, tws_path, s73_plot_path, tws_plot_path])

git_pull.(all_paths)
map(all_paths) do path
    Pkg.activate(path)
    Pkg.update()
end

Pkg.activate(rel_path("_research"))
Pkg.update()

Pkg.instantiate()
Pkg.build()
