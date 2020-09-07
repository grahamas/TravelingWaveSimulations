using Pkg

rel_path(dirs...) = joinpath(@__DIR__, dirs...)
git_pull(path) = run(`git -C $(path) pull`)
function init_julia_git(path)
    Pkg.activate(path)
    Pkg.update()
end
function init_julia_git(path, dev_dep_paths)
    Pkg.activate(path)
    Pkg.develop([PackageSpec(path=x) for x in dev_dep_paths if x != ""])
    Pkg.update()
end

s73_path = rel_path("..", "Simulation73")
nm_path = rel_path("..", "NeuralModels")
wcm_path = rel_path("..", "WilsonCowanModel")
pde_path = ""# rel_path("..", "FindPDE")
s73_plot_path = rel_path("..", "Simulation73Plotting")
tws_plot_path = rel_path("..", "TravelingWaveSimulationsPlotting")
tws_path = rel_path(".")

all_paths = filter(!=(""), [s73_path, nm_path, wcm_path, pde_path, tws_path, s73_plot_path, tws_plot_path])

for path in all_paths
	if !ispath(path)
		run(`git clone git@github.com:grahamas/$(basename(path)).git $(path)`)
	end
end

git_pull.(all_paths)

#init_julia_git(pde_path)
init_julia_git(s73_path)
init_julia_git(nm_path, [s73_path])
init_julia_git(wcm_path, [s73_path, nm_path])
init_julia_git(tws_path, [s73_path, nm_path, wcm_path])
init_julia_git(s73_plot_path, [s73_path])
init_julia_git(tws_plot_path, [s73_plot_path, s73_path, nm_path, wcm_path, tws_path])
init_julia_git(rel_path("_research"), [s73_path, nm_path, wcm_path, pde_path, tws_path, s73_plot_path, tws_plot_path])

Pkg.build()
