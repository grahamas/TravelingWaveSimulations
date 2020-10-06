using Pkg
using GitCommand

devdir(dirs...) = joinpath(homedir(), ".julia", "dev", dirs...)
rel_path(dirs...) = joinpath(@__DIR__, dirs...)
function git_pull(git, path)
    run(`$git -C $(path) pull`)
end
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

git() do git
    git_pull.(Ref(git), all_paths)
end

#init_julia_git(pde_path)
for path in [s73_plot_path, s73_path, nm_path, tws_path, tws_plot_path]
    Pkg.activate(path)
    Pkg.develop("AxisIndices")
end
git() do git
    axisindices_path = joinpath(homedir(), ".julia", "dev", "AxisIndices")
    try
        run(`$git -C $axisindices_path remote add grahamas git@github.com:grahamas/AxisIndices.jl.git`)
    catch e
        @show e
    end
	run(`$git -C $axisindices_path fetch grahamas`)
	run(`$git -C $axisindices_path checkout logical_indexing`)
end
init_julia_git(s73_path)
init_julia_git(nm_path, [s73_path])
init_julia_git(wcm_path, [s73_path, nm_path])
init_julia_git(tws_path, [s73_path, nm_path, wcm_path])
init_julia_git(s73_plot_path, [s73_path])
init_julia_git(tws_plot_path, [s73_plot_path, s73_path, nm_path, wcm_path, tws_path])

init_julia_git(rel_path("_research"), [s73_path, nm_path, wcm_path, pde_path, tws_path, s73_plot_path, tws_plot_path])

Pkg.build()
Pkg.resolve()
