using Pkg

rel_path(dirs...) = joinpath(@__DIR__, dirs...)
git_pull(path) = run(`git -C $(path) pull`)

s73_path = rel_path("..", "Simulation73")
nm_path = rel_path("..", "NeuralModels")
wcm_path = rel_path("..", "WilsonCowanModel")
pde_path = rel_path("..", "FindPDE")
tws_path = rel_path(".")

for path in [s73_path, nm_path, wcm_path, pde_path]
	if !ispath(path)
		run(`git clone git@github.com:grahamas/$(basename(path)).git $(path)`)
	end
end

git_pull.([s73_path, nm_path, wcm_path, pde_path, tws_path])

Pkg.activate(pde_path)
Pkg.update()

Pkg.activate(s73_path)
Pkg.update()

Pkg.activate(nm_path)
Pkg.develop(PackageSpec(path=s73_path))
Pkg.update()

Pkg.activate(wcm_path)
Pkg.develop([PackageSpec(path=x) for x in [s73_path, nm_path]])
Pkg.update()

Pkg.activate(tws_path)
Pkg.develop([PackageSpec(path=x) for x in [s73_path, nm_path, wcm_path, pde_path]])
Pkg.update()

Pkg.activate(rel_path("_research"))
Pkg.develop([PackageSpec(path=x) for x in [s73_path, nm_path, wcm_path, pde_path, tws_path]])
Pkg.update()

Pkg.build()
