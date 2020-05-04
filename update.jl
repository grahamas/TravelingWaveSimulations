using Pkg

rel_path(dirs...) = joinpath(@__DIR__, dirs...)
git_pull(path) = run(`git -C $(path) pull`)

s73_path = rel_path("..", "Simulation73")
nm_path = rel_path("..", "NeuralModels")
wcm_path = rel_path("..", "WilsonCowanModel")
pde_path = rel_path("..", "FindPDE")
tws_path = rel_path(".")

git_pull.([s73_path, nm_path, wcm_path, pde_path, tws_path])

Pkg.activate(s73_path)
Pkg.update()
Pkg.activate(joinpath(s73_path, "test"))
Pkg.update()

Pkg.activate(nm_path)
Pkg.update()
Pkg.activate(joinpath(nm_path, "test"))
Pkg.update()

Pkg.activate(wcm_path)
Pkg.update()
Pkg.activate(joinpath(wcm_path, "test"))
Pkg.update()

Pkg.activate(pde_path)
Pkg.update()

Pkg.activate(tws_path)
Pkg.update()
Pkg.activate(joinpath(tws_path, "test"))
Pkg.update()

Pkg.instantiate()
Pkg.build()
