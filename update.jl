using Pkg

rel_path(dirs...) = joinpath(@__DIR__, dirs...)
git_pull(path) = run(`git pull --git-dir=$(path)`)

s73_path = rel_path("..", "Simulation73")
nm_path = rel_path("..", "NeuralModels")
wcm_path = rel_path("..", "WilsonCowanModel")
tws_path = rel_path(".")

git_pull.(s73_path, nm_path, wcm_path, tws_path)

Pkg.activate(s73_path)
Pkg.update()

Pkg.activate(nm_path)
Pkg.develop(PackageSpec(path=s73_path))
Pkg.update()

Pkg.activate(wcm_path)
Pkg.develop([PackageSpec(path=x) for x in [s73_path, nm_path]])
Pkg.update()

Pkg.activate(tws_path)
Pkg.develop([PackageSpec(path=x) for x in [s73_path, nm_path, wcm_path]])
Pkg.update()

Pkg.build()
