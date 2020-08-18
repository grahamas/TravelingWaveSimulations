using Documenter, Simulation73
DocMeta.setdocmeta!(Simulation73, :DocTestSetup, :(using Simulation73); recursive=true)
makedocs(
    sitename="Simulation73",
    modules=[Simulation73],
    pages = [
       "Home" => "index.md",
       "Manual" => Any[
           "man/simulation.md"
       ],
       "Library" => Any[
           "Public" => "lib/public.md",
           hide("Internals" => "lib/internals.md",
            Any[joinpath("lib/internals", f) for f in readdir("docs/src/lib/internals")]
           )
       ]
   ]
)

deploydocs(
    repo = "github.com/grahamas/Simulation73.git"
)
