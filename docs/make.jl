using Documenter, TravelingWaveSimulations

DocMeta.setdocmeta!(TravelingWaveSimulations, :DocTestSetup, :(using TravelingWaveSimulations); recursive=true)

makedocs(
    sitename="TravelingWaveSimulations",
    modules=[TravelingWaveSimulations],
)

deploydocs(
    repo = "github.com/grahamas/TravelingWaveSimulations.git"
)
