using Documenter, TravelingWaveSimulations

makedocs(
    sitename="TravelingWaveSimulations",
    modules=[TravelingWaveSimulations],
)

deploydocs(
    repo = "github.com/grahamas/TravelingWaveSimulations.git"
)
