# TravelingWaveSimulations
Ensemble simulations of neural activity with analyses of traveling wave behavior (Ph.D. thesis project)

## Installation

```julia
using Pkg
pkg"registry add https://github.com/grahamas/SmithThesisRegistry.git"
pkg"add TravelingWaveSimulations"
```

For visualising the simulations, see [grahamas/TravelingWavePlotting](https://github.com/grahamas/TravelingWaveSimulationsPlotting)


## Usage

See scripts in `_research/repl` for examples. Note that they are intended to be run within the julia REPL and mostly require a global variable `data_root` to be set.

Typically I use the function `iterate_prototype` which takes the name of a prototype defined in `src/prototypes/prototypes.jl` and a specification of modifications to be made to that prototype's parameters.

`iterate_prototype` will automatically parallelize multiple simulations specified by an array of modifications over workers or threads, if available. Does not currently support parallelizing over both, though the DifferentialEquations solver should use multithreading.
