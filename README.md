# TravelingWaveSimulations
Ensemble simulations of neural activity with analyses of traveling wave behavior (Ph.D. thesis project)

## Installation

Clone repository and run `setup.jl` in julia, e.g. `julia -e 'include("setup.jl")'` from the command line with working directory at the top level of this repo.

For visualising the simulations, see [grahamas/TravelingWavePlotting](https://github.com/grahamas/TravelingWaveSimulationsPlotting)


## Usage

See scripts in `_research/repl` for examples. Note that they are intended to be run within the julia REPL and mostly require a global variable `data_root` to be set.

Typically I use the function `iterate_prototype` which takes the name of a prototype defined in `src/prototypes/prototypes.jl` and a specification of modifications to be made to that prototype's parameters.

`iterate_prototype` will automatically parallelize multiple simulations specified by an array of modifications over workers or threads, if available. Does not currently support parallelizing over both, though the DifferentialEquations solver should use multithreading.
