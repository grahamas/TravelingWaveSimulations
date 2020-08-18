# Simulation

```@meta
DocTestSetup = quote
    using Simulation73
end
```


Run a simulation from a `.jl` specification file that defines [`Simulation`](@ref), [`Output`](@ref), and [`Analysis`](@ref) objects. Call [`run_simulation`](@ref) on the specification file, which will `include` the specification file, copy the file to the output directory for reproducibility, [`execute`](@ref) the simulation, and run the analyses.

A `Simulation` is specified as

```
simulation = Simulation(;
  model = ...,
  space = ...,
  tspan = ...,
  initial_value = ...,
  solver = ...
  )
```
The first parameter, `model` specifies the differential equation to be solved, and the remainder specify the solving environment:
- `space` defines the space of each solution point (a.k.a. mesh. For typical ODE this may be unnecessary).
- `tspan` is a `Tuple{T,T}` defining the start and end time of the solution.
- `initial_value` must be of the same size as `space`.
- `solver` is one of the solvers from DifferentialEquations or as accepted by a similar imported `solve` function.

## Specifying a Model

Begin any simulation by defining a model, say of type `ToyModel <: AbstractModel`. The model must extend [`make_system_mutator`](@ref) to return the differential function to be solved.

To implement parameter exploration, the model must additionally extend [`update_from_p!`](@ref).
