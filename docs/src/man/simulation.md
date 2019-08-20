# Simulation

Run a simulation from a `.jl` specification file that defines [`Simulation`](@ref), [`Output`](@ref), and [`Analysis`](@ref) objects. Call [`run_simulation`](@ref) on the specification file, which will `include` the specification file, copy the file to the output directory for reproducibility, [`execute`](@ref) the simulation, and run the analyses.

A `Simulation` is specified as

```
simulation = Simulation(;
  model = ...,
  solver = ...
  )
```
The `solver` parameter is of type [`Solver`](@ref) as defined by this module, whereas the `model` parameter is a user-defined subtype of [`AbstractModel`](@ref).

`solver` specifies the argument to the differential equation solver.

`model` specifies the differential equation to solve.

## Specifying a Model

Begin any simulation by defining a model, say of type `ToyModel <: AbstractModel`. Minimally, if `model isa ToyModel` then `model.space isa AbstractSpace` and `initial_value(model::ToyModel)` must be defined. Finally, the model must extend [`make_system_mutator`](@ref) to return the differential function to be solved.

To implement parameter exploration, the model must additionally extend [`update_from_p!`](@ref).
