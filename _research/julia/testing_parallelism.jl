# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all
#     notebook_metadata_filter: all,-language_info
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Julia 1.3.0
#     language: julia
#     name: julia-1.3
# ---

# %%
using TravelingWaveSimulations, Simulation73, Plots, JuliaDB, DifferentialEquations
using DiffEqBase: AbstractTimeseriesSolution

# %%
example=TravelingWaveSimulations.examples_dict["sigmoid_normal_fft"]

# %%
execution = execute(example(n=512, x=700.0, amplitude=([25.0 -25.2; 35.0 -4.0]), stop_time=30.0));

# %%
anim = TravelingWaveSimulations.custom_animate(execution)
mp4(anim, "tmp/_tmp.mp4")

# %%
failures = based_on_example(; example_name="sigmoid_normal_fft", data_root="data", modifications=["See=1.0:10.0"], batch=5)

# %%
using Distributed;  @show nworkers(); addprocs(1); @show nworkers()
@everywhere using Pkg
@everywhere Pkg.activate("..")
@everywhere using TravelingWaveSimulations
@show nworkers()


# %%
# small parallel example
failures = based_on_example(; example_name="sigmoid_normal_fft", data_root="data", modifications=["See=1.0:10.0"], batch=5)

# %%
function most_recent_subdir(datadir)
    subdirs = joinpath.(Ref(datadir), readdir(datadir))
    sort(subdirs, by=mtime, rev=true)[1]
end

# %%
# abstract type AbstractSolution{S} end
# const SomeSolution{S} = Union{AbstractTimeseriesSolution{S}, AbstractSolution{S}}
struct BareSolution{S,N,U<:Array{<:Array{<:S,N}},X,T} <: AbstractTimeseriesSolution{S,N}
    u::U
    x::X
    t::T
end
BareSolution(; u::U, x::X, t::T) where {S,N,U<:Array{<:Array{<:S,N}},X,T} = BareSolution{S,N,U,X,T}(u,x,t)
timepoints(bs::BareSolution) = bs.t
coordinates(bs::BareSolution) = bs.x
function Base.show(io::IO, A::BareSolution)
  print(io,"t: ")
  show(io, A.t)
  println(io)
  print(io,"u: ")
  show(io, A.u)
end
function Base.show(io::IO, m::MIME"text/plain", A::BareSolution)
  print(io,"t: ")
  show(io,m,A.t)
  println(io)
  print(io,"u: ")
  show(io,m,A.u)
end
@generated function Simulation73.population_timepoint(solution::BareSolution{T,NP}, pop_dx::Int, time_dx::Int) where {T,NP}
    N = NP - 1 # N + pops
    colons = [:(:) for i in 1:N]
    :(solution[time_dx][$(colons...), pop_dx])
end
function TravelingWaveSimulations.custom_animate(execution::Execution{T,<:Simulation{T},<:BareSolution}; kwargs...) where T
    solution = execution.solution
    pop_names = execution.simulation.model.pop_names
    x = space(execution)
    t = Simulation73.timepoints(execution)
    max_val = maximum(solution)
    min_val = minimum(solution)
    @animate for time_dx in 1:length(t) # TODO @views
        plot([
                plot(
                    x, population_timepoint(solution, 1, time_dx); label=pop_names[1],
                    val_lim=(min_val,max_val), title="t = $(round(t[time_dx], digits=4))",
                    xlab = "Space (a.u. approx. um)",kwargs...
                    )
                for i_pop in 1:length(pop_names)
            ]...)
    end
end

# %%
recent_dir = most_recent_subdir("data/sigmoid_normal_fft")
d1 = load(joinpath(recent_dir, "1.jdb"))
d2 = load(joinpath(recent_dir, "2.jdb"))

# %% jupyter={"outputs_hidden": true} collapsed=true
recent_dir = most_recent_subdir("$(ENV["HOME"])/sims")
recent_db = load(joinpath(recent_dir, "1.jdb"))

# %%
recent_dir = most_recent_subdir("$(ENV["HOME"])/sim_data")
@show recent_dir
recent_db = load(joinpath(recent_dir, "1.jdb"))

# %%
all_sols = select(recent_db, Not(Keys()))
all_keys = select(recent_db, Keys());

# %%
idx = 3000
bs = BareSolution(;all_sols[idx]...);
this_key = all_keys[idx]
ex = Execution(example(; this_key...), bs)
mp4(custom_animate(ex; title = string(this_key)), "tmp.mp4")

# %%
longer_execution = execute(example(; stop_time=ex.simulation.tspan[2]*2, this_key...))
mp4(custom_animate(longer_execution, title=string(this_key)), "tmp_longer.mp4")
