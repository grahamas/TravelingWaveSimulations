# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.3.1
#     language: julia
#     name: julia-1.3
# ---

# %%
using TravelingWaveSimulations, Simulation73, Plots, JuliaDB, DifferentialEquations, Revise, Dates, DrWatson, NeuralModels, WilsonCowanModel
using DiffEqBase: AbstractTimeseriesSolution

# %%


function most_recent_subdir(datadir)
    subdirs = joinpath.(Ref(datadir), readdir(datadir))
    sort(subdirs, by=mtime, rev=true)[1]
end

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

function mod_example(example; data_root, no_save_raw=false, example_name, modifications, analyses, batch, max_sims_in_mem)
    modifications, modifications_prefix = parse_modifications_argument(modifications)
    analyses, analyses_prefix = parse_analyses_argument(analyses)
    @show analyses

    # Initialize saving paths
    if length(analyses) > 0
        analyses_path = joinpath(plotsdir(), example_name, "$(modifications_prefix)$(analyses_prefix)_$(Dates.now())_$(gitdescribe())")
        mkpath(analyses_path)
    else
        analyses_path = ""
    end
    if !no_save_raw
        data_path = joinpath(data_root, example_name, "$(modifications_prefix)$(Dates.now())_$(gitdescribe())")
    else
        data_path = nothing
    end
    
    if !(@isdefined nprocs) || nprocs() == 1
        @warn "Not parallelizing parameter sweep."
        failures = TravelingWaveSimulations.execute_modifications_serial(example, modifications, analyses, data_path, analyses_path, no_save_raw)
        @show failures
    else
        @warn "Parallelizing parameter sweep."
        if no_save_raw
            failures = TravelingWaveSimulations.execute_modifications_parallel_nosaving(example, modifications, analyses, data_path, analyses_path, batch)
        else
            failures = TravelingWaveSimulations.execute_modifications_parallel_saving(example, modifications, analyses, data_path, analyses_path, batch, max_sims_in_mem)
        end
        @show failures
    end
end



# %%
# Aee: E->E amplitude
# See: E->E spread

# caS: cross-auto spread ratio
# caA: cross-auto amplitude ratio

# ioA: inhibitory output amplitude scale
# ioS: inhibitory output spread scale
# iiA: inhibitory input amplitude scale
# iiS: inhibitory input spread scale
ABS_STOP=300.0
my_example = TravelingWaveSimulations.@EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; SNR_scale=80.0, stop_time=ABS_STOP,
                                                     Aee=280.0, See=70.0,
                                                     Aii=1.4, Sii=70.0,
                                                     Aie=270.0, Sie=90.0,
                                                     Aei=-297.0, Sei=90.0,
                                                     n=71, x=500.0)
  simulation = Simulation(
    WCMSpatial(;
      pop_names = ("E", "I"),
      α = (1.0, 1.0),
      β = (1.1, 1.1),
      τ = (10.0, 10.0),
      nonlinearity = pops(GaussianNonlinearity;
        sd = [6.7, sqrt(3.2)],
        θ = [18.0, 10.0]),
      stimulus = pops([SharpBumpStimulusParameter(;
          strength = 10.0,
          width = 100.0,
          time_windows = [(0.0, 10.0)]),
          NoStimulusParameter{Float64}()]),
      connectivity = FFTParameter(pops(GaussianConnectivityParameter;
          amplitude = [Aee -Aei;
                       Aie -Aii],
          spread = [(See,See) (Sei,Sei);
                    (Sie,Sie) (Sii,Sii)]
         ))
      );
      space = PeriodicLattice{Float64,N_ARR}(; n_points=(n,n), extent=(x,x)),
      save_idxs = RadialSlice(),
      tspan = (0.0,stop_time),
      dt = 1.0,
      algorithm=Euler(),
      callback=DiscreteCallback(if !(save_idxs === nothing)
        (u,t,integrator) -> begin
                    sub_u = u[integrator.opts.save_idxs];
                    (all(sub_u .≈ 0.0) || (sub_u[end] > 0.01)) && t > 5
                end
    else
        (u,t,integrator) -> begin
                    pop = population(u,1)
                    (all(u .≈ 0.0) || (sum(pop[:,end]) / size(pop,1) > 0.01)) && t > 5
            end
    end, terminate!)
  )
end




# %%
execution = execute(my_example());

# %% collapsed=true jupyter={"outputs_hidden": true}
length(execution.solution)
@show execution.solution.t

# %% collapsed=true jupyter={"outputs_hidden": true}
anim = TravelingWaveSimulations.custom_animate(execution)
mp4(anim, "tmp/_tmp.mp4")

# %%
N_ATT = 128; SAVE=RadialSlice()
failures = mod_example(example; example_name="sigmoid_normal_fft", data_root="data", modifications=["See=1.0","n=$N_ATT", "save_idxs=$SAVE"], analyses=[], batch=5, max_sims_in_mem=nothing)

# %% jupyter={"source_hidden": true}
using Distributed;  @show nworkers(); addprocs(1); @show nworkers()
@everywhere using Pkg
@everywhere Pkg.activate("..")
@everywhere using TravelingWaveSimulations
@show nworkers()


# %% collapsed=true jupyter={"outputs_hidden": true, "source_hidden": true}
# small parallel example
failures = based_on_example(; example_name="sigmoid_normal_fft", data_root="data", modifications=["See=1.0:10.0", "stop_time=90.0"], batch=20)

# %%
recent_dir = most_recent_subdir("data/sigmoid_normal_fft")
d1 = JuliaDB.load(joinpath(recent_dir, "raw_data.csv"))

# %%
param_db = select(d1, Keys())
soln_db = select(d1, Not(Keys()))
@show map(length, select(d1, :t))
@show map((x) -> x[end], select(d1, :t))
exec = Execution(example(; save_idxs=SAVE, n=N_ATT, param_db[1]...), BareSolution(; soln_db[1]...))
@recipe function f(lattice::AbstractLattice{T,2}, values::Array{T,2}; val_lim=nothing) where T
    (x, y) = coordinate_axes(lattice) .|> collect
    seriestype := :heatmap
    if val_lim != nothing
        clim := val_lim
        zlim := val_lim
    end
    aspect_ratio := :equal
    (x,y,values)
end
plot(space(exec), population_timepoint(exec.solution, 1, 1)) |>  display
@show population(exec.solution.u[1], 1) |> size

# %%
for n_row in 1:length(d1)
    anim = custom_animate(Execution(example(;  param_db[n_row]..., n=N_ATT, save_idxs=SAVE), BareSolution(; soln_db[n_row]...)), title=(param_db[n_row] |> string))
    mp4(anim, "tmp/$(n_row)_tmp.mp4") |> display
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
