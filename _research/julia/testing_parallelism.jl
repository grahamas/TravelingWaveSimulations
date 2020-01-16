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
using TravelingWaveSimulations, Simulation73, Plots, JuliaDB, DifferentialEquations, Dates, DrWatson, NeuralModels, WilsonCowanModel
using DiffEqBase: AbstractTimeseriesSolution

# %%
##############

### Difference of Sigmoids
struct DifferenceOfSigmoids{T} <: AbstractNonlinearity{T}
    firing_sigmoid::SigmoidNonlinearity{T}
    blocking_sigmoid::SigmoidNonlinearity{T}
    DifferenceOfSigmoids(fsig::SigmoidNonlinearity{T},bsig::SigmoidNonlinearity{T}) where T = new{T}(fsig,bsig)
end
DifferenceOfSigmoids(; firing_a, firing_θ, blocking_a, blocking_θ) where T = DifferenceOfSigmoids(SigmoidNonlinearity(; θ=firing_θ,a=firing_a), SigmoidNonlinearity(; θ=blocking_θ,a=blocking_a))
(dos::DifferenceOfSigmoids)(output::AbstractArray, ignored_source, ignored_t) = output .= NeuralModels.rectified_sigmoid_fn.(output, dos.firing_sigmoid.a, dos.firing_sigmoid.θ) - NeuralModels.rectified_sigmoid_fn.(output, dos.blocking_sigmoid.a, dos.blocking_sigmoid.θ)

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
dos_example = TravelingWaveSimulations.@EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; SNR_scale=80.0, stop_time=ABS_STOP,
                                                     Aee=280.0, See=70.0,
                                                     Aii=1.4, Sii=70.0,
                                                     Aie=270.0, Sie=90.0,
                                                     Aei=-297.0, Sei=90.0,
                                                     n=71, x=500.0)
  simulation = Simulation(
    WCMSpatial(;
      pop_names = ("E", "I"),
      α = (1.0, 1.0),
      β = (1.0, 1.0),
      τ = (10.0, 10.0),
      nonlinearity = pops(DifferenceOfSigmoids;
        #sd = [6.7, sqrt(3.2)],
        #θ = [18.0, 10.0]),
        firing_θ = [10.0, 5.0],
        firing_a = [1.2, 1.0],
        blocking_θ = [18.0, 13.0],
        blocking_a = [1.2, 1.0]),
      stimulus = pops(SharpBumpStimulusParameter;
          strength = [10.0, 0.0],
          width = [28.1, 28.1],
          time_windows = [[(0.0, 10.0)], [(0.0, 10.0)]],
          baseline=[0.0, 0.0]),
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
                    (all(isapprox.(sub_u, 0.0, atol=1e-4)) || (sub_u[end] > 0.01)) && t > 5
                end
    else
        (u,t,integrator) -> begin
                    pop = population(u,1)
                    (all(isapprox.(u, 0.0, atol=1e-4)) || (sum(pop[:,end]) / size(pop,1) > 0.01)) && t > 5
            end
    end, terminate!)
  )
end



# %%
execution = execute(dos_example(; n=256, x=700.0, See=25.0, Sii=25.0, Sie=27.0, Sei=27.0,
                                Aee=200.0, Aei=75.0, Aie=50.0, Aii=10.0, strengthE=10.0, widthE=50.0));
anim = TravelingWaveSimulations.custom_animate(execution)
mp4(anim, "tmp/dos_tmp.mp4")

# %%

function point_timecourse(exec::Execution, location_proportion)
    us = exec.solution.u
    t = timepoints(exec)
    xs = space(exec)
    n_xs = size(xs)
    x_dx = floor.(Int, n_xs .* location_proportion)
    @show x_dx
    x = xs[x_dx...]
    u = [[u[x_dx...]] for u in us]
    return BareSolution(u=u,x=[x],t=t)    
end
function point_average(u::AbstractArray, pt::Tuple, xs, σ)
    pt_dist = map(xs) do x
        -sum((x .- pt) ./ σ) .^ 2) / 2
    end
    unscaled = exp.(pt_dist)
    return unscaled ./ sum(unscaled)
end

function point_average_timecourse(exec::Execution, location_proportion, σ)
    us = exec.solution.u
    t = timepoints(exec)
    xs = space(exec)
    n_xs = size(xs)
    x_dx = floor.(Int, n_xs .* location_proportion)
    @show x_dx
    x = xs[x_dx...]
    @show x
    u = [[point_average(u, x, xs.arr, σ)] for u in us]
    return BareSolution(u=u,x=[x],t=t)  
end
Base.getindex(lat::CompactLattice, dx) = getindex(lat.arr, dx)
function Base.getproperty(bs::BareSolution, sym::Symbol)
    if sym ∈ [:dense, :prob]
        return false
    elseif sym == :tslocation
        return 0
    else
        return Base.getfield(bs, sym)
    end
end



# %%
pt_over_time = point_timecourse(execution, 0.2)
plot(pt_over_time, labels=:E) |> display
pt_avg_over_time(point_average_timecourse(execution, 0.5, 25))
plot(pt_avg_over_time, labels=:E) |> display

# %%
aees = 1.0:3.0:40.0
maxes = [maximum(execute(my_example(; n=256, x=700.0, See=25.0, Sii=25.0, Sie=27.0, Sei=27.0,
                                Aee=Aee, Aei=0.0, Aie=0.0, Aii=0.1, strengthE=6.0)).solution) for Aee in aees]

# %%
plot(aees, maxes)

# %%

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
