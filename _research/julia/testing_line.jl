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
using Revise
using Simulation73, NeuralModels, TravelingWaveSimulations, Plots, Optim, LinearAlgebra, Distances, Statistics,
    IterTools, Combinatorics, DataFrames, GLM, JuliaDB, DifferentialEquations, WilsonCowanModel

equals_str(key,val) = "$key=$val"
equals_strs(mods) = [equals_str(p...) for p in pairs(mods)]
mods_filename(x) = join(equals_strs(x), "_")

function velocity_results(results)
    (coef(results[:apex_loc])[1], deviance(results[:apex_loc]))
end

function mean_skip_missing(A; dims)
    missings = ismissing.(A)
    zeroed = copy(A)
    zeroed[missings] .= 0
    nonmissingsum = sum(zeroed; dims=dims)
    nonmissingmean = nonmissingsum ./ sum(.!missings; dims=dims)
    return nonmissingmean
end

function find_first_satisfying_execution(mdb, example, dict_min=Dict(), dict_max=Dict())
    function filter_fn(row)
        above_mins = [row[key] >= val for (key, val) in pairs(dict_min)]
        below_maxes = [row[key] <= val for (key, val) in pairs(dict_max)]
        return all(above_mins) && all(below_maxes)
    end     
    for db in mdb
        satisfactory_rows = filter(filter_fn, db)
        if length(satisfactory_rows) > 0
            mods = JuliaDB.select(satisfactory_rows, Keys())
            @show "found $(mods[1])"
            sols = JuliaDB.select(satisfactory_rows, JuliaDB.Not(Keys()))
            return (mods[1], Execution(example(;mods[1]...), BareSolution(; pairs(sols[1])...)))
        end
    end
    return nothing    
end

# %%
function is_traveling(l_frame_fronts::Array{<:Array,1}, t, min_motion)
    rightmost_locations = [frame_fronts[end].slope.loc for frame_fronts in l_frame_fronts]
    rightward_motion = diff(rightmost_locations) ./ diff(t)
    return all(rightward_motion[(length(rightward_motion) ÷ 2):end-5] .> min_motion)
end

function persistent_activation(l_frames, t, min_activation)
    # Test if the left-most value is ever high and non-decreasing
    leftmost_value = [frame[1] for frame in l_frames[(length(t) ÷ 2):end]]
    d_leftmost_value = diff(leftmost_value) ./ diff(t[(length(t) ÷ 2):end])
    # test leftmost is not ONLY decreasing
    low_enough = leftmost_value[2:end] .< min_activation
    decreasing = d_leftmost_value .< 0.0
    return !all(low_enough .| decreasing)
end



function is_epileptic_radial_slice(exec, min_motion=1e-2, min_activation=5e-2)
    # We'll call it epileptic if:
    #  1. There is a traveling front
    #  2. There is persistently elevated activity following the front.
    # How to deal with oscillatory activity? 
    # Want persistent oscillation to be considered epileptic, even if sometimes zero
    # In practice we'll ignore this for now --- FIXME
    l_frames = exec.solution.u
    t = timepoints(exec)
    x = [x[1] for x in space(exec).arr]
    l_frame_fronts = TravelingWaveSimulations.substantial_fronts.(l_frames, Ref(x))
    persistent_activation(l_frames, t, min_activation) && (mean(population(l_frames[end],1)) > 0.2 || is_traveling(l_frame_fronts, t, min_motion)) 
end

function is_traveling_solitary_radial_slice(exec, min_motion=1e-2, max_activation=5e-2)
    # We'll call it traveling solitary if:
    #  1. There is a traveling front
    #  2. No elevated activity trails the front
    l_frames = exec.solution.u
    t = timepoints(exec)
    x = [x[1] for x in space(exec).arr]
    l_frame_fronts = TravelingWaveSimulations.substantial_fronts.(l_frames, Ref(x))
    !persistent_activation(l_frames, t, max_activation) && is_traveling(l_frame_fronts, t, min_motion)
end

# %%
ABS_STOP=300.0
line_example = TravelingWaveSimulations.@EI_kw_example function example(N_ARR=1,N_CDT=1,P=2; SNR_scale=80.0, stop_time=ABS_STOP,
                                                     Aee=70.0, See=25.0,
                                                     Aii=2.0, Sii=27.0,
                                                     Aie=35.0, Sie=25.0,
                                                     Aei=70.0, Sei=27.0,
                                                     n=256, x=700.0, 
                                                     stim_strength=6.0,
                                                     stim_width=28.1,
                                                     stim_duration=7.0,
                                                     save_idxs_arg=[IndexSubsampler((5,)), Simulation73.RightCutFrom((0.0,))])
  simulation = Simulation(
    WCMSpatial(;
      pop_names = ("E", "I"),
      α = (1.0, 1.0),
      β = (1.0, 1.0),
      τ = (3.0, 3.0),
      nonlinearity = pops(DifferenceOfSigmoids;
        firing_θ = [6.0, 11.4],
        firing_a = [1.2, 1.0],
        blocking_θ = [30.0, 30.0],
        blocking_a = [1.2, 1.0]),
      stimulus = pops(SharpBumpStimulusParameter;
          strength = [stim_strength, stim_strength],
          width = [stim_width, stim_width],
          time_windows = [[(0.0, stim_duration)], [(0.0, stim_duration)]],
          baseline=[0.0, 0.0]),
      connectivity = FFTParameter(pops(GaussianConnectivityParameter;
          amplitude = [Aee -Aei;
                       Aie -Aii],
          spread = [(See,) (Sei,);
                    (Sie,) (Sii,)]
         ))
      );
      space = PeriodicLattice{Float64,N_ARR}(; n_points=(n,), extent=(x,)),
      tspan = (0.0,stop_time),
      dt = 0.1,
      algorithm=Tsit5(),
      save_idxs = save_idxs_arg,
      callback=DiscreteCallback(if !(save_idxs_arg === nothing)
        (u,t,integrator) -> begin
                    sub_u = u[integrator.opts.save_idxs];
                    t > 5 && ((all(isapprox.(sub_u, 0.0, atol=1e-4)) || (sub_u[end] > 0.005)))
                end
    else
        (u,t,integrator) -> begin
                    pop = population(u,1)
                    t > 5 && ((all(isapprox.(u, 0.0, atol=1e-4)) || (pop[end] > 0.005)))
            end
    end, terminate!)
  )
end



# %%
line_exec = execute(line_example())

# %%
@show is_traveling_solitary_radial_slice(line_exec)
line_anim = custom_animate(line_exec)
mp4(noisy_anim, "tmp/line_anim.mp4")

# %%
