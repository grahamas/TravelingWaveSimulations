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
    IterTools, Combinatorics, DataFrames, GLM

# %% jupyter={"source_hidden": true}
multi_factor_name(name_combo) = join(name_combo, "_")

function calculate_factor_matrix(mdb, max_order)
    mods = TravelingWaveSimulations.get_mods(mdb)
    mod_names = keys(mods) |> collect # These must be same order v
    mod_values = values(mods) |> collect # These must be same order ^
    dx_combos = cat((Combinatorics.with_replacement_combinations.(Ref(1:length(mod_names)), 1:max_order) .|> collect)..., dims=1)
    name_combos = map((dx_combo) -> mod_names[dx_combo], dx_combos)
    col_combos = map((dx_combo) -> mod_names[dx_combo], dx_combos)
    
    factors = Array{Float64}(undef, length.(mod_values)..., length(dx_combos))
    for factor_dx in CartesianIndices(factors)
        dxs = Tuple(factor_dx)
        all_mod_dxs = dxs[1:end-1]
        combo_dx = dxs[end]
        these_mod_vals = mod_values[dx_combos[combo_dx]]
        these_mod_dxs = all_mod_dxs[dx_combos[combo_dx]]
        mod_vals = getindex.(these_mod_vals, these_mod_dxs)
        factors[factor_dx] .= prod(mod_vals)
    end
    
    factor_names = multi_factor_name.(name_combos)
    return (factor_names, factors)
end



# %%
data_root = joinpath(homedir(), "sim_data")
(example, mdb) = TravelingWaveSimulations.load_data(data_root, "sigmoid_normal_fft",4);
example_name = TravelingWaveSimulations.get_example_name(mdb.fns[1])
sim_name = TravelingWaveSimulations.get_sim_name(mdb.fns[1])

# %%
# Analyse and extract twscore
# mdb_execs = MultiDBExecIter(example, dbs, ())
GC.gc()
mods = TravelingWaveSimulations.get_mods(mdb)
mod_names = keys(mods) |> collect
mod_values = values(mods) |> collect
A_tws = Array{Float64}(undef, length.(values(mods))...)
A_velocity= Array{Union{Float64,Missing}}(undef, length.(values(mods))...)
A_velocity_errors= Array{Union{Float64,Missing}}(undef, length.(values(mods))...)
for db in mdb
    for (this_mod, exec) in DBExecIter(example, db, ())
        this_mod_key = keys(this_mod)
        this_mod_val = values(this_mod)
        A_idx = TravelingWaveSimulations.mod_idx(this_mod_key, this_mod_val, mod_names, mod_values)
        tws = TravelingWaveStats(exec);
        if tws === nothing
            A_tws[A_idx] = 0.0
            A_velocity[A_idx] = missing
            A_velocity_errors[A_idx] = missing
        else
            A_tws[A_idx] = tws.score
            A_velocity[A_idx] = TravelingWaveSimulations.velocity(tws)
            A_velocity_errors[A_idx] = tws.center.err
        end
    end
end
@show sum(ismissing.(A_velocity))
@show prod(size(A_velocity))

# %%
factor_names, factors = calculate_factor_matrix(mdb, 2);
parameters_mx = reshape(factors, (:,size(factors)[end]));
df = DataFrame(Dict(zip(factor_names, [parameters_mx[:,i] for i in 1:size(parameters_mx,2)])))
df.vel = A_velocity[:]
df.tws = A_tws[:]
wts = dropmissing(df).tws

# %%
bics = map(1:3) do i
    factor_names, factors = calculate_factor_matrix(mdb, i)
    parameters_mx = reshape(factors, (:, size(factors)[end]))
    df = DataFrame(Dict(zip(factor_names, [parameters_mx[:,i] for i in 1:size(parameters_mx, 2)])))
    df.vel = A_velocity[:]
    df.tws = A_tws[:]
    wts = dropmissing(df).tws
    fmla = Term(:vel) ~ sum(Term.(Symbol.(factor_names))) + ConstantTerm(1);
    lm_fit = lm(fmla, dropmissing(df))
    glm_fit = glm(fmla, dropmissing(df), Normal(), IdentityLink(); wts=wts)
    return StatsModels.bic.([lm_fit, glm_fit])
end
lm_bics, glm_bics = zip(bics...)
plot([lm_bics...], title="LM BIC") |> display
plot([glm_bics...], title="GLM BIC") |> display


# %%

# %% jupyter={"outputs_hidden": true} collapsed=true
# tests

test_xs = 0.0:0.1:20.0
test_sigmoid = 1.0 .- NeuralModels.simple_sigmoid_fn.(test_xs, 4.0, 10.0)
test_sech2 = NeuralModels.sech2_fn.(test_xs, 0.7, 10.0)
test_double_sigmoid = NeuralModels.simple_sigmoid_fn.(test_xs, 4.0, 10.0) + NeuralModels.simple_sigmoid_fn.(test_xs, 4.0, 15.0)
# plot([plot(test_xs, test_double_sigmoid), plot(test_xs, test_sech2)]...) |> display
# plot(plot(test_xs[2:end], -diff(test_double_sigmoid)), plot(test_xs[2:end], -diff(test_sech2)))|> display
# plot(plot(test_xs[3:end], diff(-diff(test_double_sigmoid))), plot(test_xs[3:end], diff(-diff(test_sech2)))) |> display
# plot(plot(test_xs[4:end], diff(diff(-diff(test_sigmoid)))), plot(test_xs[4:end], diff(diff(-diff(test_sech2))))) |> display

scored_rightmost_wavefront(test_sech2, test_xs)

# %% jupyter={"outputs_hidden": true} collapsed=true
## Difference of Sigmoids example

ABS_STOP=300.0
using WilsonCowanModel
dos_example = TravelingWaveSimulations.@EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; SNR_scale=80.0, stop_time=ABS_STOP,
                                                     Aee=280.0, See=70.0,
                                                     Aii=1.4, Sii=70.0,
                                                     Aie=270.0, Sie=90.0,
                                                     Aei=-297.0, Sei=90.0,
                                                     n=71, x=500.0)
  simulation = Simulation(
    WilsonCowanModel.WCMSpatial(;
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

# %% jupyter={"outputs_hidden": true} collapsed=true
using DifferentialEquations
execution = execute(dos_example(; n=256, x=700.0, See=25.0, Sii=25.0, Sie=27.0, Sei=27.0,
                                Aee=250.0, Aei=75.0, Aie=50.0, Aii=10.0, strengthE=10.0, widthE=50.0,
                                algorithm=Tsit5()));
anim = TravelingWaveSimulations.custom_animate(execution)
mp4(anim, "tmp/dos_tmp.mp4")

# %%
results, scores, df = tw_metrics(execution)

# %%
plot(df.t, df.front_slope)

# %%
exec = execution    
u = exec.solution.u
t = timepoints(exec)
x = [x[1] for x in space(exec).arr]
complex_wave = u[30][:,1]

vs = ValuedSpace(complex_wave,x)
wfs = detect_all_fronts(vs)
# plot(vs)
# plot!(wfs, vs) |> display
# plot(diff(vs))
# plot!(wfs, ylim=[-0.03,0.03]) |> display
# plot!(diff(diff(vs)), ylim=[-0.003, 0.003])|> display
# plot(getslice(vs, (-Inf,35.0)), xlim=[0.0,35.0], legend=false, seriestype=:scatter, markersize=20)
# plot!(wfs, vs) |> display
# plot(diff(vs), xlim=[0.0,35.0], legend=false, seriestype=:scatter)
# plot!([wfs[1].right.loc], [diff(vs)[wfs[1].right.loc].val], seriestype=:scatter) |> display
# plot!(diff(diff(vs)))
# # plot!(float_index.(Ref(x), (3:length(x)-1) .- 0.5), diff(complex_wave)|> diff |> diff) |> display
plot(vs)
cwfs = consolidate_fronts(wfs, vs)
plot!(cwfs) |> display
plot(x[2:end], diff(complex_wave)) |> display 
plot(x[3:end], diff(diff(complex_wave))) |> display

# %%
Base.Broadcast.broadcasted(f, vs::ValuedSpace, x) = Broadcast.broadcasted(Broadcast.ArrayStyle{ValuedSpace}(), f, vs, x);

x = vs .+ 1.0;
typeof(x)

# %%
tw_metrics(WavefrontMetrics, exec)

# %%
# Aee: E->E amplitude
# See: E->E spread

# caS: cross-auto spread ratio
# caA: cross-auto amplitude ratio

# ioA: inhibitory output amplitude scale
# ioS: inhibitory output spread scale
# iiA: inhibitory input amplitude scale
# iiS: inhibitory input spread scale
ABS_STOP = 300.0
sigmoid_example = TravelingWaveSimulations.@EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; stop_time=ABS_STOP,
                                                     Aee=24.0, See=25.0,
                                                     Aii=4.0, Sii=27.0,
                                                     Aie=27.0, Sie=25.0,
                                                     Aei=18.2, Sei=27.0,
                                                     n=128, x=700.0, stim_strength=1.2,
                                                     stim_width=28.1)
  simulation = Simulation(
    WCMSpatial(;
      pop_names = ("E", "I"),
      α = (1.0, 1.0),
      β = (1.0, 1.0),
      τ = (3.0, 3.0),
      nonlinearity = pops(SigmoidNonlinearity;
        a = [1.2, 1.0],
        θ = [2.6, 8.0]),
      stimulus =  pops(SharpBumpStimulusParameter;
          strength = [stim_strength,stim_strength],
          width = [stim_width, stim_width],
          time_windows = [[(0.0, 5.0)], [(0.0, 5.0)]]),
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


using DifferentialEquations
execution = execute(sigmoid_example(; algorithm=Tsit5()));
anim = TravelingWaveSimulations.custom_animate(execution)
mp4(anim, "tmp/sigmoid_tmp.mp4")

# %%
(results, scores, df) = tw_metrics(SolitaryWaveformMetrics, execution);

# %%
using StatsPlots
@df df plot(:t, :right_slope_loc)

# %%
results[:right_slope_loc]
