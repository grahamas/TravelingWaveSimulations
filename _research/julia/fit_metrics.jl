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
using Simulation73, NeuralModels, TravelingWaveSimulations, Plots, Optim, LinearAlgebra, Distances, Statistics,
    IterTools, Combinatorics, DataFrames, GLM

# %%
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

# %% collapsed=true jupyter={"outputs_hidden": true}
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
lm_bics

# %%
fit = TravelingWaveSimulations.linreg_dropmissing(A_velocity[:], parameters_mx, wts)

# %%
fmla = Term(:vel) ~ sum(Term.(Symbol.(factor_names))) + ConstantTerm(1);

# %% collapsed=true jupyter={"outputs_hidden": true}
lm_fit = lm(fmla, dropmissing(df))

# %%
glm_fit = glm(fmla, dropmissing(df), Normal(), IdentityLink(); wts=wts)

# %%
deviance(lm_fit)

# %%
deviance(glm_fit)

# %%
Value(loc::Int, val::T) where T = (@assert loc > 0; Value{Int,T}(loc,val))

const MaybeData{T} = Union{T,Missing}
using TravelingWaveSimulations: Value, linear_interpolate

function translate(val::Value{<:Union{Int,CartesianIndex}}, new_arr::AbstractArray, offset=0)
    Value(val.loc + offset, new_arr[val.loc + offset])
end
from_idx_to_space(val::Value, space) = Value(space[val.loc], val.val)

struct Scored{OBJ,SCR}
    obj::OBJ
    score::SCR
end
Scored(obj::OBJ) where OBJ = Scored(obj, score(obj))
abstract type AbstractWaveform{T_LOC,T_VAL} end

drop_proportion(apex::Value, nadir::Value) = ((apex.val - nadir.val) / apex.val)
drop_duration(apex::Value, nadir::Value) = abs(apex.loc - nadir.loc)
sigmoid(x::Real) = one(x) / (one(x) + exp(-x))
drop_score(apex::Value{T}, nadir::Value{T}, θ) where T = drop_proportion(apex, nadir) * sigmoid((drop_duration(apex, nadir) - θ))
valid_score(score) = if isnan(score)
        return 0.0
    elseif score < 0
        # Confirm score range
        @warn "negative score: left: $left, apex: $apex, right: $right"
        return 0.0
    else
        return score
    end 
struct SolitaryWaveform{T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}} <: AbstractWaveform{T_LOC,T_VAL}
    left::V
    apex::V
    right::V
    width::MaybeData{T_VAL} # Never in index units
end
function score(swf::SolitaryWaveform, width_θ::T=10) where T
    left_score = drop_score(swf.apex, swf.left, width_θ)
    right_score = drop_score(swf.apex, swf.right, width_θ)
    score = left_score * right_score
    return valid_score(score)
end
# RIGHT WAVEFRONT FIXME
struct Wavefront{T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}} <: AbstractWaveform{T_LOC,T_VAL}
    left::V
    front::V
    front_slope::T_VAL # Actually T_VAL / T_LOC
    apex::V
    right::V
end
function score(wf::Wavefront, width_θ::T=5) where T
    left_score = drop_score(wf.apex, wf.left, width_θ)
    right_score = drop_score(wf.apex, wf.right, width_θ)
    score = right_score - left_score
    return valid_score(score)
end

struct WaveformMetrics{T}
    left_baseline::T
    right_baseline::T
    apex_height::T
    apex_loc::T
    width::MaybeData{T}
    front_slope::T
    front_loc::T
end
function metrics(swf::SolitaryWaveform)
    WaveformMetrics(
        swf.left.val,
        swf.right.val,
        swf.apex.val,
        swf.apex.loc,
        swf.width,
        missing,
        missing
    )
end
function metrics(wf::Wavefront)
    WaveformMetrics(
        wf.left.val,
        wf.right.val,
        wf.apex.val,
        wf.apex.loc,
        missing,
        wf.front_slope,
        wf.front.loc        
    ) 
end
    
choose_best(::Nothing, s::Scored) = s
choose_best(s1::Sd, s2::Sd) where {OBJ,Sd <: Scored{OBJ}} = (s1.score >= s2.score) ? s1 : s2
function choose_best(arr1::AbstractArray{<:Sd}, arr2::AbstractArray{<:Sd}) where {OBJ,Sd <: Scored{OBJ}}
    sum((elt) -> elt.score, arr1) > sum((elt) -> elt.score, arr2) ? arr1 : arr2
end

function wave_width(wave_left, wave_apex, wave_right, wave_val, wave_space)
    half_max = (wave_apex.val - min(wave_left.val, wave_right.val)) / 2
    width = if wave_val[1] >= half_max || wave_val[end] >= half_max
        missing
    else
        first = findfirst(wave_val .>= half_max) # must be greater than 1
        last = findlast(wave_val .>= half_max) # must be less than end
        left = linear_interpolate(wave_val[[first-1,first]], wave_space[[first-1,first]], half_max)
        right = linear_interpolate(wave_val[[last,last+1]], wave_space[[last,last+1]], half_max)
        right - left
    end
    return width
end

# Divide between peaks of second derivative (inflection points)
# Center on peaks of first derivative

# TODO: principled way not relying on arbitrary atol
function is_extremum(i, vec::Vector, atol=1e-16)
    (vec[i-2] < vec[i-1] <= vec[i] >= vec[i+1] > vec[i+2] || vec[i-2] > vec[i-1] >= vec[i] <= vec[i+1] < vec[i+2]) 
end

function detect_all_fronts(signal::Vector)
    @assert length(signal) > 2
    d_signal = diff(signal)
    dd_signal = diff(d_signal)
    fronts = Wavefront[]
    left_bd = Value(1, signal[1])
    center = nothing
    first = true
    # TODO: indices are fucked up
    for i=3:length(signal)-4
        if is_extremum(i, d_signal)
            center = Value(i, signal[i])
            @assert !is_extremum(i, dd_signal)
            first = false
        elseif center === nothing && first && is_extremum(i, dd_signal)
            left_bd = Value(i, signal[i])
            first = false
        elseif center !== nothing && is_extremum(i, dd_signal)
            here = Value(i, signal[i])
            (max_val, max_dx) = findmax(signal[left_bd.loc:i])
            apex = Value(max_dx, max_val)
            push!(fronts, Wavefront(left_bd, center, d_signal[center.loc], apex, here))
            center = nothing
            left_bd = here
        end
    end
    if center !== nothing
        (max_val, max_dx) = findmax(signal[left_bd.loc:end])
        apex = Value(max_dx, max_val)
        push!(fronts, Wavefront(left_bd, center, d_signal[center.loc], apex, Value(length(signal), signal[end])))
    end
    return fronts
end

# Steps:
# 1. Detect local maxima and minima in derivative
# 2. Categorize as maxima (left facing, LF) and minima (right facing, RF)
# 3. Return list
# 4. Simple approach: consider only movement of right-most RF. Call "solitary waveform" if preceded by LF

function best_waveform(::Type{<:SolitaryWaveform}, frame::AbstractArray{T,1}, xs::AbstractArray, min_wing_dxs::Int=10) where T
    best_peak = detect_peak_with_wings(frame, xs, min_wing_dxs)
    return best_peak
end
function best_waveform(::Type{<:Wavefront}, frame::AbstractArray{T,1}, xs::AbstractArray, min_wing_dxs::Int=5) where T
    scored_max_d_peak = detect_peak_with_wings(-diff(frame), xs, min_wing_dxs)
    max_d_peak = scored_max_d_peak.obj
    front_left = translate(max_d_peak.left, frame, 1)
    front = translate(max_d_peak.apex, frame, 1)
    front_right = translate(max_d_peak.right, frame, 1)
    apex_val, apex_locp1 = findmax(frame[front_left.loc:front_right.loc])
    apex_loc = apex_locp1 + front_left.loc - 1
    apex = Value(apex_loc, apex_val)
    max_wavefront = Wavefront(
            from_idx_to_space(front_left, xs),
            from_idx_to_space(front, xs),
            max_d_peak.apex.val,
            from_idx_to_space(apex, xs),
            from_idx_to_space(front_right, xs)
        ) |> Scored
    return max_wavefront
end

function traveling_wave_metrics_linear_regression(scored_wave_arr::AbstractArray{<:Scored{WAVE}}, t) where {WAVE <: AbstractWaveform}
    waveform_metrics = [metrics(s.obj) for s in scored_wave_arr]
    scores = [s.score for s in scored_wave_arr]
    wave_metric_syms = fieldnames(WAVE)
    wave_metrics_df = DataFrame(Dict(zip(metric_syms, [[getproperty(wave, sym) for wave in waves] for sym in metric_syms])))
    wave_metrics_df.t = t
    # TODO: don't attempt to do bad fit
    fomulae = [(Term(metric) ~ :t + ConstantTerm(1)) for metric in wave_metric_syms]
    results = map(wave_metric_syms) do metric_sym
        fmla = Term(metric_sym) ~ :t + ConstantTerm(1)
        glm_fit = glm(fmla, dropmissing(wave_metrics_df), Normal(), IdentityLink(); wts=scores)
        return metric_sym => glm_fit
    end |> Dict
    return results
end


# vvvvvvvvvvvv Maybe... vvvvvvvvvvvvv

function TravelingWave(frames, x)
    choose_best(best_waveform.(SolitaryWaveform, frames, Ref(x)), best_waveform.(Wavefront, frames, Ref(x)))
end

function stats(frames, t, x)
    stats(stats(TravelingWave(frames), x), t)
end

function stats(exec::Execution)
    u = exec.solution.u
    t = exec.solution.t
    x = [x[1] for x in exec.solution.x]
    stats(u, t, x)
end

# %%
# tests

test_xs = 0.0:0.1:20.0
test_sigmoid = 1.0 .- NeuralModels.simple_sigmoid_fn.(test_xs, 4.0, 10.0)
test_sech2 = NeuralModels.sech2_fn.(test_xs, 0.7, 10.0)
test_double_sigmoid = NeuralModels.simple_sigmoid_fn.(test_xs, 4.0, 10.0) - NeuralModels.simple_sigmoid_fn.(test_xs, 4.0, 15.0)
plot([plot(test_xs, test_double_sigmoid), plot(test_xs, test_sech2)]...) |> display
plot(plot(test_xs[2:end], -diff(test_double_sigmoid)), plot(test_xs[2:end], -diff(test_sech2)))|> display
plot(plot(test_xs[3:end], diff(-diff(test_double_sigmoid))), plot(test_xs[3:end], diff(-diff(test_sech2)))) |> display
plot(plot(test_xs[4:end], diff(diff(-diff(test_sigmoid)))), plot(test_xs[4:end], diff(diff(-diff(test_sech2))))) |> display

detect_all_fronts(test_double_sigmoid)

# %%
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
using DifferentialEquations
execution = execute(dos_example(; n=256, x=700.0, See=25.0, Sii=25.0, Sie=27.0, Sei=27.0,
                                Aee=250.0, Aei=75.0, Aie=50.0, Aii=10.0, strengthE=10.0, widthE=50.0,
                                algorithm=Tsit5()));
anim = TravelingWaveSimulations.custom_animate(execution)
mp4(anim, "tmp/dos_tmp.mp4")
