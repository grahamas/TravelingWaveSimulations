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
using Simulation73, TravelingWaveSimulations, Plots, Optim, LinearAlgebra, Distances, Statistics,
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
# Define an activation as a local D1 peak
abstract type AbstractActivation{T_LOC,T_SCR} <: AbstractScored{T_LOC} end


struct RightmostActivation{T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}} <: AbstractActivation{T_LOC,T_VAL}
    left::V
    max_derivative::V
    right::V
    score::T_VAL
    function RightmostActivation{T_LOC,T_VAL,V}(l::V,m::V,r::V) where {T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}}
        new{T_LOC,T_VAL,V}(l,m,r,score(RightmostActivation,l,m,r))
    end
end
RightmostActivation(l::V,m::V,r::V) where {T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}} = RightmostActivation{T_LOC,T_VAL,V}(l,m,r)
const RightmostActivationIdx = RightmostActivation{Int}
from_idx_to_space(idx::RightmostActivationIdx, space) = RightmostActivation(from_idx_to_space(idx.location,space), idx.score)

function score(::Type{<:RightmostActivation},l::V,m::V,r::V,θ=5.0) where {T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}}
    left_score = drop_score(l, m, θ)
    right_score = drop_score(m, r, θ)
    score = max(left_score, right_score)
    if isnan(score)
        return 0
    elseif score < 0
        @warn "negative score: left: $left, apex: $apex, right: $right"
        return 0
    else
        return score
    end
end

# Ughhhh needs to be robust to noise... (incl numerical noise)
function best_activation(frame::AbstractArray{T,1}, min_wing_dxs::Int=5) where {T<:Number}
    d_frame = frame[2:end] - frame[1:end-1]
    left = if d_frame[1] < 0
        Value(1, d_frame[1])
    else
        nothing
    end
    potential_d_mins = if d_frame[1] < d_frame[2]
        [left]
    else
        []
    end
    right = nothing
    for idx in 2:(length(frame)-1)
        if left === nothing
            if d_frame[idx] <= 0
                left = Value(idx, d_frame[idx])
            end
        end
        if d_frame[idx] < 0
            if d_frame[idx-1] >= d_frame[idx] < d_frame
                push!(potential_d_mins, Value(idx, d_frame[idx]))
            end
        else
            right = Value(idx, d_frame[idx])
            if d_frame[idx] >= 0
                right = Value(idx, d_frame[idx])
                this_activation = reduce(potential_d_mins) do d_min1, d_min2
                    choose_best(RightActivation(left,dmin1,right),
                                RightActivation(left,dmin2,right))
                end
                best_activation = choose_best(best_activation, this_activation)
                right = left = nothing
                potential_d_mins = []
            end
        end
    end
    if left !== nothing
        if d_frame[end] < d_frame[end-1]
            push!(potential_d_mins, Value(length(frame), d_frame[end]))
        end
        final_activation = reduce(potential_d_mins) do d_min1, d_min2
            choose_best(RightActivation(left,dmin1,right),
                        RightActivation(left,dmin2,right))
        end
        best_activation = choose_best(best_activation, final_activation)
    end
    return best_activation
end

function RightmostActivationStats(ra::RightmostActivation{T_LOC,T_VAL}, ra_val, ra_space)
    
end

function RightmostActivationStats(frame::AbstractArray{T,1}, space::AbstractArray{T,1})
    ra_idx_obj = best_activation(frame)
    ra_val = frame[ra_idx_obj.left.loc:ra_idx_obj.right.loc] 
    ra_space = space[ra_idx_obj.left.loc:ra_idx_obj.right.loc]
    ra_obj = from_idx_to_space(ra_idx_obj, space)
    RightmostActivationStats(ra_obj, ra_val, ra_space)
end

# %%

choose_best(::Nothing, s::AbstractScored) = s
choose_best(s1::AbstractScored, s2::AbstractScored) = (s1.score >= s2.score) ? s1 : s2
function choose_best(arr1::AbstractArray{<:Sd}, arr2::AbstractArray{<:Sd}) where {Sd <: AbstractScored}
    sum((elt) -> elt.score, arr1) > sum((elt) -> elt.score, arr2) ? arr1 : arr2
end

# Ughhhh needs to be robust to noise... (incl numerical noise)
function detect_peak_with_wings(frame::AbstractArray{T,1}, min_wing_dxs::Int=5,
        peak_type::Type{<:AbstractPeak}=StaticPeak) where {T<:Number}
    # Find peak, with surrounding left and right d_frame zeros
    apex = right = max_scoring_peak = nothing
    left = Value(1, frame[1])
    prev_val = frame[1]
    idx = 2
    for val in frame[2:end]
        change = val - prev_val
        if apex === nothing
            if change <= 0
                apex = Value(idx-1, prev_val)
            end
        else
            if change >= 0
                right = Value(idx-1, prev_val)
                peak = peak_type(left, apex, right)
                max_scoring_peak = choose_max_scoring(max_scoring_peak, peak)
                apex = nothing
                left = right
            end
        end
        prev_val = val
        idx += 1
    end
    rightmost = Value(idx-1, prev_val)
    apex = apex === nothing ? rightmost : apex
    final_peak = peak_type(left, apex, rightmost)
    max_peak = choose_peakiest(max_peak, final_peak)
    return max_peak
end

function WaveSolitary(frame::AbstractArray{T,1}, min_wing_dxs::Int=10) where T
    detect_peak_with_wings(frame, min_wing_dxs)
end
function WaveFront(frame::AbstractArray{T,1}, min_wing_dxs::Int=5) where T
    detect_peak_with_wings(diff(frame))
end

function TravelingWave(frames)
    choose_best(WaveSolitary.(frames), WaveFront.(frames))
end

function stats(wave_arr::AbstractArray{<:WaveSolitary})
    widths = [st.width for st in wave_arr]
    centers = [st.center for st in wave_arr]
    amplitudes = [st.amplitude for st in wave_arr]
    scores = [st.score for st in wave_arr]
    
    if sum(valid_metric.(widths, scores)) < 5 || sum(valid_metric.(centers, scores)) < 5 || sum(valid_metric.(amplitudes,scores)) < 5 #At least five wave frames
        return nothing
    end

    if mean(scores) < 1e-2 || mean(amplitudes) < 1e-3 # roughly, less than 1% of the run contains TW
        return nothing # can't try fitting; singular
    end
    width_linfit = linreg_dropmissing(widths, t, scores)
    center_linfit = linreg_dropmissing(centers, t, scores)
    amplitude_linfit = linreg_dropmissing(amplitudes, t, scores)
    StatsWaveSolitary(FitErr(widths, width_linfit, t, scores), 
        FitErr(centers, center_linfit, t, scores),
        FitErr(amplitudes, amplitude_linfit, t, scores),
        norm(scores))
end

function stats(wave_arr::AbstractArray{<:WaveFront}, t)
    front_locs = [st.front for st in wave_arr]
    amplitudes = [st.amplitude for st in wave_arr]
    left_baselines = [st.left_baseline for st in wave_arr]
    right_baselines = [st.right_baseline for st in wave_arr]
    scores = [st.score for st in wave_arr]
    
    if sum(valid_metric.(right_baselines, scores)) < 5 || sum(valid_metric.(front_locs, scores)) < 5 || sum(valid_metric.(amplitudes, scores)) < 5 || sum(valid_metric.(left_baselines,scores)) < 5 #At least five wave frames
        return nothing
    end
    
    if mean(scores) < 1e-2 || mean(amplitudes) < 1e-3 # roughly, less than 1% of the run contains TW
        return nothing # can't try fitting; singular
    end
    front_locs_lf = linreg_dropmissing(front_locs, t, scores)
    amplitudes_lf = linreg_dropmissing(amplitudes, t, scores)
    left_baselines_lf = linreg_dropmissing(left_baselines, t, scores)
    right_baselines_lf = linreg_dropmissing(right_baselines, t, scores)
    StatsWaveFront(
        FitErr(front_locs, t, scores),
        FitErr(amplitudes, t, scores),
        FitErr(left_baselines, t, scores),
        FitErr(right_baselines, t, scores)
    )
end

function stats(frames)
    stats(TravelingWave(frames))
end
    
    


