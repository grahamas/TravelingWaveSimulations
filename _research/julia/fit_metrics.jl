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
using TravelingWaveSimulations: Value, linear_interpolate
const MaybeData{T} = Union{T,Missing}

TravelingWaveSimulations.Value(loc::Int, val::T) where T = (@assert loc > 0; Value{Int,T}(loc,val))
Base.:(==)(a::Value, b::Value) = a.val == b.val
Base.:(==)(a::Value, b::Number) = a.val == b
Base.isless(a::Value, b::Value) = a.val < b.val
Base.isless(a::Value, b::Number) = a.val < b
Base.isless(a::Number, b::Value) = a < b.val

# function Base.minimum(arr::Array{<:Value})
#     val, idx = findmin([elt.val for elt in arr]) 
#     return arr[idx]
# end


function linear_interpolate_loc((left,right)::Tuple{Value,Value}, val)
    loc = linear_interpolate((left.val,right.val), (left.loc,right.loc), val)
    Value(loc, val)
end
function linear_interpolate_val((left,right)::Tuple{Value,Value}, loc)
    val = linear_interpolate((left.loc,right.loc), (left.val,right.val), loc)
    Value(loc, val)
end

all_but_last(itr, n=1) = Iterators.take(itr, length(itr)-n)

struct ValuedSpace{T,C,N,AT<:AbstractArray{T,N},AC<:AbstractArray{C,N}} <: AbstractArray{T,N}
    values::AT
    coordinates::AC
end

function ValuedSpace(values::AT, coordinates::AC) where {T,C,N,AT<:AbstractArray{T,N},AC<:AbstractArray{C,N}}
    @assert size(values) == size(coordinates)
    return ValuedSpace{T,C,N,AT,AC}(values, coordinates)
end

Base.getindex(vs::ValuedSpace, I::AbstractArray) = [vs[i] for i in I]
Base.getindex(vs::ValuedSpace, idx::Union{CartesianIndex,Int}) = Value(vs.coordinates[idx], vs.values[idx])
function Base.getindex(vs::ValuedSpace, fidx::AbstractFloat)
    lidx = findlast(vs.coordinates .< fidx)
    if lidx == length(vs)
        return vs[lidx]
    end
    val = linear_interpolate_val((vs[lidx], vs[lidx+1]), fidx)
    return val
end
Base.size(vs::ValuedSpace) = size(vs.values)
function zero_crossing(A::AbstractArray, i, j, atol=1e-4)
    if (A[i] <= 0 && A[j] >= 0) || (A[i] >= 0 && A[j] <= 0)
        return linear_interpolate_loc((A[i],A[j]), 0.0)
    elseif (A[i] < -atol && A[j] >= -atol) || (A[i] > atol && A[j] <= atol)
        return A[j]
    else
        return nothing
    end
end

Base.diff(vs::ValuedSpace) = ValuedSpace(diff(vs.values) ./ diff(vs.coordinates), 0.5 .* diff(vs.coordinates) .+ collect(all_but_last(vs.coordinates)))

function translate(val::Value, vs::ValuedSpace)
    vs[val.loc]
end

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
# RIGHT WAVEFRONT FIXME
struct Wavefront{T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}} <: AbstractWaveform{T_LOC,T_VAL}
    left::V
    slope::V
    right::V
end
function translate(wf::Wavefront, args...)
    Wavefront(
        translate(wf.left, args...),
        translate(wf.slope, args...),
        translate(wf.right, args...)
        )
end

function score(wf::Wavefront, width_θ::T=5) where T
    apex = max(wf.left, wf.right)
    left_score = drop_score(apex, wf.left, width_θ)
    right_score = drop_score(apex, wf.right, width_θ)
    score = max(right_score, left_score) # TODO account for slope
    return valid_score(score)
end

struct WaveformMetrics{T}
    left_height::T
    right_height::T
    slope::T
    slope_loc::T
    width::T
end
function metrics(wf::Wavefront)
    WaveformMetrics(
        wf.left.val,
        wf.right.val,
        wf.slope.val,
        wf.slope.loc,
        wf.right.loc - wf.left.loc        
    ) 
end
    
choose_best(::Nothing, s::Scored) = s
choose_best(s1::Sd, s2::Sd) where {OBJ,Sd <: Scored{OBJ}} = (s1.score >= s2.score) ? s1 : s2
function choose_best(arr1::AbstractArray{<:Sd}, arr2::AbstractArray{<:Sd}) where {OBJ,Sd <: Scored{OBJ}}
    sum((elt) -> elt.score, arr1) > sum((elt) -> elt.score, arr2) ? arr1 : arr2
end


function find_min_between((left,right)::Tuple{Value,Value}, vs::ValuedSpace)
    xs = vs.coordinates
    maybe_left = findlast(xs .< left.loc)
    left_idx = if maybe_left === nothing
        findfirst(xs .>= left.loc)
    else
        maybe_left + 1
    end
    maybe_right = findfirst(xs .> right.loc)
    right_idx = if maybe_right === nothing
        findlast(xs .<= right.loc)
    else
        maybe_right - 1
    end
    return minimum(vs[left_idx:right_idx])
end

# Divide between peaks of second derivative (inflection points)
# Center on peaks of first derivative

function detect_all_fronts(valued_space::ValuedSpace)
    d_values = diff(valued_space)
    dd_values = diff(valued_space)
    fronts = Wavefront{Float64,Float64}[]
    left_boundary = valued_space[CartesianIndex(1)]
    d_steepest_slope = nothing
    for idx=all_but_last(eachindex(d_values))
        this_d_steepest_slope = zero_crossing(dd_values, idx, idx+1)
        if this_d_steepest_slope !== nothing
            if d_steepest_slope === nothing
                d_steepest_slope = this_d_steepest_slope
            elseif this_d_steepest_slope.val > d_steepest_slope.val
                d_steepest_slope = this_d_steepest_slope
            end
        end
        d_right_boundary = zero_crossing(d_values, idx, idx+1, 1e-4)
        if d_right_boundary !== nothing
            steepest_slope = if d_steepest_slope === nothing
                find_min_between((left_boundary, d_right_boundary), d_values)
            else
                translate(d_steepest_slope, d_values)
            end
            right_boundary = translate(d_right_boundary, valued_space)
            push!(fronts, Wavefront(left_boundary,
                                    steepest_slope,
                                    right_boundary)
            )
            d_steepest_slope = nothing
            left_boundary = right_boundary
            d_right_boundary = nothing
        end
    end
    right_boundary = valued_space[end]
    steepest_slope = if d_steepest_slope === nothing
        find_min_between((left_boundary, right_boundary), d_values)
    else
        translate(d_steepest_slope, d_values)
    end
    push!(fronts, Wavefront(left_boundary,
                            steepest_slope,
                            right_boundary)
    )
    return fronts
end

# function consolidate_fronts(fronts::AbstractVector{<:Wavefront{Float64}}, frame, slope_min, xs)
#     if length(fronts) == 0
#         return fronts
#     end
#     new_fronts = Wavefront[]
#     first_suff_dx = 1
#     while abs(fronts[first_suff_dx].front_slope) < slope_min
#         if first_suff_dx < length(fronts)
#             first_suff_dx += 1
#         else
#             break
#         end
#     end
#     current_front = if fronts[first_suff_dx].left.loc > 1
#         max_val, max_dx = findmax(frame[1:fronts[first_suff_dx].right.loc])
#         apex = value_from_idx(max_dx, xs, frame)
#         Wavefront(value_from_idx(1, xs, frame), fronts[first_suff_dx].front, fronts[first_suff_dx].front_slope, apex, fronts[first_suff_dx].right)
#     else
#         fronts[first_suff_dx]
#     end
#     if length(fronts) > first_suff_dx
#         for i in first_suff_dx+1:length(fronts)
#             next_front = fronts[i]
#             if next_front.front_slope < slope_min
#                 idxs = space_index((current_front.left.loc,next_front.right.loc), xs)
#                 max_val, max_dx = findmax(frame[idxs])
#                 apex = Value(xs[idxs[max_dx]], max_val)
#                 current_front = Wavefront(current_front.left, current_front.front, current_front.front_slope, apex, next_front.right)
#             else
#                 push!(new_fronts, current_front)
#                 current_front = next_front
#             end
#         end
#     end
#     if current_front.right.loc < length(frame)
#         idxs = space_index((left_bd.loc, rightmost.loc), xs)
#         (max_val, max_idx) = findmax(signal[idxs])
#         apex = Value(xs[idxs[max_idx]], signal[max_idx])
#         current_front = Wavefront(current_front.left, current_front.front, current_front.front_slope, apex, Value(length(frame), frame[end]))
#     end
#     push!(new_fronts, current_front)
#     return new_fronts
# end

function substantial_fronts(frame, xs::AbstractVector, slope_min=1e-4)
    vs = ValuedSpace(frame, xs)
    all_fronts = detect_all_fronts(vs)
    @show all_fronts
    return all_fronts
#     consolidated = consolidate_fronts(all_fronts, frame, slope_min, xs)
#     @show consolidated
#     return consolidated
end

function scored_rightmost_wavefront(frame::AbstractVector{T}, xs::AbstractVector)::Scored{<:Wavefront{T,T},T} where T
    fronts = substantial_fronts(frame, xs)
    if length(fronts) == 0
        return fronts
    else
        return Scored(fronts[end])
    end
end
scored_rightmost_wavefront(multipop::AbstractArray{T,2}, xs) where T = scored_rightmost_wavefront(multipop[:,1], xs)

function metrics_df(scored_wave_arr::AbstractArray{<:Scored{WAVE}}) where {WAVE <: AbstractWaveform}
    waveform_metrics = [metrics(s.obj) for s in scored_wave_arr]
    scores = [s.score for s in scored_wave_arr]
    wave_metric_syms = fieldnames(WaveformMetrics)
    wave_metrics_df = DataFrame(Dict(zip(wave_metric_syms, [[getproperty(wave, sym) for wave in waveform_metrics] for sym in wave_metric_syms])))
end

function traveling_wave_metrics_linear_regression(wave_metrics_df::DataFrame, t)
    wave_metrics_syms = names(wave_metrics_df)
    mostly_extant_cols = filter(c -> count(ismissing, wave_metrics_df[:,c])/size(wave_metrics_df,1) < 0.01, wave_metrics_syms)
    wave_metrics_df.t = t
    results = map(mostly_extant_cols) do metric_sym
        fmla = Term(metric_sym) ~ Term(:t) + ConstantTerm(1)
        glm_fit = glm(fmla, wave_metrics_df, Normal(), IdentityLink(); wts=scores)
        return metric_sym => glm_fit
    end |> Dict
    return (results, scores, wave_metrics_df)
end

function tw_metrics(frames, t, x)
    traveling_wave_metrics_linear_regression(metrics_df(scored_rightmost_wavefront.(frames, Ref(x))), t)
end

function tw_metrics(exec::Execution)
    u = exec.solution.u
    t = timepoints(exec)
    x = [x[1] for x in space(exec).arr]
    tw_metrics(u, t, x)
end

# %%
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

# %%
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

# %%
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

@recipe function f(vs::ValuedSpace)
    (vs.coordinates, vs.values)
end

@recipe function f(wf_arr::Array{<:Wavefront})
    seriestype := :scatter
    @series begin
        color := :green
        marker := :star
        markersize := 10
        [wf.left.loc for wf in wf_arr], [wf.left.val for wf in wf_arr]
    end
    @series begin
        color := :red
        [wf.right.loc for wf in wf_arr], [wf.right.val for wf in wf_arr]
    end
    @series begin
        color := :blue
        [wf.slope.loc for wf in wf_arr], [wf.slope.val for wf in wf_arr]
    end
end
vs = ValuedSpace(complex_wave,x)
wfs = detect_all_fronts(vs)
plot(vs)
plot!(wfs) |> display
plot(diff(vs))
plot!(wfs, ylim=[-0.05,0.05]) |> display
# plot!(diff(diff(vs)))
# # plot!(float_index.(Ref(x), (3:length(x)-1) .- 0.5), diff(complex_wave)|> diff |> diff) |> display
# plot(vs)
# cwfs = consolidate_fronts(detect_all_fronts(complex_wave,x), complex_wave, 1e-4, x)
# @show x
# plot!(cwfs) |> display
# plot(x[2:end], diff(complex_wave)) |> display 
# plot(x[3:end], diff(diff(complex_wave))) |> display

# %%
diff(vs)[162.0]
