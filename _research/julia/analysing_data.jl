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
#     display_name: Julia 1.3.0
#     language: julia
#     name: julia-1.3
# ---

# %% [markdown]
# TODO: Needs parsing of filename for non-pkey params

# %%
using Simulation73, TravelingWaveSimulations, Plots, Optim, LinearAlgebra, Distances
using DiffEqBase: AbstractTimeseriesSolution

# %%
using Lazy, JuliaDB

data_root = joinpath(homedir(), "sim_data")

function _load_data(sim_path)
    @warn "loading"
   JuliaDB.load(sim_path)
end
"Load most recent simulation"
function load_data(data_root, example_name)
    nsims = length(readdir(joinpath(data_root, example_name)))
    load_data(data_root, example_name, nsims)
end
"Load nth simulation, ordered by time"
function load_data(data_root, example_name, nth::Int)
    ex_path = joinpath(data_root, example_name)
    sims = readdir(ex_path)
    sorted_sims = sort(joinpath.(Ref(ex_path), sims), by=mtime)
    sim_path = sorted_sims[nth]
    return (get_example(example_name), @lazy map(_load_data, joinpath.(Ref(sim_path), readdir(sim_path))))
end

macro ifsomething(ex)
    quote
        result = $(esc(ex))
        result === nothing && return nothing
        result
    end
end
struct MultiDBRowIter
    dbs
end
function _get_row(dbs, db, dbs_state, row_state...)
    row_tuple = iterate(db, row_state...)
    while row_tuple === nothing
        db, dbs_state = @ifsomething iterate(dbs)
        row_tuple = iterate(db)
    end
    (row, row_state) = row_tuple
    return (row, (dbs_state, db, row_state))
end
function iterate(it::MultiDBRowIter)
    (db, dbs_state) = @ifsomething iterate(it.dbs)
    return _get_row(it.dbs, db, dbs_state)
end
function iterate(it::MultiDBRowIter, (db_state, db, row_state))
    return _get_row(it.dbs, db, dbs_state, row_state)
end
    


# %%
(example, dbs) = load_data(data_root, "sigmoid_normal_fft");

# %%
n_db = 5
keys1 = select(dbs[n_db], Keys())
vals1 = select(dbs[n_db], Not(Keys()))

# %%
n=5
mdl1 = example(; keys1[n]...)
example_exec = Execution(mdl1, BareSolution(; vals1[n]...));

# %% collapsed=true jupyter={"outputs_hidden": true}
anim1 = custom_animate(example_exec)
mp4(anim1, "tmp/tmp.mp4")

# %% jupyter={"source_hidden": true}
# working code

function solitary_peak(frame::AbstractArray{T,1}, n_dxs_per_regime=10) where {T<:Number}
    prev_val = frame[1]
    left_min = frame[1]
    uphill_count = 0
    have_peaked = false
    downhill_count = 0 # downhill can only come after long enough uphill
    cur_peak = nothing
    max_solitary_peakiness = 0
    for val in frame[2:end]
        if !have_peaked # still looking...
            if val - prev_val > 0 #continuing uphill
                uphill_count += 1
                downhill_count = 0
            else
                if uphill_count >= n_dxs_per_regime
                    have_peaked = true
                    cur_peak = prev_val
                    downhill_count = 1
                else
                    uphill_count = 0
                    left_min = val
                end
            end
        else
            if val - prev_val < 0 # continuing downhill
                downhill_count += 1
            else
                if downhill_count >= n_dxs_per_regime
                    right_drop_proportion = (cur_peak - prev_val) / cur_peak
                    left_drop_proportion = (cur_peak - left_min) / cur_peak
                    max_solitary_peakiness = max(max_solitary_peakiness, right_drop_proportion * left_drop_proportion)
                end
                left_min = prev_val
                have_peaked = false
                downhill_count = 0
                uphill_count = 1
            end
        end
        prev_val = val
    end
    if downhill_count >= n_dxs_per_regime
        right_drop_proportion = (cur_peak - prev_val) / cur_peak
        left_drop_proportion = (cur_peak - left_min) / cur_peak
        max_solitary_peakiness = max(max_solitary_peakiness, right_drop_proportion * left_drop_proportion)
    end
    return max_solitary_peakiness
end

# %%
struct Value{T_LOC,T_VAL}
    loc::T_LOC
    val::T_VAL
end
from_idx_to_space(val::Value, space) = Value(space[val.loc], val.val)
struct StaticPeak{T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}}
    left::V
    apex::V
    right::V
    score::T_VAL
end
StaticPeak(left::V, apex::V, right::V) where {T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}} = StaticPeak{T_LOC,T_VAL,V}(left,apex,right,peakiness(left,apex,right))
StaticPeak(left::V, apex::V, right::V,score) where {T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}} = StaticPeak{T_LOC,T_VAL,V}(left,apex,right,score)
const StaticPeakIdx = StaticPeak{Int}
from_idx_to_space(sp::StaticPeakIdx, space) = StaticPeak(from_idx_to_space(sp.left,space), from_idx_to_space(sp.apex,space), from_idx_to_space(sp.right,space), sp.score)
struct TravelingPeak{T_LOC,T_VAL,V<:Value{T_LOC,T_VAL},A_TIME<:AbstractArray{<:T_VAL}, A_PEAK<:AbstractArray{<:StaticPeak{<:V}}}
    times::A_TIME
    peaks::A_PEAK
end

function linear_interpolate((x1,x2), (y1,y2), x)
    @assert x1 <= x <= x2 || x2 <= x <= x1
    y1 + (x - x1) * ((y2 - y1) / (x2 - x1)) 
end
struct StaticWaveStats{T_LOC,T_VAL}
    baseline::T_VAL
    amplitude::T_VAL
    width::Union{T_LOC,Missing}
    center::T_LOC
    score::Float64
end
function StaticWaveStats(sp::StaticPeak{T_LOC,T_VAL}, wave_val, wave_space) where {T_LOC,T_VAL}
    baseline = min(sp.left.val, sp.right.val)
    amplitude = sp.apex.val - baseline
    half_max = amplitude / 2
    width = if wave_val[1] >= half_max || wave_val[end] >= half_max
        missing
    else
        first = findfirst(wave_val .>= half_max) # must be greater than 1
        last = findlast(wave_val .>= half_max) # must be less than end
        left = linear_interpolate(wave_val[[first-1,first]], wave_space[[first-1,first]], half_max)
        right = linear_interpolate(wave_val[[last,last+1]], wave_space[[last,last+1]], half_max)
        right - left
    end
    return StaticWaveStats{T_LOC,T_VAL}(baseline, amplitude, width, sp.apex.loc, sp.score)
end
function StaticWaveStats(frame::AbstractArray{T,1}, space::AbstractArray{T,1}) where T
    peak_idx_obj = peakiest_peak(frame)
    peak_val = frame[peak_idx_obj.left.loc:peak_idx_obj.right.loc]
    peak_space = space[peak_idx_obj.left.loc:peak_idx_obj.right.loc]
    peak_obj = from_idx_to_space(peak_idx_obj, space)
    StaticWaveStats(peak_obj, peak_val, peak_space)
end

function calc_err(data, fit, x, w)
    notmissing = .!ismissing.(data)
    weuclidean(data[notmissing], fit(x[notmissing]), w[notmissing]) / sum(notmissing)
end
struct FitErr{T}
    fit::LinearFit
    err::T
end
FitErr(y, lf::LinearFit, x, w::AbstractArray{T}) where T = FitErr{T}(lf, calc_err(y, lf, x, w))
struct TravelingWaveStats{T_LOC,T_VAL}
    width::FitErr{T_LOC}
    center::FitErr{T_LOC}
    amplitude::FitErr{T_VAL}
end
function TravelingWaveStats(stats_arr::AbstractArray{<:StaticWaveStats{T_LOC,T_VAL}}, t) where {T_LOC,T_VAL}
    widths = [st.width for st in stats_arr]
    centers = [st.center for st in stats_arr]
    amplitudes = [st.amplitude for st in stats_arr]

    scores = [st.score for st in stats_arr]
    width_linfit = linreg_dropmissing(widths, t, scores)
    center_linfit = linreg_dropmissing(centers, t, scores)
    amplitude_linfit = linreg_dropmissing(amplitudes, t, scores)
    TravelingWaveStats(FitErr(widths, width_linfit, t, scores), 
        FitErr(centers, center_linfit, t, scores),
        FitErr(amplitudes, amplitude_linfit, t, scores))
end
# TODO: don't rely on 1D traveling wave
function TravelingWaveStats(exec::Execution)
    u = exec.solution.u
    t = exec.solution.t
    x = [x[1] for x in exec.solution.x]
    
    static_stats = static_wave_stats.(population.(u,1), Ref(x))
    TravelingWaveStats(static_stats, t)
end

struct LinearFit{X}
    x::X
end
(lr::LinearFit)(a) = make_matrix(LinearFit, a) * lr.x
make_matrix(::Type{LinearFit}, a) = [a ones(size(a,1))]
function linreg_dropmissing(b_with_missing, A_with_missing, weights)
    A_with_missing = make_matrix(LinearFit, A_with_missing)
    notmissing = .!ismissing.(b_with_missing) .& (b_with_missing .> 0.0001)
    A = A_with_missing[notmissing,:]
    b = b_with_missing[notmissing]
    W = diagm(weights[notmissing])
    x = (A' * W * A) \ (A' * W * b)
    return LinearFit(x)
end
    

# Peakiness
drop_proportion(apex::Value, nadir::Value) = ((apex.val - nadir.val) / apex.val)
drop_duration(apex::Value, nadir::Value) = abs(apex.loc - nadir.loc)
sigmoid(x::Real) = one(x) / (one(x) + exp(-x))
drop_score(apex::Value{T}, nadir::Value{T}, θ::T) where T = drop_proportion(apex, nadir) * sigmoid((drop_duration(apex, nadir) - θ))
function peakiness(left::Value{T}, apex::Value{T}, right::Value{T}, duration_θ::T=10) where T
    left_score = drop_score(apex, left, duration_θ)
    right_score = drop_score(apex, right, duration_θ)
    score = left_score * right_score
    if isnan(score)
        return 0
    else
        return score
    end
end
choose_peakiest(::Nothing, peak::StaticPeak) = peak
choose_peakiest(peak1::StaticPeak, peak2::StaticPeak) = (peak1.score >= peak2.score) ? peak1 : peak2
    
function peakiest_peak(frame::AbstractArray{T,1}, n_dxs_per_regime=10) where {T<:Number}
    left = nothing
    apex = nothing
    right = nothing
    peakiest = nothing
    max_peak = nothing
    
    prev_val = frame[1]
    left = Value(1, frame[1])
    idx = 2
    for val in frame[2:end]
        change = val - prev_val
        if apex === nothing
            if change <= 0 # = so that plateaus ruin peak score.
                apex = Value(idx-1,prev_val)
            end
        else
            if change >= 0
                right = Value(idx-1,prev_val)
                peak = StaticPeak(left, apex, right)
                max_peak = choose_peakiest(max_peak, peak)
                apex = nothing
                left = right
            end
        end
        prev_val = val
        idx +=1
    end
    rightmost = Value(idx-1,prev_val)
    apex = apex === nothing ? rightmost : apex
    final_peak = StaticPeak(left,apex,rightmost)
    max_peak = choose_peakiest(max_peak, final_peak)
    return max_peak
end

# %% jupyter={"source_hidden": true}
using LsqFit, LinearAlgebra
function sech2(x::AbstractArray{T}, p) where {T<:NTuple{1}}
    x1s = [y[1] for y in x]
    @. p[1] * sech(p[2] * (x1s + p[3]))^p[4]
end
function traveling_sech2(x::AbstractArray{T},t::Number,p) where {T<:NTuple{1}}
    x1s = [y[1] for y in x]
    @. p[1] * sech(p[2] * (x1s + p[3] - p[4] * t)) ^ p[5]
end
function Dx(x)
    x[2:end] - x[1:end-1]
end
sigmoid(x::Real) = one(x) / (one(x) + exp(-x))

function fit_single_sech2_optim(frame, coords)
    #xt = [coords..., t]
    p_guess = [0.5, 0.01, -75.0,2.0]
    p_lb = [0.0, 0.0, -150.0,0.5]
    p_ub = [1.0, 1.0, 150.0,10.0]
    cost = let y = frame, x = coords
        function (p)
            norm(y .- sech2(x,p))
        end
    end
    f = OnceDifferentiable(cost, zeros(length(p_guess)), autodiff=:forward)
    Optim.optimize(f, p_lb, p_ub, p_guess, Fminbox(LBFGS()))
    #curve_fit(sech2, coords, frame, p_guess)#, lower=p_lb, upper=p_ub)
end

function fit_traveling_sech2(frames, coords, ts) where {T}
    p_guess = [0.5, 0.01, -75.0, 1.0, 2.0]
    p_lb = [0.0, 0.0, -150.0,0.0,0.5]
    p_ub = [1.0, 1.0, 150.0,100.0,10.0]
    cost = let samples = [(frames[idx], coords, ts[idx]) for idx in 1:length(ts)]
        function (p)
            norm(
                map(samples) do (frame, coord, t)
                    frame .- traveling_sech2(coord, t, p)
                end
            )
        end
    end
    f = OnceDifferentiable(cost, zeros(length(p_guess)), autodiff=:forward)
    result = Optim.optimize(f, p_lb, p_ub, p_guess, Fminbox(LBFGS()))
end
                    

# TODO: Don't ignore other populations
function fit_traveling_wave_subset(exec, n_dxs_per_regime::Int=10, is_solitary_peak_threshold::Float64=0.99)
    e_pop = population.(exec.solution.u, 1)
    solitary_peakiness = solitary_peak.(e_pop, n_dxs_per_regime)
    wave_dxs = solitary_peakiness .> is_solitary_peak_threshold
    wave_frames = e_pop[wave_dxs]
    p_optim = fit_traveling_sech2(wave_frames, coordinates(exec), exec.solution.t[wave_dxs])    
end
        
# using Optim
# function fit_single_sech2(frame, coords, loss=(x,y) -> norm(x .- y))
#     #xt = [coords..., t]
#     p_guess = [1.0, 1.0, 0.0]
#     p_lb = [0.0, -Inf, -Inf]
#     p_ub = [Inf, Inf, Inf]
#     optimize((p) -> loss(frame, sech2(coords, p)), p_guess, BFGS())#, lower=p_lb, upper=p_ub)
# end


# %% jupyter={"source_hidden": true}
example_frame = population(example_exec.solution.u[10], 1)
example_coords = [x[1] for x in example_exec.solution.x]
plot(space(example_exec), example_frame) |> display
static_wave_stats(example_frame, example_coords)
# plot!(space(example_exec), [0.0, (Dx(example_frame) .< 0)...])
#plot!(space(example_exec), [0.0, (Dx(Dx(example_frame)) .< 0)...,0.0])
# plot!(space(example_exec), [0.0, (Dx(Dx(Dx(example_frame))) .< 0)..., 0.0, 0.0])
#plot!(space(example_exec), [0.0, 0.0, (Dx(example_frame) .< 0)..., 0.0, 0.0])
#@show solitary_peak(example_frame)

# %% collapsed=true jupyter={"outputs_hidden": true, "source_hidden": true}
example_frame = population(example_exec.solution.u[35], 1)
@show example_frame[end-30]
coord_tuples = example_exec.solution.x
fit = fit_single_sech2_optim(example_frame, coord_tuples)
@show fit
@show Optim.minimizer(fit)
plot(space(example_exec), example_frame) |> display
plot!(space(example_exec), sech2(coord_tuples, Optim.minimizer(fit))) |> display
@show solitary_peak(example_frame, 10)

# %%
# fits = map(example_exec.solution.u) do u
#     fit_single_sech2_optim(population(u,1), example_coords)
# end
# plot([fit.param[1] for fit in fits]) |> display
# plot([fit.param[2] for fit in fits]) |> display
# plot([fit.param[3] for fit in fits]) |> display
# plot([norm(fit.resid) for fit in fits]) |> display
# plot([log(norm(fit.resid) / abs(fit.param[1])) for fit in fits]) |> display
stats_arr = StaticWaveStats.(population.(example_exec.solution.u,1), Ref(example_coords))
t = example_exec.solution.t
widths = [st.width for st in stats_arr]
amplitudes = [st.amplitude for st in stats_arr]
centers = [st.center for st in stats_arr]

scores = [st.score for st in stats_arr]
width_linfit = linreg_dropmissing(widths, t, scores)
amplitude_linfit = linreg_dropmissing(amplitudes, t, scores)
center_linfit = linreg_dropmissing(centers, t, scores)

# plot([scores, [score .> 0.0001 ? score : missing for score in scores]]) |> display
# plot([widths, width_linfit(t)]) |> display
# plot([amplitudes, amplitude_linfit(t)]) |> display
# plot([centers, center_linfit(t)]) |> display

traveling_stats = TravelingWaveStats(example_exec)


# %%
result = fit_traveling_wave_subset(example_exec)
p_optim = Optim.minimizer(result)

# %%
e_pop = [population(u, 1) for u in example_exec.solution.u]
anim = animate(e_pop[solitary_peak.(e_pop) .> 0.99])# traveling_sech2(example_coords, example_exec.solution.t[30], p_optim))

# %%
function custom_animate_with_fit(execution::Execution{T,<:Simulation{T}}; kwargs...) where T
    solution = execution.solution
    pop_names = execution.simulation.model.pop_names
    x = space(execution)
    t = timepoints(execution)
    max_val = maximum(solution)
	min_val = minimum(solution)
    fitted_result = fit_traveling_wave_subset(execution)
    fitted_wave = traveling_sech2.(Ref(coordinates(x)), t, Ref(Optim.minimizer(fitted_result)))
    @animate for time_dx in 1:length(t) # TODO @views
        plot(
            x, population_timepoint(solution, 1, time_dx); 
            val_lim=(min_val,max_val), title="t = $(round(t[time_dx], digits=4))",
            xlab = "Space (a.u. approx. um)",kwargs...
            )
        plot!(x,  fitted_wave[time_dx])
    end
end

anim = custom_animate_with_fit(example_exec) 
mp4(anim, "tmp/tmp.mp4")
