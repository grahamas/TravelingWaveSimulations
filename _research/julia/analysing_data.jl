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

# %% [markdown]
# TODO: Needs parsing of filename for non-pkey params

# %%
using Simulation73, TravelingWaveSimulations, Plots, Optim, LinearAlgebra, Distances, Statistics,
    IterTools
using DiffEqBase: AbstractTimeseriesSolution

# %%
valid_metric(metric, score, min_score=0.0001) = !ismissing(metric) && (!ismissing(score) && (score > min_score))
function TravelingWaveSimulations.TravelingWaveStats(stats_arr::AbstractArray{<:StaticWaveStats{T_LOC,T_VAL}}, t) where {T_LOC,T_VAL}
    widths = [st.width for st in stats_arr]
    centers = [st.center for st in stats_arr]
    amplitudes = [st.amplitude for st in stats_arr]
    scores = [st.score for st in stats_arr]
    
    if sum(valid_metric.(widths, scores)) < 5 || sum(valid_metric.(centers, scores)) < 5 || sum(valid_metric.(amplitudes,scores)) < 5 #At least five wave frames
        return nothing
    end

    if mean(scores) < 1e-2 || mean(amplitudes) < 1e-3 # roughly, less than 1% of the run contains TW
        return nothing # can't try fitting; singular
    end
    width_linfit = TravelingWaveSimulations.linreg_dropmissing(widths, t, scores)
    center_linfit = TravelingWaveSimulations.linreg_dropmissing(centers, t, scores)
    amplitude_linfit = TravelingWaveSimulations.linreg_dropmissing(amplitudes, t, scores)
    TravelingWaveStats(TravelingWaveSimulations.FitErr(widths, width_linfit, t, scores), 
        TravelingWaveSimulations.FitErr(centers, center_linfit, t, scores),
        TravelingWaveSimulations.FitErr(amplitudes, amplitude_linfit, t, scores),
        norm(scores))
end
equals_str(key,val) = "$key=$val"
equals_strs(mods) = [equals_str(p...) for p in pairs(mods)]
mods_filename(x) = join(equals_strs(x), "_")

# %%
function linreg_dropmissing(b_with_missing, A_with_missing, weights)
    # Ax = b
    A_with_missing = TravelingWaveSimulations.make_matrix(LinearFit, A_with_missing)
    notmissing = valid_metric.(b_with_missing, weights)
    A = A_with_missing[notmissing,:]
    b = b_with_missing[notmissing]
    W = diagm(weights[notmissing])
    x = nothing
    try
        x = (A' * W * A) \ (A' * W * b)
    catch e
        @show b_with_missing
        @show A_with_missing
        @show A
    end
    return LinearFit(x)
end

# %%
# Load most recent simulation data
data_root = joinpath(homedir(), "sim_data")
(example, mdb) = TravelingWaveSimulations.load_data(data_root, "sigmoid_normal_fft");
example_name = TravelingWaveSimulations.get_example_name(mdb.fns[1])
sim_name = TravelingWaveSimulations.get_sim_name(mdb.fns[1])

# %%
# Analyse and extract twscore
# mdb_execs = MultiDBExecIter(example, dbs, ())
GC.gc()
mods = TravelingWaveSimulations.get_mods(mdb)
mod_names = keys(mods) |> collect
mod_values = values(mods) |> collect
A_is_traveling = Array{Bool}(undef, length.(values(mods))...)
A_velocity= Array{Union{Float64,Missing}}(undef, length.(values(mods))...)
A_velocity_errors= Array{Union{Float64,Missing}}(undef, length.(values(mods))...)
for db in mdb
    for (this_mod, exec) in DBExecIter(example, db, ())
        this_mod_key = keys(this_mod)
        this_mod_val = values(this_mod)
        A_idx = TravelingWaveSimulations.mod_idx(this_mod_key, this_mod_val, mod_names, mod_values)
        tws = TravelingWaveStats(exec);
        if tws === nothing
            A_is_traveling[A_idx] = false
            A_velocity[A_idx] = missing
            A_velocity_errors[A_idx] = missing
        else
            A_is_traveling[A_idx] = true
            A_velocity[A_idx] = TravelingWaveSimulations.velocity(tws)
            A_velocity_errors[A_idx] = tws.center.err
        end
    end
end
@show sum(ismissing.(A_velocity))
@show prod(size(A_velocity))

# %%
function mean_skip_missing(A; dims)
    missings = ismissing.(A)
    zeroed = copy(A)
    zeroed[missings] .= 0
    nonmissingsum = sum(zeroed; dims=dims)
    nonmissingmean = nonmissingsum ./ sum(.!missings; dims=dims)
    return nonmissingmean
end
all_dims = 1:length(mod_names)
for (x,y) in IterTools.subsets(all_dims, Val{2}())
    collapsed_dims = Tuple(setdiff(all_dims, (x,y)))
    is_traveling = dropdims(mean_skip_missing(A_is_traveling, dims=collapsed_dims), dims=collapsed_dims)
    velocities = dropdims(mean_skip_missing(A_velocity, dims=collapsed_dims), dims=collapsed_dims)
    velocity_errors = dropdims(mean_skip_missing(A_velocity_errors, dims=collapsed_dims), dims=collapsed_dims)
    prop_notmissing = dropdims(mean(.!ismissing.(A_velocity), dims=collapsed_dims), dims=collapsed_dims)
    plot(
        #heatmap(mod_values[x], mod_values[y], is_traveling, xlab=mod_names[x], ylab=mod_names[y], title="\"peakiness\" avgd across other spreads"),
        heatmap(mod_values[y], mod_values[x], velocities, xlab=mod_names[y], ylab=mod_names[x], title="velocity avgd"),
        heatmap(mod_values[y], mod_values[x], velocity_errors, xlab=mod_names[y], ylab=mod_names[x], title="error"),
        heatmap(mod_values[y], mod_values[x], prop_notmissing, xlab=mod_names[y], ylab=mod_names[x], title="prop not missing")
        ) |> display
    path = "tmp/$(example_name)/$(sim_name)/$(mod_names[x])_$(mod_names[y])_centerfiterror.png"
    mkpath(dirname(path))
    png(path)
end

# %%
mods = Dict(:See=>20.0, :Sei=>30.0)
anim = custom_animate(execute(example(;mods...)))
mp4(anim, "tmp/$(example_name)/$(sim_name)/anim_$(mods_filename(mods)).mp4")

# %% collapsed=true jupyter={"outputs_hidden": true}
function satisfies_criteria(exec, mod_keys, mod_vals, (param_criteria, tws_criterion))
    satisfies_params = all( param_criterion(mod_vals[findfirst(mod_keys .== param)])
        for (param, param_criterion) in param_criteria
    )
    if satisfies_params
        tws = TravelingWaveStats(exec)
        if tws_criterion(tws)
            return true
        end
    end
    return false
end
    

criteria = ([(:Sei, (x) -> x >= 34),(:See, (x) -> x <= 16)], (tws) -> velocity(tws) >= 15)
for (this_mod, exec) in MultiDBExecIter(example, mdb, ())
    this_mod_key = keys(this_mod)
    this_mod_val = values(this_mod)
    if satisfies_criteria(exec, this_mod_key, this_mod_val, criteria)
        tws = TravelingWaveStats(exec)
        params = join(["$name=$val" for (name,val) in zip(this_mod_key, this_mod_val)], "_")
        anim = custom_animate(exec; title="$(params); vel=$(velocity(tws))")
        mp4(anim, "tmp/$(example_name)/$(sim_name)/$(params).mp4")
        break
    end
end

# %% jupyter={"outputs_hidden": true} collapsed=true
# SECOND TIME
macro summarysizeall()
    allnames = names(Main)
    allsizes = [:(Base.summarysize($name)) for name in allnames]
    quote
        @show $(allsizes...)
    end
end
@summarysizeall

# %% collapsed=true jupyter={"outputs_hidden": true}
anim1 = custom_animate(example_exec)
mp4(anim1, "tmp/tmp.mp4")

# %% jupyter={"outputs_hidden": true} collapsed=true
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


# %%
example_frame = population(example_exec.solution.u[10], 1)
example_coords = [x[1] for x in example_exec.solution.x]
plot(space(example_exec), example_frame) |> display
StaticWaveStats(example_frame, example_coords)
# plot!(space(example_exec), [0.0, (Dx(example_frame) .< 0)...])
#plot!(space(example_exec), [0.0, (Dx(Dx(example_frame)) .< 0)...,0.0])
# plot!(space(example_exec), [0.0, (Dx(Dx(Dx(example_frame))) .< 0)..., 0.0, 0.0])
#plot!(space(example_exec), [0.0, 0.0, (Dx(example_frame) .< 0)..., 0.0, 0.0])
#@show solitary_peak(example_frame)

# %% collapsed=true jupyter={"outputs_hidden": true}
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

plot([scores, [score .> 0.0001 ? score : missing for score in scores]], title="\"peakiness\" score") |> display
plot([widths, width_linfit(t)], title="width") |> display
plot([amplitudes, amplitude_linfit(t)], title="amplitude") |> display
plot([centers, center_linfit(t)], title="center") |> display

traveling_stats = TravelingWaveStats(example_exec)
@show traveling_stats

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
