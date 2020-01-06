using StaticArrays

export StaticWaveStats, TravelingWaveStats

# Analysis functions
struct Value{T_LOC,T_VAL}
    loc::T_LOC
    val::T_VAL
end
from_idx_to_space(val::Value, space) = Value(space[val.loc], val.val)
struct LinearFit{T}
    coefs::SVector{2,T}
end
LinearFit(coefs::AbstractArray{T}) where T = LinearFit{T}(SVector{2,T}(coefs))
linear_coef(lf::LinearFit) = lf.coefs[1]
const_coef(lf::LinearFit) = lf.coefs[2]
(lf::LinearFit)(a) = make_matrix(LinearFit, a) * lf.coefs
make_matrix(::Type{LinearFit}, a) = [a ones(size(a,1))]
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
function Base.print(io::IO, fe::FitErr)
    print(io, "$(fe.fit); error: $(fe.err)")
end
struct TravelingWaveStats{T_LOC,T_VAL}
    width::FitErr{T_LOC}
    center::FitErr{T_LOC}
    amplitude::FitErr{T_VAL}
    score::T_VAL
end
valid_metric(metric, score, min_score=0.0001) = !ismissing(metric) && (!ismissing(score) && (score > min_score))
function TravelingWaveStats(stats_arr::AbstractArray{<:StaticWaveStats{T_LOC,T_VAL}}, t) where {T_LOC,T_VAL}
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
    width_linfit = linreg_dropmissing(widths, t, scores)
    center_linfit = linreg_dropmissing(centers, t, scores)
    amplitude_linfit = linreg_dropmissing(amplitudes, t, scores)
    TravelingWaveStats(FitErr(widths, width_linfit, t, scores), 
        FitErr(centers, center_linfit, t, scores),
        FitErr(amplitudes, amplitude_linfit, t, scores),
        norm(scores))
end
# TODO: don't rely on 1D traveling wave
function TravelingWaveStats(exec::Execution)
    u = exec.solution.u
    t = exec.solution.t
    x = [x[1] for x in exec.solution.x]
    
    static_stats = StaticWaveStats.(population.(u,1), Ref(x))
    TravelingWaveStats(static_stats, t)
end
function Base.show(io::IO, ::MIME"text/plain", tws::TravelingWaveStats)
    println(io, "$(typeof(tws)):")
    for fname in fieldnames(TravelingWaveStats)
        println(io, "    $fname = $(getfield(tws,fname))")
    end
end
velocity(tws::TravelingWaveStats) = linear_coef(tws.center.fit)

function linreg_dropmissing(b_with_missing, A_with_missing, weights)
    # Ax = b
    A_with_missing = make_matrix(LinearFit, A_with_missing)
    notmissing = valid_metric.(b_with_missing, weights)
    A = A_with_missing[notmissing,:]
    b = b_with_missing[notmissing]
    W = diagm(weights[notmissing])
    x = nothing
    try
        x = (A' * W * A) \ (A' * W * b)
    catch e
        @show e
        @show b_with_missing
        @show A_with_missing
        @show A
    end
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
    elseif score < 0
        @warn "negative score: left: $left, apex: $apex, right: $right"
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