export tw_metrics,
    SolitaryWaveformMetrics, WavefrontMetrics,
    ValuedSpace

const MaybeData{T} = Union{T,Missing}
abstract type AbstractWaveform{T_LOC,T_VAL} end
abstract type AbstractWaveformMetrics{T} end

struct Scored{OBJ,SCR}
    obj::OBJ
    score::SCR
end
Scored(obj::OBJ) where OBJ = Scored(obj, score(obj))

sigmoid(x::Real) = one(x) / (one(x) + exp(-x))
drop_proportion(apex::Value, nadir::Value) = min(abs(apex.val - nadir.val) / abs(apex.val), 1.0)
drop_duration(apex::Value, nadir::Value) = abs(apex.loc - nadir.loc)
drop_score(apex::Value{T}, nadir::Value{T}, θ) where T = drop_proportion(apex, nadir) * sigmoid((drop_duration(apex, nadir) - θ))
valid_score(score) = if isnan(score)
        return 0.0
    elseif score < 0
        # Confirm score range
        @warn "negative score"
        return 0.0
    else
        return score
    end 
# RIGHT WAVEFRONT FIXME -- It doesn't look like it's right-side?
struct Wavefront{T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}} <: AbstractWaveform{T_LOC,T_VAL}
    left::V
    slope::V
    right::V
end
slope_loc(wv::Wavefront) = wv.slope.loc
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
struct WavefrontMetrics{T} <: AbstractWaveformMetrics{T}
    left_height::T
    right_height::T
    slope::T
    slope_loc::T
    width::T
end
function metrics(wf::Wavefront)
    WavefrontMetrics(
        wf.left.val,
        wf.right.val,
        wf.slope.val,
        wf.slope.loc,
        wf.right.loc - wf.left.loc        
    ) 
end




# Divide between peaks of second derivative (inflection points)
# Center on peaks of first derivative

function detect_all_fronts(valued_space::ValuedSpace)
    "Partition space at extrema"
    d_values = diff(valued_space)
    dd_values = diff(d_values)
    fronts = Wavefront{Float64,Float64,Value{Float64,Float64}}[]
    left_boundary = getvalue(valued_space, 1)
    steepest_slope = nothing
    for idx=collect(eachindex(d_values))[2:end-2]
        right_boundary = translate(zero_crossing(d_values, idx, idx+1, 1e-4), valued_space)
        if right_boundary !== nothing
            this_front = getslice(d_values, (left_boundary, right_boundary))
            _, midx = findmax(abs.(this_front))
            steepest_slope = getvalue(this_front, midx)
            push!(fronts, Wavefront(left_boundary,
                                    steepest_slope,
                                    right_boundary)
            )
            steepest_slope = nothing
            left_boundary = right_boundary
            right_boundary = nothing
        end
    end
    right_boundary = getvalue(valued_space, length(valued_space))
    this_front = getslice(d_values, (left_boundary, right_boundary))
    _, midx = findmax(abs.(this_front))
    steepest_slope = getvalue(this_front, midx)
    push!(fronts, Wavefront(left_boundary,
                            steepest_slope,
                            right_boundary)
    )
    return fronts
end

eat_left(left, right) = Wavefront(left.left, right.slope, right.right)
eat_right(left, right) = Wavefront(left.left, left.slope, right.right)
function consolidate_fronts(fronts::AbstractVector{WF}, vs::ValuedSpace, slope_min=1e-4) where {WF <: Wavefront{Float64}}
    "Insist that all fronts have a minimum slope, otherwise consolidate"
    if length(fronts) == 0
        return fronts
    end
    new_fronts = WF[]
    first_suff_dx = 1
    while abs(fronts[first_suff_dx].slope.val) < slope_min
        if first_suff_dx < length(fronts)
            first_suff_dx += 1
        else
            return WF[]
        end
    end
    
    current_front = eat_left(fronts[1], fronts[first_suff_dx])
    if length(fronts) > first_suff_dx
        for i in first_suff_dx+1:length(fronts)
            next_front = fronts[i]
            if abs(next_front.slope.val) < slope_min
                current_front = eat_right(current_front, next_front)
            else
                push!(new_fronts, current_front)
                current_front = next_front
            end
        end
    end
    push!(new_fronts, current_front)
    return new_fronts
end

# TODO deal with multipop
function substantial_fronts(multipop::AbstractMatrix, xs::AbstractVector, slope_min=1e-4)
    substantial_fronts(population(multipop,1), xs, slope_min)
end
function substantial_fronts(frame::AbstractVector, xs::AbstractVector, slope_min=1e-4)
    vs = ValuedSpace(frame, xs)
    all_fronts = detect_all_fronts(vs)
    consolidated = consolidate_fronts(all_fronts, vs, slope_min)
    return consolidated
end
