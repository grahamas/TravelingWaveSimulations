
abstract type AbstractWaveform{T_LOC,T_VAL} end

#################
### Wavefront ###

# RIGHT WAVEFRONT FIXME -- It doesn't look like it's only right-side?
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

# Center on peaks of first derivative


function substantial_fronts(exec::AbstractFullExecution)
    substantial_fronts.(exec.solution.u, Ref([x[1] for x in space(exec).arr])) 
end

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


############################
### Persistent Wavefront ###
############################

struct Persistent{WAVE<:TravelingWaveSimulations.AbstractWaveform,
        T,ARR_WAVE<:AbstractArray{WAVE,1},ARR_TIME<:AbstractArray{T,1}}
    waveforms::ARR_WAVE
    t::ARR_TIME
end
approx_sign(x, eps=1e-5) = abs(x) > eps ? sign(x) : 0
get_stop_place(p::Persistent) = p.waveforms[end].slope.loc
get_start_time(p::Persistent) = p.t[1]
get_stop_time(p::Persistent) = p.t[end]
get_slope_sign(p::Persistent) = approx_sign(p.waveforms[1].slope.val)
Base.length(p::Persistent) = length(p.waveforms)
function get_velocities(p::Persistent)
    # Discard the last frame because it's a repeat
    dxs = diff([wf.slope.loc for wf in p.waveforms])[1:end-1]
    dts = diff(p.t)[1:end-1]
    return dxs ./ dts
end
estimate_velocity(::Nothing) = (missing, missing)
function estimate_velocity(p::Persistent)
    vels = get_velocities(p)
    est = mean(vels)
    err = mean(abs.(vels .- est))
    return (est, err)
end
function Base.push!(persistent::Persistent{WAVE,T}, (wf, t)::Tuple{WAVE,T}) where {WAVE,T}
    push!(persistent.waveforms, wf) 
    push!(persistent.t, t)
end

function waveform_identity_distance(front1::WF, front2::WF) where {WF <: TravelingWaveSimulations.Wavefront}
    if sign(front1.slope.val) != sign(front2.slope.val)
        return Inf
    end
    abs(TravelingWaveSimulations.slope_loc(front1) - TravelingWaveSimulations.slope_loc(front2))
end


function has_traveling_front(frame_fronts::AbstractArray{<:AbstractArray{WF}}, ts::AbstractArray{<:Number}, max_vel=20, min_vel=1e-3, min_traveling_frames=5) where {WF <:TravelingWaveSimulations.Wavefront}
    p_fronts = persistent_fronts(frame_fronts, ts, max_vel)
    any(is_traveling.(p_fronts, min_vel, min_traveling_frames)) 
end

function is_traveling(persistent::Persistent{<:TravelingWaveSimulations.Wavefront}, min_vel, min_traveling_frames)
    if length(persistent.waveforms) > min_traveling_frames
        vels = get_velocities(persistent)
        num_traveling_frames = sum(vels .> min_vel)
        return num_traveling_frames .> min_traveling_frames
    else
        return false
    end
    return false
end

# Constructs list of all persistent fronts within array of arrays of fronts
# How to handle when new front appears near old front?
function persistent_fronts(frame_fronts::AbstractArray{<:AbstractArray{WF}}, ts::AbstractArray{T}, max_vel=200) where {T<:Number, WF<:TravelingWaveSimulations.Wavefront{T,T,Value{T,T}}}
    PARRTYPE = Persistent{WF,T,Array{WF,1},Array{T,1}}
    prev_fronts = frame_fronts[1]
    prev_t = ts[1]
    inactive_travelers = PARRTYPE[]
    active_travelers = PARRTYPE[PARRTYPE([front], [t]) for (front, t) in zip(prev_fronts, prev_t)]
    for (fronts, t) in zip(frame_fronts[2:end], ts[2:end])
        dt = t - prev_t
        first_possible_dx = 1
        slope_dists = [
            waveform_identity_distance(front, traveler.waveforms[end]) for (front, traveler) in Iterators.product(fronts, active_travelers)
        ]
        matched_fronts = BitSet()
        matched_travelers = BitSet()
        i_matches = 1
        while i_matches <= min(length(fronts), length(active_travelers))
            dist, idx = findmin(slope_dists)
            if (dist / dt) > max_vel
                #@info "I broke!"
                break
            end
            front_idx, traveler_idx = Tuple(idx)
            if (front_idx in matched_fronts) || (traveler_idx in matched_travelers)
                #@warn "Ambiguous traveling waves."
                slope_dists[idx] .= Inf
                continue
            end
            push!(matched_fronts, front_idx)
            push!(matched_travelers, traveler_idx)
            push!(active_travelers[traveler_idx], (fronts[front_idx], t))
            slope_dists[idx] .= Inf
            i_matches += 1
        end
        unmatched_fronts = setdiff(BitSet(1:length(fronts)), matched_fronts)
        unmatched_travelers = setdiff(BitSet(1:length(active_travelers)), matched_travelers)
        append!(inactive_travelers, active_travelers[unmatched_travelers |> collect])
        active_travelers::Array{PARRTYPE,1} = active_travelers[matched_travelers |> collect]
        new_active = PARRTYPE[Persistent([front], [t]) for front in fronts[unmatched_fronts |> collect]]
        append!(active_travelers, new_active)
    end
    return PARRTYPE[inactive_travelers..., active_travelers...]
end
