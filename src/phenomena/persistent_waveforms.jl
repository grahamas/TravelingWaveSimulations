############################
### Persistent Waveform ###
############################

struct Persistent{T,AT<:AbstractVector{T}}#{WAVE<:TravelingWaveSimulations.AbstractWaveform,T,ARR_WAVE<:AbstractArray{WAVE,1},ARR_TIME<:AbstractArray{T,1}}
    waveforms
    t::AT
end
approx_sign(x, eps=1e-5) = abs(x) > eps ? sign(x) : 0
get_stop_place(p::Persistent) = slope_loc(p.waveforms[end])
get_start_time(p::Persistent) = p.t[1]
get_stop_time(p::Persistent) = p.t[end]
get_slope_sign(p::Persistent) = approx_sign(slope_val(p.waveforms[1]))
Base.length(p::Persistent) = length(p.waveforms)
function get_velocities(p::Persistent)
    # Discard the last frame because it's a repeat
    dxs = diff([slope_loc(wf) for wf in p.waveforms])[1:end-1]
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
function Base.push!(persistent::Persistent, (wf, t)::Tuple{WAVE,T}) where {WAVE,T}
    push!(persistent.waveforms, wf) 
    push!(persistent.t, t)
end

function waveform_identity_distance(front1::WF, front2::WF) where {WF <: TravelingWaveSimulations.Wavefront}
    if sign(slope_val(front1)) != sign(slope_val(front2))
        return Inf
    end
    abs(slope_loc(front1) - slope_loc(front2))
end

# Constructs list of all persistent fronts within array of arrays of fronts
# How to handle when new front appears near old front?
function link_persistent_fronts(frame_fronts::AbstractArray{<:AbstractArray}, ts::AbstractArray{T}, max_vel=200) where {T<:Number}#, WF<:TravelingWaveSimulations.Wavefront{T}}
    PARRTYPE = Persistent
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
                slope_dists[idx] = Inf
                continue
            end
            push!(matched_fronts, front_idx)
            push!(matched_travelers, traveler_idx)
            push!(active_travelers[traveler_idx], (fronts[front_idx], t))
            slope_dists[idx] = Inf
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
