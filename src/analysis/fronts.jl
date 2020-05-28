
struct Persistent{WAVE<:TravelingWaveSimulations.AbstractWaveform,T,ARR_WAVE<:AbstractArray{WAVE,1},ARR_TIME<:AbstractArray{T,1}}
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