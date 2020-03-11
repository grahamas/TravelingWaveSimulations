struct Persistent{WAVE<:TravelingWaveSimulations.AbstractWaveform,T,ARR_WAVE<:AbstractArray{WAVE,1},ARR_TIME<:AbstractArray{T,1}}
    waveforms::ARR_WAVE
    t::ARR_TIME
end
approx_sign(x, eps=1e-5) = abs(x) > eps ? sign(x) : 0
get_stop_place(p::Persistent) = p.waveforms[end].slope.loc
get_stop_time(p::Persistent) = p.t[end]
get_slope_sign(p::Persistent) = approx_sign(p.waveforms[1].slope.val)
function get_velocities(p::Persistent)
    dxs = diff([wf.slope.loc for wf in p.waveforms])
    dts = diff(p.t)
    return dxs ./ dts
end
function Base.push!(persistent::Persistent{WAVE,T}, (wf, t)::Tuple{WAVE,T}) where {WAVE,T}
    push!(persistent.waveforms, wf) 
    push!(persistent.t, t)
end
function is_traveling(persistent::Persistent{<:TravelingWaveSimulations.Wavefront}, min_vel=1e-3, min_traveling_frames=5)
    if length(persistent.waveforms) > min_traveling_frames
        dx = diff(TravelingWaveSimulations.slope_loc.(persistent.waveforms))
        dt = diff(persistent.t)
        vel = dx ./ dt
        num_traveling_frames = sum(vel .> min_vel)
        return num_traveling_frames .> min_traveling_frames
    else
        return false
    end
end

function has_traveling_front(frame_fronts::AbstractArray{<:AbstractArray{<:TravelingWaveSimulations.Wavefront}}, ts::AbstractArray{<:Number}, max_vel=20, min_vel=1e-3, min_traveling_frames=5)
    p_fronts = persistent_fronts(frame_fronts, ts, max_vel)
    any(is_traveling.(p_fronts, min_vel, min_traveling_frames)) 
end

function waveform_identity_distance(front1::WF, front2::WF) where {WF <: TravelingWaveSimulations.Wavefront}
    if sign(front1.slope.val) != sign(front2.slope.val)
        return Inf
    end
    abs(TravelingWaveSimulations.slope_loc(front1) - TravelingWaveSimulations.slope_loc(front2))
end

# How to handle when new front appears near old front?
function persistent_fronts(frame_fronts::AbstractArray{<:AbstractArray{<:TravelingWaveSimulations.Wavefront}}, ts::AbstractArray{<:Number}, max_vel=20)
    prev_fronts = frame_fronts[1]
    prev_t = ts[1]
    inactive_travelers = []
    active_travelers = [Persistent([front], [t]) for (front, t) in zip(prev_fronts, prev_t)]
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
        active_travelers = active_travelers[matched_travelers |> collect]
        append!(active_travelers, [Persistent([front], [t]) for front in fronts[unmatched_fronts |> collect]])
    end
    return [inactive_travelers..., active_travelers...]
end

function fronts(exec::AbstractExecution, slope_min=1e-2)
    l_frames = exec.solution.u
    t = timepoints(exec)
    x = [x[1] for x in space(exec).arr]
    # A list of lists
    # For each frame, a list of fronts
    # Plateaus can be absorbed into fronts, but as long as the whole plateau is moving, then that's fine. 
    l_frames_fronts = TravelingWaveSimulations.substantial_fronts.(l_frames, Ref(x), slope_min)
    
end

function TravelingWaveSimulations.custom_animate(execution::Execution{T,<:Simulation{T}}, fronts::AbstractArray{<:AbstractArray{<:TravelingWaveSimulations.AbstractWaveform}}; kwargs...) where T
    solution = execution.solution
    pop_names = execution.simulation.model.pop_names
    x = space(execution)
    x1 = [xcoord[1] for xcoord in x.arr]
    t = timepoints(execution)
    max_val = maximum(solution)
	min_val = minimum(solution)
    i_pop = 1
    @animate for time_dx in 1:length(t) # TODO @views
        frame = population_timepoint(solution, i_pop, time_dx)
        vs = ValuedSpace(frame, x1)
        plot(
            x, frame; label=pop_names[i_pop],
            val_lim=(min_val,max_val), title="t = $(round(t[time_dx], digits=4))",
            xlab = "Space (a.u. approx. um)",kwargs...
            )
        wf_arr = fronts[time_dx]
        scatter!([wf.slope.loc for wf in wf_arr], [vs[wf.slope.loc].val for wf in wf_arr])
    end
end


function is_activating_front(pf::Persistent)
    # Either the vel and slope have opposing signs
    # or it's a static front
    vel_sgn = approx_sign(mean(get_velocities(pf)))
    slope_sgn = get_slope_sign(pf)
    if vel_sgn * slope_sgn <= 0
        return true
    else
        return false
    end
end

function will_be_overtaken(persistent_front::Persistent, contemporary_fronts::AbstractArray{<:Persistent})
    pf_velocity = mean(get_velocities(persistent_front))
    pf_velocity_sign = approx_sign(pf_velocity)
    pf_stop_place = get_stop_place(persistent_front)
    for front in contemporary_fronts
        front_velocity = mean(get_velocities(front))
        if ((abs(front_velocity) > abs(pf_velocity)) # front is moving faster
                && (sign(pf_stop_place - get_stop_place(front)) == approx_sign(front_velocity))) # front velocity dir matches relative displacement
            return true # the persistent_front will be deactivated by front overtaking
        end
    end
    return false
end

function will_be_cancelled(persistent_front::Persistent, contemporary_fronts::AbstractArray{<:Persistent})
    pf_velocity = mean(get_velocities(persistent_front))
    pf_velocity_sign = approx_sign(pf_velocity)
    pf_stop_place = get_stop_place(persistent_front)
    for front in contemporary_fronts
        front_velocity = mean(get_velocities(front))
        if (isapprox(abs(front_velocity), abs(pf_velocity), atol=1e-5) # front is moving faster
                && (sign(pf_stop_place - get_stop_place(front)) == approx_sign(front_velocity))) # front velocity dir matches relative displacement
            return true # the persistent_front will be deactivated by front overtaking
        end
    end
    return false
end

function has_solitary_traveling_wave(l_activating_fronts, l_final_fronts)
    any(map(filter(is_traveling, l_activating_fronts)) do traveling_front
             will_be_cancelled(traveling_front, l_final_fronts)
        end
    )
end

function is_decaying(persistent_front::Persistent)
    maxes = map(persistent_front.waveforms) do wf
        max(wf.left.val, wf.right.val)
    end
    num_decreasing = 0
    for delta in diff(maxes)[end:-1:1]
        if delta < 0
            num_decreasing += 1
        else
            break
        end
    end
    if num_decreasing > 5
        return true
    end
    return false
end

function will_be_deactivated(persistent_front::Persistent, contemporary_fronts::AbstractArray{<:Persistent})
    # Check if it's decaying
    if is_decaying(persistent_front)
        return true
    end
    
    # Find a deactivate-to-zero that will overtake it
    if will_be_overtaken(persistent_front, contemporary_fronts)
        return true
    end
    
    # Neither decaying nor will be overtaken
    return false
end

function persistent_activation(l_frames, t, min_activation)
    # Test if the left-most value is ever high and non-decreasing
    leftmost_value = [frame[1] for frame in l_frames[(length(t) ÷ 2):end]]
    d_leftmost_value = diff(leftmost_value) ./ diff(t[(length(t) ÷ 2):end])
    # test leftmost is not ONLY decreasing
    low_enough = leftmost_value[2:end] .< min_activation
    decreasing = d_leftmost_value .< 0.0
    not_activated = low_enough .| decreasing
    not_activated[ismissing.(not_activated)] .= false
    return !all(not_activated)
end

 
function is_epileptic_radial_slice(exec, max_vel=50.0, min_vel=1e-2, min_traveling_frames=5, min_activation=5e-2, end_snip_idxs=3)
    # We'll call it epileptic if:
    #  1. There is a traveling front
    #  2. There is persistently elevated activity following the front.
    # How to deal with oscillatory activity? 
    # Want persistent oscillation to be considered epileptic, even if sometimes zero
    # In practice we'll ignore this for now --- FIXME
    l_frames = exec.solution.u
    ts = timepoints(exec)
    xs = [x[1] for x in space(exec).arr]
    l_frame_fronts = TravelingWaveSimulations.substantial_fronts.(l_frames, Ref(xs))
    l_persistent_fronts = persistent_fronts(l_frame_fronts, ts, max_vel)
    l_final_fronts = filter(l_persistent_fronts) do pf
        get_stop_time(pf) >= ts[end-end_snip_idxs]
    end
    E = population.(l_frames, 1)
    if (length(l_final_fronts) == 0) && all(E[end] .- E[end-1] .>= -1e-3) && all(E[end] .>= 0.05)
        # FIXME should use this to account for non-traveling case
        return true
    end
    l_activating_fronts = filter(l_final_fronts) do pf
        is_activating_front(pf)#, min_vel, min_traveling_frames)
    end
    # At least one front is not deactivated implies epilepsy
    return !all(will_be_deactivated.(l_activating_fronts, Ref(l_final_fronts)))
end

function is_traveling_solitary_radial_slice(exec, max_vel=50.0, min_vel=1e-2, min_traveling_frames=5, max_activation=5e-2, end_snip_idxs=3)
    # We'll call it traveling solitary if:
    #  1. There is a traveling front
    #  2. No elevated activity trails the front
    l_frames = exec.solution.u
    ts = timepoints(exec)
    xs = [x[1] for x in space(exec).arr]
    l_frame_fronts = TravelingWaveSimulations.substantial_fronts.(l_frames, Ref(xs))
    l_persistent_fronts = persistent_fronts(l_frame_fronts, ts, max_vel)
    l_final_fronts = filter(l_persistent_fronts) do pf
        get_stop_time(pf) >= ts[end-end_snip_idxs]
    end
    l_activating_fronts = filter(l_final_fronts) do pf
        is_activating_front(pf)#, min_vel, min_traveling_frames)
    end
    all(will_be_deactivated.(l_activating_fronts, Ref(l_final_fronts))) && has_solitary_traveling_wave(l_activating_fronts, l_final_fronts)
end 

export is_epileptic_radial_slice, is_traveling_solitary_radial_slice