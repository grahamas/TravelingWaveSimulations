

struct Persistent{WAVE<:TravelingWaveSimulations.AbstractWaveform,T,ARR_WAVE<:AbstractArray{WAVE,1},ARR_TIME<:AbstractArray{T,1}}
    waveforms::ARR_WAVE
    t::ARR_TIME
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

function is_epileptic_radial_slice(exec, max_vel=50.0, min_vel=1e-2, min_traveling_frames=5, min_activation=5e-2)
    # We'll call it epileptic if:
    #  1. There is a traveling front
    #  2. There is persistently elevated activity following the front.
    # How to deal with oscillatory activity? 
    # Want persistent oscillation to be considered epileptic, even if sometimes zero
    # In practice we'll ignore this for now --- FIXME
    l_frames = exec.solution.u
    t = timepoints(exec)
    x = [x[1] for x in space(exec).arr]
    l_frame_fronts = TravelingWaveSimulations.substantial_fronts.(l_frames, Ref(x))
    (persistent_activation(l_frames, t, min_activation) 
        && (
            mean(population(l_frames[end],1)) > 0.2 
            || has_traveling_front(l_frame_fronts, t, max_vel, min_vel, min_traveling_frames)
        )
    )
end

function is_traveling_solitary_radial_slice(exec, max_vel=50.0, min_vel=1e-2, min_traveling_frames=5, max_activation=5e-2)
    # We'll call it traveling solitary if:
    #  1. There is a traveling front
    #  2. No elevated activity trails the front
    l_frames = exec.solution.u
    t = timepoints(exec)
    x = [x[1] for x in space(exec).arr]
    l_frame_fronts = TravelingWaveSimulations.substantial_fronts.(l_frames, Ref(x))
    !persistent_activation(l_frames, t, max_activation) && has_traveling_front(l_frame_fronts, t, max_vel, min_vel, min_traveling_frames)
end