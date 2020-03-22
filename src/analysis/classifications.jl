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


function has_traveling_front(frame_fronts::AbstractArray{<:AbstractArray{WF}}, ts::AbstractArray{<:Number}, max_vel=20, min_vel=1e-3, min_traveling_frames=5) where {WF <:TravelingWaveSimulations.Wavefront}
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

function all_fronts(exec::AbstractFullExecution)
    TravelingWaveSimulations.substantial_fronts.(exec.solution.u, Ref([x[1] for x in space(exec).arr])) 
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

struct ID{OBJ}
    obj::OBJ
    id::Int
end
get_id(id_obj::ID) = id_obj.id

using Colors, AxisIndices
function TravelingWaveSimulations.custom_animate(execution::Execution{T,<:Simulation{T}}, p_fronts::AbstractArray{<:Persistent{WF}}; kwargs...) where {T,WF}
    ts = timepoints(execution)
    if ts[end-1] == ts[end]
        ts = ts[1:end-1]
    end
    n_p_fronts = length(p_fronts)
    cmap = distinguishable_colors(n_p_fronts+1, RGB(1,1,1), dropseed=true)
    fronts_by_time = AxisIndicesArray(Array{ID{WF},1}[ID{WF}[] for _ in ts], (ts,))
    for (i_p_front, p_front) in enumerate(p_fronts)
         for (wf, t) in zip(p_front.waveforms, p_front.t)
            push!(fronts_by_time[t], ID{WF}(wf, i_p_front))
        end
    end
    solution = execution.solution
    pop_names = execution.simulation.model.pop_names
    xs = space(execution)
    x1 = [xcoord[1] for xcoord in xs.arr]
    max_val = maximum(solution)
	min_val = minimum(solution)
    i_pop = 1
    
    @animate for time_dx in 1:length(ts) # TODO @views
        frame = population_timepoint(solution, i_pop, time_dx)
        vs = ValuedSpace(frame, x1)
        plot(
            xs, frame; label=pop_names[i_pop],
            val_lim=(min_val,max_val), title="t = $(round(ts[time_dx], digits=4))",
            xlab = "Space (a.u. approx. um)",kwargs...
            )
        wf_arr = fronts_by_time[time_dx]
        if length(wf_arr) > 0
            @show get_id.(wf_arr)
            scatter!([wf.obj.slope.loc for wf in wf_arr], [vs[wf.obj.slope.loc].val for wf in wf_arr], color=[cmap[get_id(wf)] for wf in wf_arr])
        end
    end

end


function is_activating_front(pf::Persistent, max_background_amp)
    # Either the vel and slope have opposing signs
    # or it's a static front
    # FIXME make sure big enough activation
    vel_sgn = approx_sign(mean(get_velocities(pf)))
    slope_sgn = get_slope_sign(pf)
    if vel_sgn * slope_sgn <= 0
        return true
    else
        return false
    end
end

function will_be_overtaken(persistent_front::Persistent, contemporary_fronts::AbstractArray{<:Persistent}, max_background_amp)
    pf_velocity, vel_err = estimate_velocity(persistent_front)
    pf_velocity_sign = approx_sign(pf_velocity)
    pf_stop_place = get_stop_place(persistent_front)
    for front in contemporary_fronts
        front_velocity = mean(get_velocities(front))
        if ((abs(front_velocity) > abs(pf_velocity)) # front is moving faster
                && (sign(pf_stop_place - get_stop_place(front)) == approx_sign(front_velocity))) # front velocity dir matches relative displacement
            if pf_velocity_sign == 1 && front.waveforms[end].left.val < max_background_amp
                return true
            elseif pf_velocity_sign == -1 && front.waveforms[end].right.val < max_background_amp
                return true # the persistent_front will be deactivated by front overtaking
            end
        end
    end
    return false
end

is_solitary(::Nothing, ::Any) = false
function is_solitary(subject::Persistent, contemporary_fronts::AbstractArray{<:Persistent}, max_background_amp, middle_portion_to_compare=0.3, noise_tol=1e-3)
    if length(subject) * middle_portion_to_compare < 4 
        return false # not enough points to compare
        # We want the *shorter-duration* front to have the minimum number of points
    end
    cut_portions = (1.0 - middle_portion_to_compare) / 2
    subject_middle_idx = floor(Int, length(subject) * cut_portions):ceil(Int, length(subject) * (cut_portions+middle_portion_to_compare))
    middle_ts = subject.t[subject_middle_idx]
    # So now we only compare other fronts that began before this front
    candidate_fronts = filter(contemporary_fronts) do front
        (get_start_time(front) <= get_start_time(subject) 
            && get_stop_time(front) >= middle_ts[end])            
    end
    subject_middle_pf = Persistent(subject.waveforms[subject_middle_idx], middle_ts)
    subject_middle_velocities = get_velocities(subject_middle_pf)
    subject_velocity_sign = approx_sign(subject_middle_velocities[1])
    if !all(subject_velocity_sign .== approx_sign.(subject_middle_velocities))
        return false
    end
    for candidate in candidate_fronts
        candidate_compare_idx = middle_ts[1] .<= candidate.t .<= middle_ts[end]
        candidate_compare_pf = Persistent(candidate.waveforms[candidate_compare_idx], middle_ts)
        candidate_compare_velocities = get_velocities(candidate_compare_pf)
        if all(isapprox.(candidate_compare_velocities, subject_middle_velocities, atol=noise_tol))
            # matching velocity
            if subject_velocity_sign == 1 # moving right
                if subject.waveforms[end].slope.loc > candidate.waveforms[end].slope.loc # subject leading
                    if candidate.waveforms[end].left.val < max_background_amp
                        return true
                    end
                else
                    if subject.waveforms[end].left.val < max_background_amp
                        return true
                    end
                end
            else # moving left
                if subject.waveforms[end].slope.loc < candidate.waveforms[end].slope.loc # subject leading
                    if candidate.waveforms[end].right.val < max_background_amp
                        return true
                    end
                else
                    if subject.waveforms[end].right.val < max_background_amp
                        return true
                    end
                end
            end
        end     
        # FIXME: should account for decreasing trailing front
    end
    return false
end

function has_solitary_traveling_wave(l_activating_fronts, l_final_fronts)
    any(map(filter(is_traveling, l_activating_fronts)) do traveling_front
             will_be_cancelled(traveling_front, l_final_fronts)
        end
    )
end

is_decaying(::Nothing) = false
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

function will_be_deactivated(persistent_front::Persistent, contemporary_fronts::AbstractArray{<:Persistent}, max_background_amp)
    # Check if it's decaying
    if is_decaying(persistent_front)
        return true
    end
    
    # Find a deactivate-to-zero that will overtake it
    if will_be_overtaken(persistent_front, contemporary_fronts, max_background_amp)
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

function get_farthest_traveling_front(arr_pfronts::Array{<:Persistent,1}, min_dist, min_vel, min_traveling_frames)
    arr_traveling_pfronts = filter(x -> is_traveling(x, min_vel, min_traveling_frames), arr_pfronts)
    arr_lengths = map(arr_traveling_pfronts) do pfront
        abs(pfront.waveforms[1].slope.loc - pfront.waveforms[end].slope.loc)
    end
    if length(arr_lengths) == 0
        return nothing
    end
    farthest_length, idx = findmax(arr_lengths)
    if farthest_length < min_dist
        return nothing
    end
    return arr_pfronts[idx]
end

function get_wave_properties(l_frame_fronts::Array{<:Array{<:Wavefront{T,T,Value{T,T}}}}, ts::Array{T,1}; max_vel=50.0, min_vel=1e-2, min_traveling_frames=5, max_background_amp=5e-2, end_snip_idxs=3) where T
    min_dist = min_vel * min_traveling_frames
    
    l_persistent_fronts = persistent_fronts(l_frame_fronts, ts, max_vel)
    l_final_fronts = filter(l_persistent_fronts) do pf
        get_stop_time(pf) >= ts[end-end_snip_idxs]
    end
    l_final_activating_fronts = filter(l_final_fronts) do pf
        is_activating_front(pf,max_background_amp)#, min_vel, min_traveling_frames)
    end
    
    all_final_fronts_will_die = all(will_be_deactivated.(l_final_activating_fronts, Ref(l_final_fronts), max_background_amp))
    
    # We'll call it traveling solitary if:
    #  1. There is a traveling front
    #  2. No elevated activity trails *that* front
    maybe_farthest_traveling_front = get_farthest_traveling_front(l_persistent_fronts, min_dist, min_vel, min_traveling_frames)
    #b_traveling_solitary = is_solitary(maybe_farthest_front, l_persistent_fronts)
    b_traveling_solitary = any(map(l_persistent_fronts) do front
        is_traveling(front,min_vel, min_traveling_frames) && is_solitary(front, l_persistent_fronts, max_background_amp)
    end)
    
    # We'll call it epileptic if:
    #  1. There is a traveling front
    #  2. There is persistently elevated activity following the front.
    # How to deal with oscillatory activity? 
    # Want persistent oscillation to be considered epileptic, even if sometimes zero
    # In practice we'll ignore this for now --- FIXME
    # FIXME Needs is_traveling
    b_epileptic = (maybe_farthest_traveling_front !== nothing) && !all_final_fronts_will_die
    
    b_decaying = is_decaying(maybe_farthest_traveling_front)
    velocity, velocity_error = estimate_velocity(maybe_farthest_traveling_front)
    
    WaveProperties(; 
        epileptic=b_epileptic, 
        traveling_solitary=b_traveling_solitary, 
        decaying=b_decaying,
        velocity=velocity,
        velocity_error=velocity_error
    )
end

function get_wave_properties(exec::Execution; params...)
    l_frames = exec.solution.u
    ts = timepoints(exec)
    xs = [x[1] for x in space(exec).arr]
    l_frame_fronts = TravelingWaveSimulations.substantial_fronts.(l_frames, Ref(xs))
    get_wave_properties(l_frame_fronts, ts; params...)
end

function get_wave_properties(exec::Union{AugmentedExecution,ReducedExecution{<:Wavefront}}; params...)
    get_wave_properties(exec.saved_values.saveval, exec.saved_values.t; params...)
end

function get_wave_properties(nt::NamedTuple; params...)
    get_wave_properties(nt.wavefronts, nt.wavefronts_t; params...)
end

struct WaveProperties
    epileptic::Bool
    traveling_solitary::Bool
    decaying::Bool
    velocity::Union{Float64,Missing}
    velocity_error::Union{ Float64,Missing}
end
WaveProperties(; epileptic, traveling_solitary, decaying, velocity, velocity_error) = WaveProperties(epileptic, traveling_solitary, decaying, velocity, velocity_error)