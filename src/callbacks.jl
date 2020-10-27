
function terminate_when_E_fully_propagates(save_idxs; proportion_full=0.95, min_time=0.1, zero_tol=1e-3)
    DiscreteCallback(E_is_fully_propagated(save_idxs, proportion_full, min_time, zero_tol), terminate!)
end

near_zero(u::T, zero_tol::T) where {T <: Number} = abs(u) < zero_tol
near_zero(u::AbstractArray, zero_tol) = all(near_zero.(u, zero_tol))

function E_is_fully_propagated(save_idxs, proportion_full, min_time, zero_tol)
    callback = if save_idxs === nothing
        (u,t,integrator) -> begin
            pop = population(u,1)
            fully_propagated_dx = floor(Int, proportion_full * length(pop))
            if t > min_time
                if (near_zero(pop, zero_tol) || !near_zero(pop[fully_propagated_dx], zero_tol))
                    return true
                end
            end
            return false
        end
    else
        (u,t,integrator) -> begin
            pop = population(u[integrator.opts.save_idxs],1)
            fully_propagated_dx = floor(Int, proportion_full * length(pop))
            if t > min_time
                if (near_zero(pop, zero_tol) || !near_zero(pop[fully_propagated_dx], zero_tol))
                    return true
                end
            end
            return false
        end
    end
    return callback
end


# is_propagated defined in phenomena/classifications
function Simulation73.handle_callback(sim::Simulation{T,M,S,IV,ALG,DT,SV_IDX,CB,GR}) where {T,M,S,IV,ALG,DT,SV_IDX,GR, CB <: Tuple{typeof(is_propagated), <:NamedTuple}}
    fn, nt = sim.callback
    min_dist_for_propagation = nt.min_dist_for_propagation
    has_traveled_dist(rf::RunningFront) = abs(rf.previous_location - rf.starting_location) >= min_dist_for_propagation
    has_traveled_dist(::Nothing) = false
    xs = frame_xs(sim)
    d1_ghost_op = make_ghost_op(T, xs, 1, sim.space isa PeriodicLattice)
    dframe_cache = AxisArray(rand(size(xs)...), xs)
    running_fronts = AxisVector{RF_RFM_NOTHING}(RF_RFM_NOTHING[nothing for _ in xs], collect(xs))
    p = merge((running_fronts=running_fronts, dframe_cache=dframe_cache, return_fn=has_traveled_dist, d1_ghost_op=d1_ghost_op, has_propagation=[false]), nt)
    cb = DiscreteCallback(fn, terminate!)
    return (p, cb)
end