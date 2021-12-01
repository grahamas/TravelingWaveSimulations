
function terminate_when_E_fully_propagates(save_idxs; proportion_full=0.95, min_time=0.1, zero_tol=1e-3)
    DiscreteCallback(E_is_fully_propagated(save_idxs, proportion_full, min_time, zero_tol), terminate!)
end

near_zero(u::T, zero_tol::T) where {T <: Number} = abs(u) < zero_tol
# near_zero(u::AbstractArray, zero_tol) = all(near_zero.(u, zero_tol))

function E_is_fully_propagated(save_idxs, proportion_full, min_time, zero_tol)
    callback = if save_idxs === nothing
        (u,t,integrator) -> begin
            pop = population(u,1)
            fully_propagated_dx = floor(Int, proportion_full * length(pop))
            if t > min_time
                if (all(near_zero.(pop, zero_tol)) || any(.!near_zero.(pop[fully_propagated_dx:end], zero_tol)))
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
                if (all(near_zero.(pop, zero_tol)) || any(.!near_zero.(pop[fully_propagated_dx:end], zero_tol)))
                    return true
                end
            end
            return false
        end
    end
    return callback
end

function is_spread(u, t, integrator)
    save_idxs = integrator.opts.save_idxs
    pop = save_idxs === nothing ? population(u, 1) : population(u[integrator.opts.save_idxs],1)
    max_spread_dx = floor(Int, integrator.p.max_spread_proportion * length(pop))
    if t > integrator.p.min_spread_time
        if (all(near_zero.(pop, integrator.p.max_spread_value)) || any(.!near_zero.(pop[max_spread_dx:end], integrator.p.max_spread_value)))
            return true
        end
    end
    return false
end
export is_spread

function Simulation73.handle_callback(sim::Simulation, cb::CB) where {CB <: Tuple{typeof(is_spread), <:NamedTuple}}
    fn, nt = cb
    return (nt, DiscreteCallback(fn, terminate!))
end

# is_propagated defined in phenomena/classifications
# function Simulation73.handle_callback(sim::Simulation{T}, cb::CB) where {T,CB <: Tuple{typeof(is_propagated), <:NamedTuple}}
#     fn, nt = cb
#     min_dist_for_propagation = nt.min_dist_for_propagation
#     has_traveled_dist(rf::RunningFront) = abs(rf.previous_location - rf.starting_location) >= min_dist_for_propagation
#     has_traveled_dist(::Nothing) = false
#     xs = frame_xs(sim)
#     d1_ghost_op = make_ghost_op(T, xs, 1, sim.space isa PeriodicLattice)
#     dframe_cache = rand(size(xs)...)
#     running_fronts = RF_possibilities(T)[nothing for _ in xs]
#     p = merge((running_fronts=running_fronts, dframe_cache=dframe_cache, return_fn=has_traveled_dist, d1_ghost_op=d1_ghost_op, has_propagation=Union{Missing,Bool}[missing]), nt)
#     cb = DiscreteCallback(fn, terminate!)
#     return (p, cb)
# end
