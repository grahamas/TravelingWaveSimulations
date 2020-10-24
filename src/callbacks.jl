
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
