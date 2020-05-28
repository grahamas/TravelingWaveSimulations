
export terminate_when_E_fully_propagates

function terminate_when_E_fully_propagates(save_idxs, min_time=5)
    DiscreteCallback(E_is_fully_propagated(save_idxs, min_time), terminate!)
end

near_zero(u::T) where {T <: Number} = isapprox(u, zero(T), atol=1e-4)
near_zero(u::AbstractArray) = all(near_zero.(u))

function E_is_fully_propagated(save_idxs, min_time)
    callback = if save_idxs === nothing
        (u,t,integrator) -> begin
            pop = population(u,1)
            if t > min_time
                if (near_zero(pop) || (pop[end-10] > 0.005))
                    return true
                end
            end
            return false
        end
    else
        (u,t,integrator) -> begin
            pop = population(u[integrator.opts.save_idxs],1)
            if t > min_time
                if (near_zero(pop) || (pop[end-10] > 0.005))
                    return true
                end
            end
            return false
        end
    end
    return callback
end
