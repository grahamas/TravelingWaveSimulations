
export front_array_type, reduce_to_fronts

const front_array_type = Array{Wavefront{Float64,Float64,Value{Float64,Float64}},1}

function reduce_to_fronts(save_idxs, space)
    #error("FIXME: this function does not agree with manual reduction.")
    return if save_idxs === nothing
        (u, t, integrator) -> begin
            # Need to get x from integrator (or simulation)
            sub_u = population(u, 1)
            sub_x = [x[1] for x in space.arr]
            fronts = TravelingWaveSimulations.substantial_fronts(sub_u, sub_x)
            return fronts
          end
    else
        (u, t, integrator) -> begin
            # Need to get x from integrator (or simulation)
            sub_idx = integrator.opts.save_idxs
            sub_u = u[sub_idx]
            sub_x = [x[1] for x in space.arr[population(sub_idx,1)]]
            fronts = TravelingWaveSimulations.substantial_fronts(sub_u, sub_x)
            return fronts
          end
    end
end
