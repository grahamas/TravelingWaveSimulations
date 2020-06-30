export reduce_to_wave_properties

function reduce_to_wave_properties(sol)
    if sol.prob.p isa Union{NamedTuple,Dict}
        (sol.prob.p..., wave_properties=ExecutionClassifications(sol),)
    else
        @show sol.prob.p
        (wave_properties=ExecutionClassifications(sol),)
    end
end
