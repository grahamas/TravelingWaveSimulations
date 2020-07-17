function reduce_to_wave_properties(sol; kwargs...)
    if sol.prob.p isa Union{NamedTuple,Dict}
        (sol.prob.p..., wave_properties=ExecutionClassifications(sol; kwargs...),)
    else
        @show sol.prob.p
        (wave_properties=ExecutionClassifications(sol; kwargs...),)
    end
end
