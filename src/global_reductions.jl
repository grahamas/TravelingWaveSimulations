function reduce_to_wave_properties(sol; classification_kwargs...)
    if sol.prob.p isa Union{NamedTuple,Dict}
        (sol.prob.p..., wave_properties=ExecutionClassifications(sol; classification_kwargs...),)
    else
        @show sol.prob.p
        (wave_properties=ExecutionClassifications(sol; classification_kwargs...),)
    end
end
