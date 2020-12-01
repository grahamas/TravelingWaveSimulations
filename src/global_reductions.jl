function reduce_to_wave_properties(sol; classification_kwargs...)
    if sol.prob.p isa Union{NamedTuple,Dict}
        (sol.prob.p..., wave_properties=ExecutionClassifications(sol; classification_kwargs...),)
    else
        (wave_properties=ExecutionClassifications(sol; classification_kwargs...),)
    end
end


function reduce_to_min_propagation_cls(sol; kwargs...)
    if sol.prob.p isa Union{NamedTuple,Dict}
        (sol.prob.p..., propagation=MinimalPropagationClassification(sol; kwargs...),)
    else
        (propagation=MinimalPropagationClassification(sol; kwargs...),)
    end
end

export already_reduced_to_min_propagation_cls
function already_reduced_to_min_propagation_cls(sol)
    (sol.prob.p..., propagation=MinimalPropagationClassification(sol.prob.p.has_propagation |> only))
end