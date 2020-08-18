function maybe_sw(prototype_name; kwargs...)
    prototype = get_prototype(prototype_name)
    name, exec = execute_single_modification(prototype, kwargs)
    classifications = ExecutionClassifications(exec)
    return classifications.has_propagation & !classifications.persistently_active_near_origin
end
