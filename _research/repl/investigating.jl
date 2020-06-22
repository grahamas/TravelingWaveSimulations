function maybe_sw(example_name; kwargs...)
    example = get_example(example_name)
    name, exec = execute_single_modification(example, kwargs)
    classifications = ExecutionClassifications(exec)
    return classifications.has_propagation & !classifications.persistently_active_near_origin
end
