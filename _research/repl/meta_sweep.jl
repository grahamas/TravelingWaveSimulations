
function meta_sweep_from_large_run(fpath,
                                   meta_params,
                                   param_axes;
                                   property_sym=:has_propagation)
    whole_run_result = TravelingWaveSimulations.load_ExecutionClassifications(
                                                     AbstractArray, fpath
                                                    )[property_sym]
    names_to_axes = Dict(name => axis for (name, axis) in 
                         zip(
                             _namedaxisarray_names(whole_run_result), 
                             axes_keys(whole_run_result)
                            )
                        )
    sigmoid_fit_results = map(product(meta_params .|> name -> (name, names_to_axes[name]))) do (name, params)
        # FIXME do stuff
    end

end

