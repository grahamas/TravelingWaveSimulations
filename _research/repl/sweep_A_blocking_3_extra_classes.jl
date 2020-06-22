classifications_A[:maybe_solitary_wave] = mod_array(Bool) .= classifications_A[:has_propagation] .& .!classifications_A[:persistently_active_near_origin]

classifications_A[:maybe_epileptic] = mod_array(Bool) .= classifications_A[:has_propagation] .& classifications_A[:persistently_active_near_origin]
