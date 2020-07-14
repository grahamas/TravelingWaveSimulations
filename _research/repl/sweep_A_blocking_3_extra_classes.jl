classifications_A[:maybe_solitary_wave] = similar(classifications_A[:has_propagation]) .= classifications_A[:has_propagation] .& .!classifications_A[:persistently_active_near_origin]

classifications_A[:maybe_epileptic] = similar(classifications_A[:has_propagation]) .= classifications_A[:has_propagation] .& classifications_A[:persistently_active_near_origin]
