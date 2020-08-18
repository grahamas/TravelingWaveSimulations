using Parameters
@with_kw struct R2S1IsotropicConnectivity{T} <: AbstractExpDecayingConnectivityParameter{T,3}
    local_amplitude::T
    long_range_amplitude::T
    local_R2_spread::NTuple{2,T}
    long_range_R2_spread::NTuple{2,T}
    long_range_S1_spread::NTuple{1,T}
end

function NeuralModels.directed_weights(conn::R2S1IsotropicConnectivity{T},
                                       locations::AbstractEmbeddedLattice{T,2,3}) where {T}
    diffs = differences(locations)
    diffs_R2 = [diff[1:2] for diff in diffs]
    diffs_S1 = [diff[[3]] for diff in diffs]
    step_size = step(locations)
    local_weights = directed_weights.(Ref(ExpSumSqDecayingConnectivity{T,2}), diffs_R2,
        conn.local_amplitude, Ref(conn.local_R2_spread), Ref(step_size[1:2]))
    S1_based_long_range_weights = directed_weights.(Ref(ExpSumSqDecayingConnectivity{T,1}), diffs_S1,
        1.0, Ref(conn.long_range_S1_spread), Ref(step_size[[3]]))
    R2_based_long_range_weights = directed_weights.(Ref(ExpSumSqDecayingConnectivity{T,2}), diffs_R2,
        conn.long_range_amplitude, Ref(conn.long_range_R2_spread), Ref(step_size[1:2]))
    long_range_weights = S1_based_long_range_weights .* R2_based_long_range_weights
    return local_weights + long_range_weights
end
