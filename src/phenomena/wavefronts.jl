struct Wavefront{T,AA<:AxisArray{T,1}} <: AbstractWaveform{T,1}
# RIGHT WAVEFRONT FIXME -- It doesn't look like it's only right-side?
    left::AA
    slope::AA
    right::AA
end

left(wf::Wavefront) = wf.left
slope(wf::Wavefront) = wf.slope
right(wf::Wavefront) = wf.right
_val(aa::AxisArray) = only(aa)
_loc(aa::AxisArray) = aa |> axes_keys |> only |> only
left_val(wf::Wavefront) = left(wf) |> _val
slope_val(wf::Wavefront) = slope(wf) |> _val
right_val(wf::Wavefront) = right(wf) |> _val
left_loc(wf::Wavefront) = left(wf) |> _loc
slope_loc(wf::Wavefront) = slope(wf) |> _loc
right_loc(wf::Wavefront) = right(wf) |> _loc
Base.max(wf::Wavefront) = max(left_val(wf), right_val(wf))

# Center on peaks of first derivative

function substantial_fronts(exec::AbstractFullExecution)
    substantial_fronts.(exec.solution.u, Ref([x[1] for x in space(exec).arr])) 
end

# TODO deal with multipop
function substantial_fronts(multipop::AbstractMatrix, periodic::Bool, slope_min=1e-4)
    substantial_fronts(population(multipop,1), periodic, slope_min)
end

periodic_or_neumann0_bc(T::DataType, ::Val{true}) = PeriodicBC(T)
periodic_or_neumann0_bc(T::DataType, ::Val{false}) = Neumann0BC(zero(T))

function make_ghost_op(T, axis, 
        deriv_degree::Int, periodic::Bool, deriv_order::Int=deriv_degree+1)
    step = axis[2] - axis[1]
    n = length(axis)
    D = CenteredDifference(deriv_degree, deriv_order, step, n)
    Q = periodic_or_neumann0_bc(T, Val(periodic))
    return D*Q
end

function deriv(arr::AxisVector{T}, deriv_degree::Int, periodic::Bool, deriv_order::Int=deriv_degree+1) where T
    axis = axes_keys(arr) |> only
    return AxisArray(make_ghost_op(T, axis, deriv_degree, periodic, deriv_order) * arr, axis)
end

function deriv!(darr::AxisVector, arr::AxisVector, ghost_op::GhostDerivativeOperator)
    LinearAlgebra.mul!(parent(darr), ghost_op.L, ghost_op.Q * parent(arr))
end

function detect_all_fronts(arr::AA, periodic) where {T, AA<:AxisVector{T}}
    "Partition space at extrema"
    if all(arr .== arr[begin])
        return Wavefront{T,AA}[]
    end
    locs = only(axes_keys(arr))
    d1_arr = deriv(arr, 1, periodic)
    d2_arr = deriv(arr, 2, periodic)
    arr_interp = interpolate(arr, Gridded(Linear()))
    slopes_begin_loc, slopes_end_loc = locs[[begin, end]]
    
    # For each slope extremum loc (already filtered for slope magnitude),
    # search remaining slope zero locs for next greater than extremum loc
    # If index by nothing, then two extrema locs greater than all zero locs
    # (which is mathematically impossible)
    prev_slope_zero_idx = 1
    left_bdry = slopes_begin_loc
    slope_zero_locs = [linear_find_zeros(d1_arr[begin+1:end-1])..., slopes_end_loc]
    wavefronts = map(slope_zero_locs) do right_bdry
        front = arr[left_bdry..right_bdry]
        slope_front = d1_arr[left_bdry..right_bdry]
        if length(slope_front) > 0
            slope_extremum_idx = argmax(abs.(slope_front))
            slope_extremum = slope_front[slope_extremum_idx]
            slope_extremum_loc = only(axes_keys(front))[slope_extremum_idx]
        else
            slope_extremum = 0.
            slope_extremum_loc = (left_bdry + right_bdry) / 2
        end
        left = point(arr_interp(left_bdry), left_bdry)
        left_bdry = right_bdry
        Wavefront(left,
                  point(slope_extremum, slope_extremum_loc),
                  point(arr_interp(right_bdry), right_bdry))
    end
    return wavefronts
end


function substantial_fronts(frame::AxisArray{T,1}, periodic, slope_min=1e-4) where T
    all_fronts = detect_all_fronts(frame, periodic)
    consolidated = consolidate_fronts(all_fronts, slope_min)
    return consolidated
end

eat_left(left, right) = Wavefront(left.left, right.slope, right.right)
eat_right(left, right) = Wavefront(left.left, left.slope, right.right)
function consolidate_fronts(fronts::AbstractVector{WF}, min_slope=1e-4)::AbstractVector{WF} where {WF <: Wavefront{Float64}}
    "Insist that all fronts have a minimum slope, otherwise consolidate"
    if length(fronts) == 0
        return fronts
    end

    fronts = filter(fronts) do front
        abs(slope_val(front)) >= min_slope
    end
    if length(fronts) == 0
        return WF[]
    end
    fronts[begin] = Wavefront(left(fronts[begin]), slope(fronts[begin]), right(fronts[begin]))
    for idx in eachindex(fronts[begin:end-1])
        this_front = fronts[idx]
        next_front = fronts[idx+1]
        if right(this_front) != left(next_front)
            fronts[idx] = Wavefront(left(this_front), slope(this_front), left(next_front))
        end
    end
    fronts[end] = Wavefront(left(fronts[end]), slope(fronts[end]), right(fronts[end]))
    return fronts
end