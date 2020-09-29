struct Wavefront{T,AA<:AbstractAxisArray{T,1}} <: AbstractWaveform{T,1}
# RIGHT WAVEFRONT FIXME -- It doesn't look like it's only right-side?
    left::AA
    slope::AA
    right::AA
end

left(wf::Wavefront) = wf.left
slope(wf::Wavefront) = wf.slope
right(wf::Wavefront) = wf.right
_val(aa::AbstractAxisArray) = only(aa)
_loc(aa::AbstractAxisArray) = aa |> axes_keys |> only |> only
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

function deriv_periodic_bc(arr::AbstractAxisArray{T,1}, deriv_degree::Int, deriv_order::Int=deriv_degree+1) where T
    axis = axes_keys(arr) |> only
    step = axis[2] - axis[1]
    n = length(axis)
    D = CenteredDifference(deriv_degree, deriv_order, step, n)
    Q = PeriodicBC(T)
    return AxisArray(D*Q*arr, axis)
end

function detect_all_fronts(values_arr::AA) where {T, AA<:AbstractAxisArray{T,1}}
    "Partition space at extrema"
    if all(values_arr .== values_arr[begin])
        return Wavefront{T,AA}[]
    end
    d_values_arr = deriv_periodic_bc(values_arr, 1)
    dd_values_arr = deriv_periodic_bc(values_arr, 2)
    values = interpolate(values_arr, Gridded(Linear()))
    d_values = interpolate(d_values_arr, Gridded(Linear()))
    slopes_begin_loc, slopes_end_loc = only(axes_keys(d_values_arr))[[begin, end]]
    slope_zero_locs = [slopes_begin_loc, linear_find_zeros(d_values_arr)..., slopes_end_loc]
    slope_extrema_locs = linear_find_zeros(dd_values_arr)
    slope_idx = 1
    no_inflection_point = 0
    # Should be able to assume slope extremum in every interval except first and last
    # Partition on slope_zero_locs, beginning and ending with boundaries
    # Need maximum slope between
    front_left_loc = slope_zero_locs[begin]
    wavefronts = map(slope_zero_locs[begin+1:end]) do front_right_loc
        slope_extremum_loc_maybe = slope_idx <= length(slope_extrema_locs) ? slope_extrema_locs[slope_idx] : NaN
        slope_extremum_loc = if !isnan(slope_extremum_loc_maybe) && front_left_loc <= slope_extremum_loc_maybe <= front_right_loc
            slope_idx += 1
            slope_extremum_loc_maybe
        else
            no_inflection_point += 1
            #@warn "Inflection point not found within putative wavefront: $(no_inflection_point)"
            [front_left_loc, front_right_loc][argmax(d_values.([front_left_loc, front_right_loc]))]
        end
        Wavefront(point(values(front_left_loc), front_left_loc), point(d_values(slope_extremum_loc), slope_extremum_loc), point(values(front_right_loc), front_right_loc))
    end
    return wavefronts
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

# TODO deal with multipop
function substantial_fronts(multipop::AbstractMatrix, slope_min=1e-4)
    substantial_fronts(population(multipop,1), slope_min)
end
function substantial_fronts(frame::AbstractAxisArray{T,1}, slope_min=1e-4) where T
    all_fronts = detect_all_fronts(frame)
    consolidated = consolidate_fronts(all_fronts, slope_min)
    return consolidated
end