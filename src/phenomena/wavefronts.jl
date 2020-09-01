struct Wavefront{T,AA<:AbstractAxisArray{T,1}} <: AbstractWaveform{T,1}
# RIGHT WAVEFRONT FIXME -- It doesn't look like it's only right-side?
    left::AA
    slope::AA
    right::AA
end



left(wf::Wavefront) = wf.left
slope(wf::Wavefront) = wf.slope
right(wf::Wavefront) = wf.right
left_val(wf::Wavefront) = left(wf) |> only
slope_val(wf::Wavefront) = slope(wf) |> only
right_val(wf::Wavefront) = right(wf) |> only
Base.max(wf::Wavefront) = max(left_val(wf), right_val(wf))

# Center on peaks of first derivative


function substantial_fronts(exec::AbstractFullExecution)
    substantial_fronts.(exec.solution.u, Ref([x[1] for x in space(exec).arr])) 
end

function detect_all_fronts(values_arr::AA) where {T, AA<:AbstractAxisArray{T,1}}
    "Partition space at extrema"
    if all(values_arr .== values_arr[begin])
        return Wavefront{T,AA}[]
    end
    d_values_arr = diff(values_arr)
    dd_values_arr = diff(d_values_arr)
    values = interpolate(values_arr, Gridded(Linear()))
    d_values = interpolate(d_values_arr, Gridded(Linear()))
    dd_values = interpolate(dd_values_arr, Gridded(Linear()))
    ax_begin, ax_end = only(axes_keys(values_arr))[[begin, end]]
    ax_value_extrema = find_zeros(d_values, only(axes_keys(d_values_arr))[begin], only(axes_keys(d_values_arr))[end])
    ax_slope_extrema = find_zeros(dd_values, only(axes_keys(dd_values_arr))[begin], only(axes_keys(dd_values_arr))[end])
    slope_idx = 1
    no_inflection_point = 0
    # Should be able to assume slope extremum in every interval except first and last
    ax_left = ax_value_extrema[begin]
    wavefronts = map(ax_value_extrema[begin+1:end]) do ax_right
        ax_slope_maybe = ax_slope_extrema[slope_idx]
        ax_slope = if ax_left <= ax_slope_maybe <= ax_right
            slope_idx += 1
            ax_slope_maybe
        else
            no_inflection_point += 1
            #@warn "Inflection point not found within putative wavefront: $(no_inflection_point)"
            [ax_left, ax_right][argmax(d_values.([ax_left, ax_right]))]
        end
        Wavefront(point(values(ax_left), ax_left), point(d_values(ax_slope), ax_slope), point(values(ax_right), ax_right))
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
        slope_val(front) >= min_slope
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