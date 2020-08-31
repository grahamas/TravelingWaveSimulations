struct Wavefront{T,AA<:AbstractAxisArray{T,1}} <: AbstractWaveform{T,1}
# RIGHT WAVEFRONT FIXME -- It doesn't look like it's only right-side?
    left::AA
    slope::AA
    right::AA
end


slope_loc(wv::Wavefront) = wv.slope.loc
function translate(wf::Wavefront, args...)
    Wavefront(
        translate(wf.left, args...),
        translate(wf.slope, args...),
        translate(wf.right, args...)
        )
end
Base.max(wf::Wavefront) = max(wf.left, wf.right)

# Center on peaks of first derivative


function substantial_fronts(exec::AbstractFullExecution)
    substantial_fronts.(exec.solution.u, Ref([x[1] for x in space(exec).arr])) 
end

function detect_all_fronts(values_arr::AbstractAxisArray{T,1}) where T
    "Partition space at extrema"
    values = interpolate(values, Gridded(Linear()))
    d_values = interpolate(diff(values_arr), Gridded(Linear()))
    dd_values = interpolate(diff(d_values), Gridded(Linear()))
    ax_begin, ax_end = only(axes_keys(values_arr))[[begin, end]]
    ax_value_extrema = find_zeros(d_values, only(axes_keys(d_values)))
    ax_slope_extrema = find_zeros(dd_values, only(axes_keys(dd_values)))
    slope_idx = 0 
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
            @warn "Inflection point not found within putative wavefront: $(no_inflection_point)"
            [ax_left, ax_right][argmax(d_values.([ax_left, ax_right]))]
        end
        Wavefront(values[ax_left..ax_left], values[ax_slope..ax_slope], values[ax_right..ax_right])
    end
    return wavefronts
end

eat_left(left, right) = Wavefront(left.left, right.slope, right.right)
eat_right(left, right) = Wavefront(left.left, left.slope, right.right)
function consolidate_fronts(fronts::AbstractVector{WF}, min_slope=1e-4) where {WF <: Wavefront{Float64}}
    "Insist that all fronts have a minimum slope, otherwise consolidate"
    val_begin = left(fronts[begin])
    val_end = right(fronts[end])

    fronts = filter(fronts) do front
        slope_val(front) >= min_slope
    end
    if length(fronts) == 0
        return WF[]
    end
    fronts[begin] .= Wavefront(val_begin, slope(fronts[begin]), right(fronts[begin]))
    for idx in eachindex(fronts[begin:end-1])
        this_front = fronts[idx]
        next_front = fronts[idx+1]
        if right(this_front) != left(next_front)
            fronts[idx] .= Wavefront(left(this_front), slope(this_front), left(next_front))
        end
    end
    fronts[end] .= Wavefront(left(fronts[end]), slope(fronts[end]), val_end)
end

# TODO deal with multipop
function substantial_fronts(multipop::AbstractMatrix, slope_min=1e-4)
    substantial_fronts(population(multipop,1), slope_min)
end
function substantial_fronts(frame::AbstractAxisArray{T,1}, slope_min=1e-4)
    all_fronts = detect_all_fronts(frame)
    consolidated = consolidate_fronts(all_fronts, slope_min)
    return consolidated
end