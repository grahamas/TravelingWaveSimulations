
function point(val::Number, coord::Number)
    AxisArray([val], [coord])
end

coordinates(ax_arr::AbstractAxisArray) = product(axes_keys(ax_arr)...)
parent_diff(arr::AbstractArray) = diff(arr)
parent_diff(ax_arr::AbstractAxisArray) = parent_diff(parent(ax_arr))
function Base.diff(ax_arr::AA) where {T, AA <: AbstractAxisArray{T,1}}
    ax_coord_1d = only(axes_keys(ax_arr))
    ax_coord_diff_1d = diff(ax_coord_1d)
    ax_midpoints = ax_coord_1d[begin:end-1] .+ (0.5 .* ax_coord_diff_1d)
    AxisArray(parent_diff(ax_arr) ./ ax_coord_diff_1d, ax_midpoints)
end

function linear_interpolate((x1,x2), (y1,y2), x)
    @assert x1 <= x <= x2 || x2 <= x <= x1
    if x1 == x2
        return y1 + (y2 - y1) / 2 # FIXME: if flat, just split the difference
    end
    y1 + (x - x1) * ((y2 - y1) / (x2 - x1)) 
end

function zero_crossing(A::AbstractAxisArray, i, j, atol=1e-4)
    if (A[i] <= 0 && A[j] >= 0) || (A[i] >= 0 && A[j] <= 0)
        return linear_interpolate_loc((getvalue(A,i),getvalue(A,j)), 0.0)
    elseif (A[i] < -atol && A[j] >= -atol) || (A[i] > atol && A[j] <= atol)
        return getvalue(A,j)
    else
        return nothing
    end
end

function linear_find_zeros(A::AbstractAxisArray)
    ((prev_ax, prev_val), remaining_ax_val_pairs) = Iterators.peel(zip(only(axes_keys(A)), A))
    zero_axs = Float64[]
    for (ax, val) in remaining_ax_val_pairs
        if val ≈ prev_val ≈ 0
            push!(zero_axs, prev_ax)
            push!(zero_axs, ax)
        elseif sign(val) != sign(prev_val)
            push!(zero_axs, linear_interpolate((prev_val, val), (prev_ax, ax), 0.0))
        end
        prev_ax, prev_val = ax, val
    end
    return zero_axs
end

function Interpolations.interpolate(data::AbstractAxisArray, args...; kwargs...)
    axs = axes_keys(data)
    interpolate(axs, data, args...; kwargs...)
end
