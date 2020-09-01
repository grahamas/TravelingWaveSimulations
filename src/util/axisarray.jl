coordinates(ax_arr::AbstractAxisArray) = product(axes_keys(ax_arr)...)
parent_diff(arr::AbstractArray) = diff(arr)
parent_diff(ax_arr::AbstractAxisArray) = parent_diff(parent(ax_arr))
function Base.diff(ax_arr::AA) where {T, AA <: AbstractAxisArray{T,1}}
    ax_vals = only(axes_keys(ax_arr))
    coord_diff_1d = diff(ax_vals)
    ax_midpoints = ax_vals[begin:end-1] .+ (0.5 .* coord_diff_1d)
    AxisArray(parent_diff(ax_arr) ./ coord_diff_1d, ax_midpoints)
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

function Interpolations.interpolate(data::AbstractAxisArray, args...; kwargs...)
    axs = axes_keys(data)
    interpolate(axs, data, args...; kwargs...)
end

function point(val::Number, coord::Number)
    AxisArray([val], [coord])
end