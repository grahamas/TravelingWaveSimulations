
using Interpolations, Optim, DiffEqOperators

function Interpolations.interpolate(data::AbstractAxisArray, args...; kwargs...)
    axs = axes_keys(data)
    interpolate(axs, data, args...; kwargs...)
end

get_coordinates(aaa::AbstractAxisArray) = product(keys.(axes(aaa))...) |> collect

l2(x,y) = sqrt(sum((x .- y) .^ 2))
function diagonal_slice(x_axis, y_axis, data::Matrix, y_intercept, slope, dx=1)
    interp = LinearInterpolation((x_axis, y_axis), data')
    calc_y(x) = slope * x + y_intercept
    sample_points = [(x, calc_y(x)) for x in x_axis if y_axis[begin] <= calc_y(x) <= y_axis[end]]
    sample_values = sample_points .|> (x) -> interp(x...)
    distance_along_line = l2.(Ref(sample_points[begin]), sample_points)
    return (distance_along_line, sample_points, sample_values)
end

function project(coord::C, line::C, line_point::C) where {N, C <: SVector{N}}
    projection = (line * line') / (line' * line)
    return projection * coord + (I - projection) * line_point
end

line_dist_to_coord(dist, line, origin) = (dist * line) + origin 
#axes_vals(data::AbstractAxisArray) = keys.(axes(data))

struct PointVectorLine{N,T,S<:SVector{N,T}}
    point::S
    vector::S
end
slope(line::PointVectorLine{2}) = line.vector[2] / line.vector[1]
# the point where line takes value val in dim dimension
function point_from_dim_val(line::PointVectorLine, val::Number, dim::Int)
    scale = (val - line.point[dim]) / line.vector[dim]
    return line.point .+ line.vector .* scale
end
x_from_y(line::PointVectorLine{2}, y::Number) = point_from_dim_val(line, y, 2)[1]
y_from_x(line::PointVectorLine{2}, x::Number) = point_from_dim_val(line, x, 1)[2]
y_intercept(line::PointVectorLine{2}) = y_from_x(line, 0.0)
x_intercept(line::PointVectorLine{2}) = x_from_y(line, 0.0) 
point_from_distance(line::PointVectorLine, dist::Number) = line.vector .* dist .+ line.point
get_orthogonal_vector(line::PointVectorLine) = -SA[line.vector[2], -line.vector[1]]

function originate_from_left(line, xs, ys)
    x_min, x_max = extrema(xs)
    y_min, y_max = extrema(ys)

    # need increasing x
    new_vector = line.vector[1] > 0 ? line.vector : -line.vector

    if y_min <= y_from_x(line, x_min) <= y_max
        return PointVectorLine(SA[x_min, y_from_x(line, x_min)], new_vector)
    elseif x_min <= x_from_y(line, y_max) <= x_from_y(line, y_min)
        @assert x_from_y(line, y_max) <= x_max
        return PointVectorLine(SA[x_from_y(line, y_max), y_max], new_vector)
    elseif x_min <= x_from_y(line, y_min) <= x_from_y(line, y_max)
        @assert x_from_y(line, y_min) <= x_max
        return PointVectorLine(SA[x_from_y(line, y_min), y_min], new_vector)
    else
        @show line
        @show x_from_y(line, y_min)
        @show x_from_y(line, y_max)
        @show y_from_x(line, x_min)
        @show y_from_x(line, x_max)
        @show extrema(xs)
        @show extrema(ys)
        error("Line does not fall within axes")
    end
end

squish(data::NamedDimsArray, args...) = (@show "squishing 2"; squish(data.data, args...))
slice(data::NamedDimsArray, args...) = (@show "Slicing 3"; slice(data.data, args...))
reduce_along_max_central_gradient(data::NamedDimsArray, args...) = (@show "Reducing along max grad"; reduce_along_max_central_gradient(data.data, args...))
function slice(data::AbstractAxisArray{T,2}, target_line::PointVectorLine{2,T,S}, step_fineness=5) where {T,S}
    xs, ys = axes_keys(data)
    interpolation = interpolate(data, Gridded(Linear()))
    return slice(interpolation, target_line, xs, ys, step_fineness)
end

function slice(interpolation::AbstractInterpolation, target_line::PointVectorLine{2,T,S}, xs, ys, step_fineness=5) where {T,S}
    # We want this line to proceed from its point through the grid (left to right for convention's and plotting's sake) 
    line = originate_from_left(target_line, xs, ys)

    # Calculate a reasonable step
    dx = abs(xs[begin+1] - xs[begin])
    dy = abs(xs[begin+1] - xs[begin])
    step = sqrt(dx^2 + dy^2) / step_fineness

    # Step along line until
    distance_along_line = 0
    dists = T[]; points = S[]; values = T[];
    next_point = point_from_distance(line, distance_along_line)
    while xs[begin] <= next_point[1] <= xs[end] && ys[begin] <= next_point[2] <= ys[end]
        next_value = interpolation(next_point...)
        # FIXME preallocate
        push!(dists, distance_along_line)
        push!(points, next_point)
        push!(values, next_value)
        distance_along_line += step
        next_point = point_from_distance(line, distance_along_line)
    end
    return (dists, values, points, line) 
end

function squish(data::AbstractAxisArray{T,2}, target_line::PointVectorLine{2,T,S}, args...) where {T,S}
    xs, ys = axes_keys(data)
    interpolation = interpolate(data, Gridded(Linear()))
    return squish(interpolation, target_line, xs, ys, args...)
end
function squish(interpolation::AbstractInterpolation, target_line::PointVectorLine{2,T,S}, xs, ys, line_fineness=5, squish_fineness=2) where {T,S}
    squished_dists, _, squished_points, squished_line = slice(interpolation, target_line, xs, ys, line_fineness)
    orthogonal_vec = get_orthogonal_vector(target_line)
    squished_vals = map(squished_points) do point
        _, crosscut_vals, _ = slice(interpolation, PointVectorLine(point, orthogonal_vec), xs, ys, squish_fineness)
        mean(crosscut_vals)
    end
    return squished_dists, squished_vals, squished_points, squished_line
end

function midpoint(arr::AbstractArray)
    @assert arr[begin] < arr[end]
    return (arr[end] - arr[begin]) / 2 + arr[begin]
end
function reduce_along_max_central_gradient(data::AbstractAxisArray{T,2}, reduction::Function=slice, line_fineness=5) where T
    xs, ys = axes_keys(data)
    center_pt = [midpoint(xs), midpoint(ys)] 
    interpolation = interpolate(data, Gridded(Linear()))
    neg_gradient_norm(coord) = -norm(Interpolations.gradient(interpolation, coord...))
    result_optim = optimize(neg_gradient_norm, [xs[begin], ys[begin]], 
                                               [xs[end], ys[end]], 
                                               center_pt, 
                                               Fminbox(GradientDescent()); 
                                               autodiff=:forward)
    max_grad_coord = Optim.minimizer(result_optim)
    max_grad = Interpolations.gradient(interpolation, max_grad_coord...)
    max_grad_line = PointVectorLine(SA[max_grad_coord...], max_grad)
    @show string(reduction)
    return reduction(interpolation, max_grad_line, xs, ys, line_fineness)
end

struct FittedSigmoid{T}
    left_val::T
    change::T
    slope::T
    threshold::T
    error::T
end

(s::FittedSigmoid)(x) = s.left_val .+ s.change .* NeuralModels.simple_sigmoid_fn(x, s.slope, s.threshold)

function fit_sigmoid(ys, xs)
    # fits sigmoid that goes from 0 to 1
    if xs[end] < xs[begin]
        xs = reverse(xs)
        ys = reverse(ys)
    end
    if xs[begin] == xs[end]
        return nothing
    end
    left_val = ys[begin]
    sigmoid_change = ys[end] - ys[begin]
    sigmoid_fn(params) = left_val .+ sigmoid_change .* NeuralModels.simple_sigmoid_fn.(xs, params[1], params[2])
    loss(params) = l2(ys, sigmoid_fn(params))
    D = CenteredDifference{1}(1, 2, xs[2] - xs[1], length(ys))
    derivative_estimate = D*ys
    slope_est, i_max = findmax(abs.(derivative_estimate))
    threshold_est = xs[i_max]

    lower_param_bounds = [slope_est/5, xs[begin]]
    upper_param_bounds = [5*slope_est, xs[end]]
    fit_result = optimize(loss, lower_param_bounds, upper_param_bounds,
                          [slope_est, threshold_est], Fminbox(NelderMead());
                          autodiff=:forward)
    slope_fit, threshold_fit = Optim.minimizer(fit_result)
    error_fit = Optim.minimum(fit_result) |> abs
    return FittedSigmoid(
                         left_val,
                         sigmoid_change,
                         slope_fit,
                         threshold_fit,
                         error_fit
                        )
end

function phase_space_sigmoid_fit(uncollapsed_phase_space, (x_axis_name, y_axis_name), reduction::Function=squish)
    (x_axis, y_axis, phase_space) = _collapse_to_axes(uncollapsed_phase_space, x_axis_name, y_axis_name)
    line_dists, line_vals, line_locs, line = reduce_along_max_central_gradient(phase_space, reduction)
    fit_sigmoid(line_vals, line_dists)
end



