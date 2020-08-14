
using Interpolations

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
function Interpolations.interpolate(data::AbstractAxisArray, opts...)
    axs = axes_keys(data)
    interpolate(axs..., data, opts...)
end

struct Line{N,T,S<:SVector{N,T}}
    point::S
    vector::S
end
slope(line::Line{2}) = line.vector[2] / line.vector[1]
y_intercept(line::Line{2}) = line.point[2] - slope(line) * line.point[1]
y_from_x(line::Line{2}, x::Number) = slope(line) * x + y_intercept(line)
point_from_distance(line::Line, dist::Number) = line.vector .* dist .+ line.point

function originate_from_left(line, xs, ys)
    x_min, x_max = extrema(xs)
    y_min, y_max = extrema(ys)

    # need increasing x
    new_vector = line.vector[1] > 0 ? line.vector : -line.vector

    if y_min <= y_from_x(line, x_min) <= y_max
        return Line(SA[x_min, y_from_x(line, x_min)], new_vector)
    elseif y_min <= y_from_x(line, x_max) <= y_max
        return Line(SA[x_max, y_from_x(line, x_max)], new_vector)
    else
        error("Line does not fall within axes")
    end
end

slice(data::NamedDimsArray, args...) = squish(data.data, args...)
function slice(data::AbstractAxisArray{T,2}, target_line::Line{2,T,S}, step_fineness=5) where {T,S}
    xs, ys = axes_keys(data)
    interpolation = LinearInterpolation(data, extrapolation_bc=NaN)
    return slice(interpolation, target_line, xs, ys, step_fineness)
end

function slice(interpolation::AbstractInterpolation, target_line::Line{2,T,S}, xs, ys, step_fineness=5) where {T,S}
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
    next_value = interpolation(next_point)
    while next_value != NaN
        push!(dists, distance_along_line)
        push!(points, next_point)
        push!(values, next_value)
        distance_along_line += step
        next_point = point_from_distance(line, distance_along_line)
        next_value = interpolation(next_point)
    end
    return (dists, values, points) 
end

function squish(data::AbstractAxisArray{T,2}, target_line::Line{2,T,S}, line_fineness=5, squish_fineness=2) where {T,S}
    xs, ys = axes_keys(data)
    interpolation = LinearInterpolation(data, extrapolation_bc=NaN)
    squished_dists, _, squished_points = slice(interpolation, target_line, xs, ys, line_fineness)
    orthogonal_vec = get_orthogonal_vector(target_line)
    squished_vals = map(squished_points) do point
        _, crosscut_vals, _ = slice(interpolation, Line(point, orthogonal_vec), xs, ys, squish_fineness)
        mean(crosscut_vals)
    end
    return squished_dists, squished_vals, squished_points
end
