
using Interpolations

l2(x,y) = sqrt(sum((x .- y) .^ 2))
function diagonal_slice(x_axis, y_axis, data::Matrix, y_intercept, slope, dx=1)
    interp = LinearInterpolation((x_axis, y_axis), data')
    calc_y(x) = slope * x + y_intercept
    sample_points = [(x, calc_y(x)) for x in x_axis if y_axis[begin] <= calc_y(x) <= y_axis[end]]
    @show first(sample_points)
    sample_values = sample_points .|> (x) -> interp(x...)
    @show (sample_points[begin], sample_points[end])
    distance_along_line = l2.(Ref(sample_points[begin]), sample_points)
    return (distance_along_line, sample_points, sample_values)
end

