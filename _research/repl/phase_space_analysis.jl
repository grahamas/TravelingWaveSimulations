
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

function move_left_boundary!(queue, new_left)
    removed_val = 0
    while !isempty(queue) && first(queue)[1] < new_left
        removed_elt = dequeue!(queue)
        removed_val += removed_elt[2]
    end
    return removed_val
end

function move_right_boundary!(queue, new_right, elements, next_elt)
    added_val = 0
    while next_elt <= length(elements) && elements[next_elt][1] <= new_right
        new_elt = elements[next_elt]
        enqueue!(queue, new_elt)
        added_val += new_elt[2]
        next_elt += 1
    end
    return (added_val, next_elt)
end

function moving_average(data_pairs::Vector{TT}, window_proportion=0.05, step_proportion=0.01,
                        extent::T1=data_pairs[end][1], 
                        window::T1=window_proportion*extent, 
                        step::T1=step_proportion*extent) where {T1 <: Number, T2 <: Number, TT<:Tuple{T1,T2}} 
    start = data_pairs[begin][1]
    stop = data_pairs[end][1]
    locs = start+(window/2):step:stop-(window/2)
    @assert length(locs) > 10
    next_pair = 1
    cur_sum = 0
    pairs_in_window = Queue{TT}()
    vals = map(locs) do loc
        left, right = loc-window/2, loc+window/2 
        cur_sum -= move_left_boundary!(pairs_in_window, left)
        inc, next_pair = move_right_boundary!(pairs_in_window, right, data_pairs, next_pair)
        cur_sum += inc
        cur_sum / length(pairs_in_window)
    end
    return locs, vals
end

function project(coord::C, line::C, line_point::C) where {N, C <: SVector{N}}
    projection = (line * line') / (line' * line)
    return projection * coord + (I - projection) * line_point
end

squish(data::NamedDimsArray, args...) = squish(data.data, args...)
function squish(data::AbstractAxisArray, target_line, target_line_origin)
    coordinates = SVector.(get_coordinates(data))
    # project coordinates onto line and convert into distance along line from y-intercept
    projected_coords = project.(coordinates, Ref(target_line), Ref(target_line_origin))
    projected_data_pairs = [(l2(coord, target_line_origin), datum) for (coord, datum) in zip(projected_coords, data)][:]
    sort!(projected_data_pairs)
    squished_line_dists, squished_line_vals = moving_average(projected_data_pairs)
    return (squished_line_dists, squished_line_vals, projected_coords)
end
