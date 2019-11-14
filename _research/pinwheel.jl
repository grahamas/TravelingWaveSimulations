using FFTW, Plots

function is_circle(n_points, r, width=0.01)
    coords = Iterators.product([LinRange(-1,1,n) for n in n_points]...) |> collect
    return map(coords) do coord
        (r-width) <= sqrt(sum(coord .^ 2)) <= (r + width)
    end
end

function make_gaussian_circle(extent, circle_radius, circle_width)
    is_circ = is_circle(extent, circle_radius, circle_width)
    gaussian_circle = zeros(size(is_circ)...)
    gaussian_circle[is_circ] .= randn(sum(is_circ))
end

function make_pinwheels(circle_arr::AbstractArray{T,2}) where T
    return mod2pi.(angle.(ifft(circle_arr)) .* 2) ./ 2
end

function make_pinwheels(extent, circle_radius = 0.2, circle_width = 0.02)
    gaussian_circle = make_gaussian_circle(extent, circle_radius, circle_width)
    make_pinwheels(gaussian_circle)
end

circle = make_gaussian_circle()
pinwheels = make_pinwheels((512,512), 0.01, 0.001)
heatmap(pinwheels, color=:colorwheel, clim=(0,π))
