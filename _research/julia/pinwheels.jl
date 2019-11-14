# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all
#     notebook_metadata_filter: all,-language_info
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Julia 1.3.0-rc4
#     language: julia
#     name: julia-1.3
# ---

# %%
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
    return gaussian_circle
end

function make_pinwheels(circle_arr::AbstractArray{T,2}) where T
    return mod2pi.(angle.(ifft(circle_arr)) .* 2) ./ 2
end

function make_pinwheels(extent::Union{Tuple,Int}, circle_radius = 0.2, circle_width = 0.02)
    gaussian_circle = make_gaussian_circle(extent, circle_radius, circle_width)
    make_pinwheels(gaussian_circle)
end

# %%
pinwheels = make_pinwheels((512, 512), 0.015, 0.001)
heatmap(pinwheels, color=:colorwheel, clim=(0,π), xlab="x-space", ylab="y-space")

# %%
savefig("pinwheels_015.png")
