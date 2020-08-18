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
#   toc-showcode: true
# ---

# %%
using Pkg; Pkg.activate(@__DIR__)
using Test
using Simulation73
using NeuralModels
using DSP
using SpecialFunctions
#using Plots

function make_testing_lattice(; 
        n_points::NTuple{N,Int}=(1001,), extent::NTuple{N,T}=(300.0,), 
        type=CompactLattice{Float64,1}
            ) where {T,N}
    lattice = type(; n_points = n_points, extent = extent)
    return lattice
end
erf_sd(x,sd) = erf(x/(sqrt(2)*sd))
function analytical_gaussian_conv(xs::NTuple{1}, halves::NTuple{1}, σs::NTuple{1})
    x = xs[1]; half = halves[1]; σ = σs[1]
    (erf_sd(half - x, σ) + erf_sd(half + x, σ)) / (2)
end
function analytical_gaussian_conv(xs::NTuple{2}, halves::NTuple{2}, σs::NTuple{2})
    (x,y), (halfx, halfy), (σx, σy) = xs, halves, σs
    ((erf_sd(x-halfx, σx) - erf_sd(x+halfx, σx)) * (erf_sd(y-halfy, σy) - erf_sd(y+halfy, σy))) / 4
end
function make_manual_bump(lattice::AbstractLattice{T,N}, attempted_half_width::NTuple{N,T}) where {N,T}
    arr = zeros(T, size(lattice)...)
    mid_point = fft_center_dx(lattice)
    half_width_dx = CartesianIndex(floor.(Ref(Int), attempted_half_width ./ step(lattice)))
    start = CartesianIndex(Tuple(mid_point) .- Tuple(half_width_dx))
    stop = CartesianIndex(Tuple(mid_point) .+ Tuple(half_width_dx))
    arr[start:stop] .= 1.0
    half_width = Tuple(half_width_dx) .* step(lattice)
    return (arr, half_width)
end
    
        

# %% {"collapsed": true, "jupyter": {"outputs_hidden": true, "source_hidden": true}}
# Original; Bad
# Principal problem: Analytical solution uses imposed half-width
# Needs to use actual discretized half-width
Simulation73.origin_idx(lattice::AbstractLattice) = CartesianIndex(round.(Int, size(lattice) ./ 2, RoundNearestTiesUp))
function NeuralModels.kernel(conn::AbstractConnectivityParameter{T,N_CDT}, lattice::AbstractSpace{T,N_CDT}) where {T,N_CDT}
    directed_weights(conn, lattice, coordinates(lattice)[origin_idx(lattice)])
end

function fft_n_points_test(n_points, atol=0.01, rtol=0.1)
    (n_points,), (extent,), (dx,), mid_point, circle, circle_zeros = make_testing_lattice(n_points=(n_points,), type=PeriodicLattice{Float64,1})
    mid_point = Tuple(mid_point)[1]
    manual_bump = copy(circle_zeros)
    half_width = 30.0
    half_width_dx = floor(Int, half_width / dx)
    manual_bump[mid_point-half_width_dx:mid_point+half_width_dx] .= 1.0
    #half_width = half_width_dx * dx
    fft_gaussian_output, naive_gaussian_output, analytical_gaussian_output = nothing, nothing, nothing
    begin#@testset "Exponentially decaying connectivity" begin
        σ = 20.0
        gaussian_conn_param = GaussianConnectivityParameter(amplitude=1.0, spread=(σ,))
        gaussian_dsp_output = conv(NeuralModels.kernel(gaussian_conn_param, circle), manual_bump)
        abs_conn_param = ExpSumAbsDecayingConnectivityParameter(amplitude=1.0, spread=(σ,))
        abs_dsp_output = conv(NeuralModels.kernel(abs_conn_param, circle), manual_bump)
        #@test all(directed_weights(gaussian_conn_param, circle, (0.0,)) .≈ directed_weights(gaussian_conn_param, circle, (extent,)))
        begin#@testset "FFT" begin
            gaussian_conn_action = gaussian_conn_param(circle)
            abs_conn_action = abs_conn_param(circle)

            naive_gaussian_output = zeros(size(manual_bump)...)
            gaussian_conn_action(naive_gaussian_output, manual_bump, 0.0)

            fft_gaussian_conn_param = FFTParameter(gaussian_conn_param)
            fft_gaussian_conn_action = fft_gaussian_conn_param(circle)

            fft_gaussian_output = zeros(size(manual_bump)...)
            fft_gaussian_conn_action(fft_gaussian_output, manual_bump, 0.0)

            # Assumes analytical form of Gaussian normalizes to 1
            #analytical_gaussian_conv(x) =  (erf((half_width-x[1])/(σ * (√(2)))) + erf((x[1]+half_width)/(σ * (√(2))))) / (2)
            analytical_gaussian_output = analytical_gaussian_conv.(coordinates(circle), Ref((half_width,)), Ref((σ,)))

            #@test all(isapprox.(fft_gaussian_output,analytical_gaussian_output, atol=atol, rtol=rtol))
            #@test all(isapprox.(naive_gaussian_output,analytical_gaussian_output, atol=atol, rtol=rtol))
        end
        begin#@testset "FFT pops" begin
            sigmas = [20.0 10.0; 30.0 40.0] .|> (x) -> (x,)
            amplitudes = [1.0 -1.0; 1.0 1.0]
            gaussian_conn_pop_params = pops(GaussianConnectivityParameter;
                amplitude=amplitudes, spread=sigmas)


            pops_bump = population_repeat(manual_bump, 2)
            manual_input = copy(pops_bump)
            manual_output = zeros(size(pops_bump))
            fft_input = copy(pops_bump)
            fft_output = zeros(size(pops_bump))
            #plot(manual_input) |> display

            for i in 1:2
                for j in 1:2
                    manual_conn_action = FFTParameter(GaussianConnectivityParameter(amplitude=amplitudes[i,j], spread=sigmas[i,j]))(circle)
                    manual_conn_action(population(manual_output,i), population(manual_input,j), 0.0)
                end
            end

            fft_gaussian_conn_pop_params = FFTParameter(gaussian_conn_pop_params)
            fft_pops_action = fft_gaussian_conn_pop_params(circle)
            #plot(fft_input) |> display
            fft_pops_action(fft_output, fft_input, 0.0)

            #plot(fft_output) |> display
            #plot(manual_output) |> display
            #@test all(fft_output .== manual_output) # Manual output just unrolls what the pops should be doing; no approximation
        end
    end
    return (analytical_gaussian_output, naive_gaussian_output, fft_gaussian_output)
end

    
    
test_ns = 71:401
fft_diffs = zeros(length(test_ns)); naive_diffs = zeros(length(test_ns));
fft_max_diffs = zeros(length(test_ns)); naive_max_diffs = zeros(length(test_ns));


for (i,n) in enumerate(test_ns)
    (theory_sq_output, naive_sq_output, fft_sq_output) = fft_n_points_test(n)
    fft_Δ = abs.(theory_sq_output .- fft_sq_output)
    fft_diffs[i] = sum(fft_Δ) / n
    fft_max_diffs[i] = maximum(fft_Δ)

    naive_Δ = abs.(theory_sq_output .- naive_sq_output)
    naive_diffs[i] = sum(naive_Δ) / n
    naive_max_diffs[i] = maximum(naive_Δ)
end

scatter(test_ns, [fft_diffs, naive_diffs], labels=[:fft, :naive], title="mean diffs",
    xlabel="# bins", ylabel="sum(abs(theory - numeric)) / # bins")

scatter(test_ns, [fft_max_diffs, naive_max_diffs], labels=[:fft, :naive], title="max diffs",
    xlabel="# bins", ylabel="max(abs(theory - numeric))")

# %% {"jupyter": {"source_hidden": true}}
# Good analytical from discretized stimulus; bad kernel
# The kernel doesn't turn out to matter that much.

Simulation73.origin_idx(lattice::AbstractLattice) = CartesianIndex(round.(Int, size(lattice) ./ 2, RoundNearestTiesUp))
function NeuralModels.kernel(conn::AbstractConnectivityParameter{T,N_CDT}, lattice::AbstractSpace{T,N_CDT}) where {T,N_CDT}
    directed_weights(conn, lattice, coordinates(lattice)[origin_idx(lattice)])
end

function fft_n_points_test(n_points, atol=0.01, rtol=0.1)
    n_points, extent, dx, mid_point, circle, circle_zeros = make_testing_lattice(n_points=n_points, type=PeriodicLattice{Float64,1})
    manual_bump = copy(circle_zeros)
    attempted_half_width = 30.0
    half_width_dx = floor(Int, attempted_half_width / dx)
    manual_bump[mid_point-half_width_dx:mid_point+half_width_dx] .= 1.0
    half_width = half_width_dx * dx
    fft_gaussian_output, naive_gaussian_output, analytical_gaussian_output = nothing, nothing, nothing
    begin#@testset "Exponentially decaying connectivity" begin
        σ = 20.0
        gaussian_conn_param = GaussianConnectivityParameter(amplitude=1.0, spread=(σ,))
        gaussian_dsp_output = conv(NeuralModels.kernel(gaussian_conn_param, circle), manual_bump)
        abs_conn_param = ExpSumAbsDecayingConnectivityParameter(amplitude=1.0, spread=(σ,))
        abs_dsp_output = conv(NeuralModels.kernel(abs_conn_param, circle), manual_bump)
        #@test all(directed_weights(gaussian_conn_param, circle, (0.0,)) .≈ directed_weights(gaussian_conn_param, circle, (extent,)))
        begin#@testset "FFT" begin
            gaussian_conn_action = gaussian_conn_param(circle)
            abs_conn_action = abs_conn_param(circle)

            naive_gaussian_output = zeros(size(manual_bump)...)
            gaussian_conn_action(naive_gaussian_output, manual_bump, 0.0)

            fft_gaussian_conn_param = FFTParameter(gaussian_conn_param)
            fft_gaussian_conn_action = fft_gaussian_conn_param(circle)

            fft_gaussian_output = zeros(size(manual_bump)...)
            fft_gaussian_conn_action(fft_gaussian_output, manual_bump, 0.0)

            # Assumes analytical form of Gaussian normalizes to 1
            analytical_gaussian_conv(x) =  (erf((half_width-x[1])/σ) + erf((x[1]+half_width)/σ)) / (2)
            analytical_gaussian_output = analytical_gaussian_conv.(coordinates(circle))

            #@test all(isapprox.(fft_gaussian_output,analytical_gaussian_output, atol=atol, rtol=rtol))
            #@test all(isapprox.(naive_gaussian_output,analytical_gaussian_output, atol=atol, rtol=rtol))
        end
        begin#@testset "FFT pops" begin
            sigmas = [20.0 10.0; 30.0 40.0] .|> (x) -> (x,)
            amplitudes = [1.0 -1.0; 1.0 1.0]
            gaussian_conn_pop_params = pops(GaussianConnectivityParameter;
                amplitude=amplitudes, spread=sigmas)


            pops_bump = population_repeat(manual_bump, 2)
            manual_input = copy(pops_bump)
            manual_output = zeros(size(pops_bump))
            fft_input = copy(pops_bump)
            fft_output = zeros(size(pops_bump))
            #plot(manual_input) |> display

            for i in 1:2
                for j in 1:2
                    manual_conn_action = FFTParameter(GaussianConnectivityParameter(amplitude=amplitudes[i,j], spread=sigmas[i,j]))(circle)
                    manual_conn_action(population(manual_output,i), population(manual_input,j), 0.0)
                end
            end

            fft_gaussian_conn_pop_params = FFTParameter(gaussian_conn_pop_params)
            fft_pops_action = fft_gaussian_conn_pop_params(circle)
            #plot(fft_input) |> display
            fft_pops_action(fft_output, fft_input, 0.0)

            #plot(fft_output) |> display
            #plot(manual_output) |> display
            #@test all(fft_output .== manual_output) # Manual output just unrolls what the pops should be doing; no approximation
        end
    end
    return (analytical_gaussian_output, naive_gaussian_output, fft_gaussian_output)
end
    
    
test_ns = 71:401
fft_diffs = zeros(length(test_ns)); naive_diffs = zeros(length(test_ns));
fft_max_diffs = zeros(length(test_ns)); naive_max_diffs = zeros(length(test_ns));


for (i,n) in enumerate(test_ns)
    (theory_sq_output, naive_sq_output, fft_sq_output) = fft_n_points_test(n)
    fft_Δ = abs.(theory_sq_output .- fft_sq_output)
    fft_diffs[i] = sum(fft_Δ) / n
    fft_max_diffs[i] = maximum(fft_Δ)

    naive_Δ = abs.(theory_sq_output .- naive_sq_output)
    naive_diffs[i] = sum(naive_Δ) / n
    naive_max_diffs[i] = maximum(naive_Δ)
end

scatter(test_ns, [fft_diffs, naive_diffs], labels=[:fft, :naive], title="mean diffs",
    xlabel="# bins", ylabel="sum(abs(theory - numeric)) / # bins")

scatter(test_ns, [fft_max_diffs, naive_max_diffs], labels=[:fft, :naive], title="max diffs",
    xlabel="# bins", ylabel="max(abs(theory - numeric))")
   

# %% {"jupyter": {"source_hidden": true}}
# Good theoretical solution; bad origin_idx; new BAD kernel
# Kernel goes from dist=0.0, but not centered correctly

Simulation73.origin_idx(lattice::AbstractLattice) = CartesianIndex(round.(Int, size(lattice) ./ 2, RoundNearestTiesUp))
function NeuralModels.kernel(conn::AbstractConnectivityParameter{T,N_CDT}, lattice::AbstractSpace{T,N_CDT}) where {T,N_CDT}
    directed_weights(conn, lattice, Tuple(zeros(T,N_CDT)))
end

function fft_n_points_test(n_points, atol=0.01, rtol=0.1)
    n_points, extent, dx, mid_point, circle, circle_zeros = make_testing_lattice(n_points=n_points, type=PeriodicLattice{Float64,1})
    manual_bump = copy(circle_zeros)
    attempted_half_width = 30.0
    half_width_dx = floor(Int, attempted_half_width / dx)
    manual_bump[mid_point-half_width_dx:mid_point+half_width_dx] .= 1.0
    half_width = half_width_dx * dx
    fft_gaussian_output, naive_gaussian_output, analytical_gaussian_output = nothing, nothing, nothing
    begin#@testset "Exponentially decaying connectivity" begin
        σ = 20.0
        gaussian_conn_param = GaussianConnectivityParameter(amplitude=1.0, spread=(σ,))
        gaussian_dsp_output = conv(NeuralModels.kernel(gaussian_conn_param, circle), manual_bump)
        abs_conn_param = ExpSumAbsDecayingConnectivityParameter(amplitude=1.0, spread=(σ,))
        abs_dsp_output = conv(NeuralModels.kernel(abs_conn_param, circle), manual_bump)
        #@test all(directed_weights(gaussian_conn_param, circle, (0.0,)) .≈ directed_weights(gaussian_conn_param, circle, (extent,)))
        begin#@testset "FFT" begin
            gaussian_conn_action = gaussian_conn_param(circle)
            abs_conn_action = abs_conn_param(circle)

            naive_gaussian_output = zeros(size(manual_bump)...)
            gaussian_conn_action(naive_gaussian_output, manual_bump, 0.0)

            fft_gaussian_conn_param = FFTParameter(gaussian_conn_param)
            fft_gaussian_conn_action = fft_gaussian_conn_param(circle)

            fft_gaussian_output = zeros(size(manual_bump)...)
            fft_gaussian_conn_action(fft_gaussian_output, manual_bump, 0.0)

            # Assumes analytical form of Gaussian normalizes to 1
            analytical_gaussian_conv(x) =  (erf((half_width-x[1])/σ) + erf((x[1]+half_width)/σ)) / (2)
            analytical_gaussian_output = analytical_gaussian_conv.(coordinates(circle))

            #@test all(isapprox.(fft_gaussian_output,analytical_gaussian_output, atol=atol, rtol=rtol))
            #@test all(isapprox.(naive_gaussian_output,analytical_gaussian_output, atol=atol, rtol=rtol))
        end
        begin#@testset "FFT pops" begin
            sigmas = [20.0 10.0; 30.0 40.0] .|> (x) -> (x,)
            amplitudes = [1.0 -1.0; 1.0 1.0]
            gaussian_conn_pop_params = pops(GaussianConnectivityParameter;
                amplitude=amplitudes, spread=sigmas)


            pops_bump = population_repeat(manual_bump, 2)
            manual_input = copy(pops_bump)
            manual_output = zeros(size(pops_bump))
            fft_input = copy(pops_bump)
            fft_output = zeros(size(pops_bump))
            #plot(manual_input) |> display

            for i in 1:2
                for j in 1:2
                    manual_conn_action = FFTParameter(GaussianConnectivityParameter(amplitude=amplitudes[i,j], spread=sigmas[i,j]))(circle)
                    manual_conn_action(population(manual_output,i), population(manual_input,j), 0.0)
                end
            end

            fft_gaussian_conn_pop_params = FFTParameter(gaussian_conn_pop_params)
            fft_pops_action = fft_gaussian_conn_pop_params(circle)
            #plot(fft_input) |> display
            fft_pops_action(fft_output, fft_input, 0.0)

            #plot(fft_output) |> display
            #plot(manual_output) |> display
            #@test all(fft_output .== manual_output) # Manual output just unrolls what the pops should be doing; no approximation
        end
    end
    return (analytical_gaussian_output, naive_gaussian_output, fft_gaussian_output)
end
    
    
test_ns = 71:401
fft_diffs = zeros(length(test_ns)); naive_diffs = zeros(length(test_ns));
fft_max_diffs = zeros(length(test_ns)); naive_max_diffs = zeros(length(test_ns));


for (i,n) in enumerate(test_ns)
    (theory_sq_output, naive_sq_output, fft_sq_output) = fft_n_points_test(n)
    fft_Δ = abs.(theory_sq_output .- fft_sq_output)
    fft_diffs[i] = sum(fft_Δ) / n
    fft_max_diffs[i] = maximum(fft_Δ)

    naive_Δ = abs.(theory_sq_output .- naive_sq_output)
    naive_diffs[i] = sum(naive_Δ) / n
    naive_max_diffs[i] = maximum(naive_Δ)
end

scatter(test_ns, [fft_diffs, naive_diffs], labels=[:fft, :naive], title="mean diffs",
    xlabel="# bins", ylabel="sum(abs(theory - numeric)) / # bins")

scatter(test_ns, [fft_max_diffs, naive_max_diffs], labels=[:fft, :naive], title="max diffs",
    xlabel="# bins", ylabel="max(abs(theory - numeric))")
   

# %% {"collapsed": true, "jupyter": {"outputs_hidden": true, "source_hidden": true}}
# Changed origin_idx AND kernel AND theory from discretization
# Now the kernel is centered correctly
# No improvement from previous
function NeuralModels.kernel(conn::AbstractConnectivityParameter{T,N_CDT}, lattice::AbstractSpace{T,N_CDT}) where {T,N_CDT}
    # Kernel has ZERO DIST at its center (or floor(extent/2) + 1)
    fft_centered_differences = differences(lattice, coordinates(lattice)[fft_center_dx(lattice)])
    NeuralModels.apply_connectivity(conn, fft_centered_differences, step(lattice), fft_centered_differences)    
end

function fft_n_points_test(n_points, atol=0.01, rtol=0.1)
    n_points, extent, dx, mid_point, circle, circle_zeros = make_testing_lattice(n_points=n_points, type=PeriodicLattice{Float64,1})
    manual_bump = copy(circle_zeros)
    attempted_half_width = 30.0
    half_width_dx = floor(Int, attempted_half_width / dx)
    manual_bump[mid_point-half_width_dx:mid_point+half_width_dx] .= 1.0
    half_width = half_width_dx * dx
    fft_gaussian_output, naive_gaussian_output, analytical_gaussian_output = nothing, nothing, nothing
    begin#@testset "Exponentially decaying connectivity" begin
        σ = 20.0
        gaussian_conn_param = GaussianConnectivityParameter(amplitude=1.0, spread=(σ,))
        gaussian_dsp_output = conv(NeuralModels.kernel(gaussian_conn_param, circle), manual_bump)
        abs_conn_param = ExpSumAbsDecayingConnectivityParameter(amplitude=1.0, spread=(σ,))
        abs_dsp_output = conv(NeuralModels.kernel(abs_conn_param, circle), manual_bump)
        #@test all(directed_weights(gaussian_conn_param, circle, (0.0,)) .≈ directed_weights(gaussian_conn_param, circle, (extent,)))
        begin#@testset "FFT" begin
            gaussian_conn_action = gaussian_conn_param(circle)
            abs_conn_action = abs_conn_param(circle)

            naive_gaussian_output = zeros(size(manual_bump)...)
            gaussian_conn_action(naive_gaussian_output, manual_bump, 0.0)

            fft_gaussian_conn_param = FFTParameter(gaussian_conn_param)
            fft_gaussian_conn_action = fft_gaussian_conn_param(circle)

            fft_gaussian_output = zeros(size(manual_bump)...)
            fft_gaussian_conn_action(fft_gaussian_output, manual_bump, 0.0)

            # Assumes analytical form of Gaussian normalizes to 1
            analytical_gaussian_conv(x) =  (erf((half_width-x[1])/σ) + erf((x[1]+half_width)/σ)) / (2)
            analytical_gaussian_output = analytical_gaussian_conv.(coordinates(circle))

            #@test all(isapprox.(fft_gaussian_output,analytical_gaussian_output, atol=atol, rtol=rtol))
            #@test all(isapprox.(naive_gaussian_output,analytical_gaussian_output, atol=atol, rtol=rtol))
        end
        begin#@testset "FFT pops" begin
            sigmas = [20.0 10.0; 30.0 40.0] .|> (x) -> (x,)
            amplitudes = [1.0 -1.0; 1.0 1.0]
            gaussian_conn_pop_params = pops(GaussianConnectivityParameter;
                amplitude=amplitudes, spread=sigmas)


            pops_bump = population_repeat(manual_bump, 2)
            manual_input = copy(pops_bump)
            manual_output = zeros(size(pops_bump))
            fft_input = copy(pops_bump)
            fft_output = zeros(size(pops_bump))
            #plot(manual_input) |> display

            for i in 1:2
                for j in 1:2
                    manual_conn_action = FFTParameter(GaussianConnectivityParameter(amplitude=amplitudes[i,j], spread=sigmas[i,j]))(circle)
                    manual_conn_action(population(manual_output,i), population(manual_input,j), 0.0)
                end
            end

            fft_gaussian_conn_pop_params = FFTParameter(gaussian_conn_pop_params)
            fft_pops_action = fft_gaussian_conn_pop_params(circle)
            #plot(fft_input) |> display
            fft_pops_action(fft_output, fft_input, 0.0)

            #plot(fft_output) |> display
            #plot(manual_output) |> display
            #@test all(fft_output .== manual_output) # Manual output just unrolls what the pops should be doing; no approximation
        end
    end
    return (analytical_gaussian_output, naive_gaussian_output, fft_gaussian_output)
end
    
    
test_ns = 71:401
fft_diffs = zeros(length(test_ns)); naive_diffs = zeros(length(test_ns));
fft_max_diffs = zeros(length(test_ns)); naive_max_diffs = zeros(length(test_ns));
using LinearAlgebra

for (i,n) in enumerate(test_ns)
    (theory_sq_output, naive_sq_output, fft_sq_output) = fft_n_points_test(n)
    fft_Δ = abs.(theory_sq_output .- fft_sq_output)
    fft_diffs[i] = norm(fft_Δ) / norm(theory_sq_output)
    fft_max_diffs[i] = maximum(fft_Δ ./ theory_sq_output)

    naive_Δ = abs.(theory_sq_output .- naive_sq_output)
    naive_diffs[i] = norm(naive_Δ) / norm(theory_sq_output)
    naive_max_diffs[i] = maximum(naive_Δ ./ theory_sq_output)
end

scatter(test_ns, [fft_diffs, naive_diffs], labels=[:fft, :naive], title="mean abs diffs",
    xlabel="# bins", ylabel="norm(theory - numeric) / norm(theory)") |> display

# scatter(test_ns, [fft_max_diffs, naive_max_diffs], labels=[:fft, :naive], title="max diffs(odd)",
#     xlabel="# bins", ylabel="max(abs(theory - numeric) / theory)") |> display
    

# %% [markdown] {"jupyter": {"source_hidden": true}}
# ## 

# %%
# Here: changing to allow multidimensional test
# Changed origin_idx AND kernel AND theory from discretization
# Disabled testing of conventional (naive) method because too slow

fft_n_points_test(n::Int; kwargs...) = fft_n_points_test((n,); kwargs...)
function fft_n_points_test_no_naive(n_points::NTuple{N,Int}, atol=0.01, rtol=0.1) where {N}
    circle = make_testing_lattice(extent=Tuple(zeros(N) .+ 300.0), n_points=n_points, type=PeriodicLattice{Float64,N})
    attempted_half_width = Tuple(zeros(N) .+ 30.0)
    manual_bump, half_width = make_manual_bump(circle, attempted_half_width)
    fft_gaussian_output, naive_gaussian_output, analytical_gaussian_output = nothing, nothing, nothing
    begin#@testset "Exponentially decaying connectivity" begin
        σ = Tuple(zeros(N) .+ 20.0)
        gaussian_conn_param = GaussianConnectivityParameter(amplitude=1.0, spread=σ)
        begin#@testset "FFT" begin
            gaussian_conn_action = gaussian_conn_param(circle)

            fft_gaussian_conn_param = FFTParameter(gaussian_conn_param)
            fft_gaussian_conn_action = fft_gaussian_conn_param(circle)

            fft_gaussian_output = zeros(size(manual_bump)...)
            fft_gaussian_conn_action(fft_gaussian_output, manual_bump, 0.0)

            # Assumes analytical form of Gaussian normalizes to 
            analytical_gaussian_output = analytical_gaussian_conv.(coordinates(circle), Ref(half_width), Ref(σ))
        end
    end
    return (analytical_gaussian_output, fft_gaussian_output)
end
    
    
test_ns = 901:201:3001
fft_diffs = zeros(length(test_ns)); naive_diffs = zeros(length(test_ns));
fft_max_diffs = zeros(length(test_ns)); naive_max_diffs = zeros(length(test_ns));
using LinearAlgebra

for (i,n) in enumerate(test_ns)
    #(theory_sq_output, naive_sq_output, fft_sq_output) = fft_n_points_test((n,n))
    (theory_sq_output, fft_sq_output) = fft_n_points_test_no_naive((n,n))
    fft_Δ = abs.(theory_sq_output .- fft_sq_output)
    fft_diffs[i] = norm(fft_Δ) / norm(theory_sq_output)
    fft_max_diffs[i] = maximum(fft_Δ ./ theory_sq_output)

#     naive_Δ = abs.(theory_sq_output .- naive_sq_output)
#     naive_diffs[i] = norm(naive_Δ) / norm(theory_sq_output)
#     naive_max_diffs[i] = maximum(naive_Δ ./ theory_sq_output)
end

scatter(test_ns, fft_diffs, labels=:fft, title="norm diffs (2D)",
    xlabel="# bins", ylabel="norm(theory - numeric) / norm(theory)") |> display

# scatter(test_ns, [fft_max_diffs, naive_max_diffs], labels=[:fft, :naive], title="max diffs(odd)",
#     xlabel="# bins", ylabel="max(abs(theory - numeric) / theory)") |> display
    

# %%
# Flip center of even
# Here: changing to allow multidimensional test
# Changed origin_idx AND kernel AND theory from discretization
# Disabled testing of conventional (naive) method because too slow

fft_n_points_test(n::Int; kwargs...) = fft_n_points_test((n,); kwargs...)
function fft_n_points_test_no_naive(n_points::NTuple{N,Int}, atol=0.01, rtol=0.1) where {N}
    circle = make_testing_lattice(extent=Tuple(zeros(N) .+ 300.0), n_points=n_points, type=PeriodicLattice{Float64,N})
    attempted_half_width = Tuple(zeros(N) .+ 30.0)
    manual_bump, half_width = make_manual_bump(circle, attempted_half_width)
    fft_gaussian_output, naive_gaussian_output, analytical_gaussian_output = nothing, nothing, nothing
    begin#@testset "Exponentially decaying connectivity" begin
        σ = Tuple(zeros(N) .+ 20.0)
        gaussian_conn_param = GaussianConnectivityParameter(amplitude=1.0, spread=σ)
        begin#@testset "FFT" begin
            gaussian_conn_action = gaussian_conn_param(circle)

            fft_gaussian_conn_param = FFTParameter(gaussian_conn_param)
            fft_gaussian_conn_action = fft_gaussian_conn_param(circle)

            fft_gaussian_output = zeros(size(manual_bump)...)
            fft_gaussian_conn_action(fft_gaussian_output, manual_bump, 0.0)

            # Assumes analytical form of Gaussian normalizes to 
            analytical_gaussian_output = analytical_gaussian_conv.(coordinates(circle), Ref(half_width), Ref(σ))
        end
    end
    return (analytical_gaussian_output, fft_gaussian_output)
end
    
    
test_ns = 201:3:501
fft_diffs = zeros(length(test_ns)); naive_diffs = zeros(length(test_ns));
fft_max_diffs = zeros(length(test_ns)); naive_max_diffs = zeros(length(test_ns));
using LinearAlgebra

for (i,n) in enumerate(test_ns)
    (theory_sq_output, fft_sq_output) = fft_n_points_test_no_naive((n,n))
    
    fft_Δ = (theory_sq_output .- fft_sq_output) .^ 2
    fft_diffs[i] = sum(fft_Δ) ./ (n^2)
    fft_max_diffs[i] = maximum(fft_Δ ./ theory_sq_output)

#     naive_Δ = abs.(theory_sq_output .- naive_sq_output)
#     naive_diffs[i] = norm(naive_Δ) / norm(theory_sq_output)
#     naive_max_diffs[i] = maximum(naive_Δ ./ theory_sq_output)
end

scatter(test_ns, fft_diffs, labels=:fft, title="mean sq error (2D)",
    xlabel="sqrt(# bins)", ylabel="mse(theory, numeric)") |> display

# scatter(test_ns, fft_diffs, labels=:fft, title="norm diffs (2D)",
#     xlabel="sqrt(# bins)", ylabel="max(sqerr / theory)") |> display


# scatter(test_ns, [fft_max_diffs, naive_max_diffs], labels=[:fft, :naive], title="max diffs(odd)",
#     xlabel="# bins", ylabel="max(abs(theory - numeric) / theory)") |> display
    

# %%
savefig("mse2d.png")

# %%
n_points = (901,901)
(theory_sq_output, fft_sq_output) = fft_n_points_test_no_naive(n_points)
normed_theory = theory_sq_output ./ norm(theory_sq_output)
normed_fft = fft_sq_output ./ norm(fft_sq_output)
plot(heatmap.([theory_sq_output, fft_sq_output])..., title="Theory and FFT output: $n_points") |> display
savefig("example901.png")
heatmap(theory_sq_output - fft_sq_output, title="Difference between Theory and FFT output: $n_points") |> display
savefig("diff901.png")
n_points = (900,900)
(theory_sq_output, fft_sq_output) = fft_n_points_test_no_naive(n_points)
heatmap(theory_sq_output - fft_sq_output, title="Difference between Theory and FFT output: $n_points") |> display
savefig("diff900.png")
n_points = (1028,1028)
(theory_sq_output, fft_sq_output) = fft_n_points_test_no_naive(n_points)
heatmap(theory_sq_output - fft_sq_output, title="Difference between Theory and FFT output: $n_points") |> display
savefig("diff1028.png")

# %%
(theory_sq_output, fft_sq_output) = fft_n_points_test_no_naive((100,))
normed_theory = theory_sq_output ./ norm(theory_sq_output)
normed_fft = fft_sq_output ./ norm(fft_sq_output)
plot([theory_sq_output, fft_sq_output]) |> display
plot([normed_theory, normed_fft]) |> display
@show norm(theory_sq_output) / norm(fft_sq_output)

# %%
(2pi)^0.5

# %% {"collapsed": true, "jupyter": {"outputs_hidden": true, "source_hidden": true}}
these_n_points = (81,81)
theory_sq_output, fft_sq_output = fft_n_points_test_no_naive(these_n_points)
n_points, extent, dx, mid_point, circle, circle_zeros = make_testing_lattice(extent=(300.0,300.0), n_points=these_n_points, type=PeriodicLattice{Float64,2})
plot([heatmap(theory_sq_output), heatmap(fft_sq_output)]...) |> display
@show norm(theory_sq_output)
@show norm(fft_sq_output)
@show norm(theory_sq_output .- fft_sq_output)
@show(step(circle) |> prod)

# %% [markdown] {"jupyter": {"source_hidden": true}}
# Wim suggests making sure the edges are actually zero.

# %% {"collapsed": true, "jupyter": {"outputs_hidden": true, "source_hidden": true}}
n_points, extent, dx, mid_point, circle, circle_zeros = make_testing_lattice(n_points=(71,), type=PeriodicLattice{Float64,1})
gaussian_conn_param = GaussianConnectivityParameter(amplitude=1.0, spread=(20.0,))
kern = NeuralModels.kernel(gaussian_conn_param, circle)
@show kern[1]
@show kern[35]
@show kern[end]
@show kern[1] < 1e-15


# %% {"collapsed": true, "jupyter": {"outputs_hidden": true, "source_hidden": true}}
# ZEROING EDGES
# Changed origin_idx AND kernel AND theory from discretization
# Now the kernel is centered correctly
# No improvement from previous
function NeuralModels.kernel(conn::AbstractConnectivityParameter{T,N_CDT}, lattice::AbstractSpace{T,N_CDT}) where {T,N_CDT}
    # Kernel has ZERO DIST at its center (or floor(extent/2) + 1)
    fft_centered_differences = differences(lattice, coordinates(lattice)[fft_center_dx(lattice)])
    kern = NeuralModels.apply_connectivity(conn, fft_centered_differences, step(lattice), fft_centered_differences)
    kern[abs.(kern) .< 1e-15] .= 0.0
    return kern
end
function fft_n_points_test(n_points::NTuple{N,Int}, atol=0.01, rtol=0.1) where {N}
    n_points, extent, dx, mid_point, circle, circle_zeros = make_testing_lattice(n_points=n_points, extent=Tuple(zeros(N) .+ 700.0), type=PeriodicLattice{Float64,1})
    manual_bump = copy(circle_zeros)
    attempted_half_width = (30.0,)
    half_width_dx = CartesianIndex(floor.(Ref(Int), attempted_half_width ./ dx))
    start = CartesianIndex(Tuple(mid_point) .- Tuple(half_width_dx))
    stop = CartesianIndex(Tuple(mid_point) .+ Tuple(half_width_dx))
    manual_bump[start:stop] .= 1.0
    half_width = Tuple(half_width_dx) .* dx
    fft_gaussian_output, naive_gaussian_output, analytical_gaussian_output = nothing, nothing, nothing
    begin#@testset "Exponentially decaying connectivity" begin
        σ = 65.0
        gaussian_conn_param = GaussianConnectivityParameter(amplitude=1.0, spread=(σ,))       
        @test NeuralModels.kernel(gaussian_conn_param, circle)[1] == 0.0
        gaussian_dsp_output = conv(NeuralModels.kernel(gaussian_conn_param, circle), manual_bump)
#         abs_conn_param = ExpSumAbsDecayingConnectivityParameter(amplitude=1.0, spread=(σ,))
#         abs_dsp_output = conv(NeuralModels.kernel(abs_conn_param, circle), manual_bump)
        #@test all(directed_weights(gaussian_conn_param, circle, (0.0,)) .≈ directed_weights(gaussian_conn_param, circle, (extent,)))
        begin#@testset "FFT" begin
            gaussian_conn_action = gaussian_conn_param(circle)
#             abs_conn_action = abs_conn_param(circle)

            naive_gaussian_output = zeros(size(manual_bump)...)
            gaussian_conn_action(naive_gaussian_output, manual_bump, 0.0)

            fft_gaussian_conn_param = FFTParameter(gaussian_conn_param)
            fft_gaussian_conn_action = fft_gaussian_conn_param(circle)

            fft_gaussian_output = zeros(size(manual_bump)...)
            fft_gaussian_conn_action(fft_gaussian_output, manual_bump, 0.0)

            # Assumes analytical form of Gaussian normalizes to 1
            analytical_gaussian_conv(x) =  sum(@. (erf((half_width-x)/σ) + erf((x+half_width)/σ)) / (2))
            analytical_gaussian_output = analytical_gaussian_conv.(coordinates(circle))

            #@test all(isapprox.(fft_gaussian_output,analytical_gaussian_output, atol=atol, rtol=rtol))
            #@test all(isapprox.(naive_gaussian_output,analytical_gaussian_output, atol=atol, rtol=rtol))
        end
        begin#@testset "FFT pops" begin
            sigmas = [20.0 10.0; 30.0 40.0] .|> (x) -> (x,)
            amplitudes = [1.0 -1.0; 1.0 1.0]
            gaussian_conn_pop_params = pops(GaussianConnectivityParameter;
                amplitude=amplitudes, spread=sigmas)


            pops_bump = population_repeat(manual_bump, 2)
            manual_input = copy(pops_bump)
            manual_output = zeros(size(pops_bump))
            fft_input = copy(pops_bump)
            fft_output = zeros(size(pops_bump))
            #plot(manual_input) |> display

            for i in 1:2
                for j in 1:2
                    manual_conn_action = FFTParameter(GaussianConnectivityParameter(amplitude=amplitudes[i,j], spread=sigmas[i,j]))(circle)
                    manual_conn_action(population(manual_output,i), population(manual_input,j), 0.0)
                end
            end

            fft_gaussian_conn_pop_params = FFTParameter(gaussian_conn_pop_params)
            fft_pops_action = fft_gaussian_conn_pop_params(circle)
            #plot(fft_input) |> display
            fft_pops_action(fft_output, fft_input, 0.0)

            #plot(fft_output) |> display
            #plot(manual_output) |> display
            #@test all(fft_output .== manual_output) # Manual output just unrolls what the pops should be doing; no approximation
        end
    end
    return (analytical_gaussian_output, naive_gaussian_output, fft_gaussian_output)
end
    
    
test_ns = 100:2:110
fft_diffs = zeros(length(test_ns)); naive_diffs = zeros(length(test_ns));
fft_max_diffs = zeros(length(test_ns)); naive_max_diffs = zeros(length(test_ns));
using LinearAlgebra

for (i,n) in enumerate(test_ns)
    (theory_sq_output, naive_sq_output, fft_sq_output) = fft_n_points_test((n,))
    fft_Δ = (theory_sq_output .- fft_sq_output)
    fft_diffs[i] = norm(fft_Δ) / norm(theory_sq_output)
    fft_max_diffs[i] = maximum(fft_Δ ./ theory_sq_output)

    naive_Δ = (theory_sq_output .- naive_sq_output)
    naive_diffs[i] = norm(naive_Δ) / norm(theory_sq_output)
    naive_max_diffs[i] = maximum(naive_Δ ./ theory_sq_output)
end

scatter(test_ns, [fft_diffs, naive_diffs], labels=[:fft, :naive], title="normalized norm diff",
    xlabel="# bins", ylabel="norm(theory - numeric) / norm(theory)") |> display

# scatter(test_ns, [fft_max_diffs, naive_max_diffs], labels=[:fft, :naive], title="max diffs(odd)",
#     xlabel="# bins", ylabel="max(abs(theory - numeric) / theory)") |> display
    

# %%
savefig("fft_three.png")

# %% {"jupyter": {"source_hidden": true}}
(theory_sq_output, naive_sq_output, fft_sq_output) = fft_n_points_test(70)
plot([theory_sq_output, naive_sq_output .+ 0.002, fft_sq_output], labels=[:theory, :naive, :fft], title = "conv output for n=70",
    xlabel = "space", ylabel="activity")
savefig("conv_output_70.png")

# %% {"collapsed": true, "jupyter": {"outputs_hidden": true, "source_hidden": true}}
function Simulation73.fft_center_dx(space::AbstractSpace)
    CartesianIndex(floor.(Ref(Int), size(space) ./ 2) .+ 1)#(size(space) .% 2))
end
(theory_sq_output, naive_sq_output, fft_sq_output) = fft_n_points_test(71)
plot([theory_sq_output, naive_sq_output, fft_sq_output], labels=[:theory, :naive, :fft], title = "conv output for n=71",
    xlabel = "space", ylabel="activity") |> display
savefig("conv_output_71.png")
