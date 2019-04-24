using DSP
using BenchmarkTools
using TensorOperations
using WilsonCowanModel
using Simulation73
using CalculatedTypes

function gaussian_dsp_conv(A::AbstractArray, gaussian_matrix)
    DSP.conv2(gaussian_matrix, A)
end

function tensor_conv!(output::Array{T,2}, conv_tensor::Array{T,4}, source::Array{T,2}) where T
    @tensor output[to_x, to_y] = conv_tensor[to_x, to_y, from_x, from_y] * source[from_x, from_y]
end

function test_tensor(full_conv_tensor::Array{T,4}, source, cutoff=1e-6) where T
    C_full = zero(source)
    C_reduced = zero(source)
    reduced_conv_tensor = full_conv_tensor
    reduced_conv_tensor[reduced_conv_tensor .< cutoff] .= zero(T)
    println("Full conv tensor")
    @btime tensor_conv!($C_full, $full_conv_tensor, $source)
    println("Reduced conv tensor")
    @btime tensor_conv!($C_reduced, $reduced_conv_tensor, $source)
    tensor_conv!(C_full, full_conv_tensor, source)
    tensor_conv!(C_reduced, reduced_conv_tensor, source)
    return (C_full, C_reduced)
end

function test_DSP(full_kernel::Array{T,2}, source, cutoff=1e-6) where T
    half_dx = ceil(Int, size(full_kernel,2) / 2)
    cutoff_dx = findfirst(full_kernel[:,half_dx] .>= cutoff)
    far_cutoff_dx = size(full_kernel,2) - cutoff_dx
    reduced_kernel = full_kernel[cutoff_dx:far_cutoff_dx, cutoff_dx:far_cutoff_dx]
    C_full = zero(source)
    C_reduced = zero(source)

    println("Full DSP convolution")
    @btime DSP.conv2($full_kernel, $source)
    println("Reduced DSP convolution")
    @btime DSP.conv2($reduced_kernel, $source)
    C_full = DSP.conv2(full_kernel, source)
    C_reduced = DSP.conv2(reduced_kernel, source)
    return (C_full, C_reduced)
end


function tests(cutoff::T=1e-6) where T
    A = 1.0
    σ = (2.5, 2.5)
    space = Calculated(Grid{Float64}(; n_points=(51,51), extent=(50.0,50.0)))
    gaussian_CONV_tensor = CalculatedTypes.calculate(GaussianConnectivity(A, σ), space)
    gaussian_matrix = WilsonCowanModel.exponential_decay_gaussian.(coordinates(space), A, Ref(σ), Ref(step(space)))
    source = rand(T, size(gaussian_matrix))

    tensor_results = test_tensor(gaussian_CONV_tensor, source)
    DSP_results = test_DSP(gaussian_matrix, source)

    # @show all.(tensor_results .≈ DSP_results)
    # @show sum(!(tensor_results[1] .≈ tensor_results[2]))

    return (tensor_results, DSP_results)
 end
