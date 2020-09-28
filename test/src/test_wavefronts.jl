
using Test, AxisIndices, TravelingWaveSimulations

# from NeuralModels/src/nonlinearity.jl
function _simple_sigmoid_fn(x, a, theta)
    1.0 / (1 + exp(-a * (x - theta)))
end

function wavefront(x, center, steepness, height, increasing::Val{true})
    height * _simple_sigmoid_fn(x, steepness, center)
end
function wavefront(x, center, steepness, height, increasing::Val{false})
    height * (1 - _simple_sigmoid_fn(x, steepness, center))
end

function make_traveling_wavefront(;velocity, x_0, steepness, height, increasing::Bool=true)
    function _traveling_wavefront(x, t)
        wavefront(x, x_0 + velocity * t, steepness, height, Val(increasing))
    end
end

function make_traveling_wavefront_pair(;velocity, x_0, width, steepness, height)
    function _traveling_wavefront_pair(x, t)
        wavefront(x, x_0 + velocity * t, steepness, height, Val(true)) - wavefront(x, x_0 + width + velocity * t, steepness, height, Val(true))
    end
end

xs = -10.0:0.1:10.0
ts = 0.0:0.1:1.0
traveling_increasing_wavefront_fn = make_traveling_wavefront(; velocity=1., x_0=0., steepness=1., height=1.)
traveling_decreasing_wavefront_fn = make_traveling_wavefront(; velocity=1., x_0=0., steepness=1., height=1., increasing=false)
traveling_pair_fn = make_traveling_wavefront_pair(; velocity=1., x_0=-3., width=5., steepness=1., height=1.)

@testset "Static front detection" begin
    single_front_increasing = AxisArray(traveling_increasing_wavefront_fn.(xs, 0.0), xs)
    single_front_decreasing = AxisArray(traveling_decreasing_wavefront_fn.(xs, 0.0), xs)
    pair_of_wavefronts = AxisArray(traveling_pair_fn.(xs, 0.0), xs)

    @testset "Raw front detection" begin
        sfi_detected_fronts = TravelingWaveSimulations.detect_all_fronts(single_front_increasing)
        @test length(sfi_detected_fronts) == 1
        @test TravelingWaveSimulations.slope_loc(sfi_detected_fronts[1]) == 0.0
        @test TravelingWaveSimulations.slope_val(sfi_detected_fronts[1]) > 0.0

        sfd_detected_fronts = TravelingWaveSimulations.detect_all_fronts(single_front_decreasing)
        @test length(sfd_detected_fronts) == 1
        @test TravelingWaveSimulations.slope_loc(sfd_detected_fronts[1]) == 0.0
        @test TravelingWaveSimulations.slope_val(sfd_detected_fronts[1]) < 0.0
    
        pair_detected_fronts = TravelingWaveSimulations.detect_all_fronts(pair_of_wavefronts)
        @test length(pair_detected_fronts) == 2
        @test TravelingWaveSimulations.slope_loc(pair_detected_fronts[1]) < 0.0
        @test TravelingWaveSimulations.slope_val(pair_detected_fronts[1]) > 0.0
        @test TravelingWaveSimulations.slope_loc(pair_detected_fronts[2]) > 2.0
        @test TravelingWaveSimulations.slope_val(pair_detected_fronts[2]) < 0.0
    end

    @testset "Substantial front detection" begin
        sfi_detected_fronts = TravelingWaveSimulations.substantial_fronts(single_front_increasing)
        @test length(sfi_detected_fronts) == 1
        @test TravelingWaveSimulations.slope_loc(sfi_detected_fronts[1]) == 0.0
        @test TravelingWaveSimulations.slope_val(sfi_detected_fronts[1]) > 0.0

        sfd_detected_fronts = TravelingWaveSimulations.substantial_fronts(single_front_decreasing)
        @test length(sfd_detected_fronts) == 1
        @test TravelingWaveSimulations.slope_loc(sfd_detected_fronts[1]) == 0.0
        @test TravelingWaveSimulations.slope_val(sfd_detected_fronts[1]) < 0.0
    
        pair_detected_fronts = TravelingWaveSimulations.substantial_fronts(pair_of_wavefronts)
        @test length(pair_detected_fronts) == 2
        @test TravelingWaveSimulations.slope_loc(pair_detected_fronts[1]) < 0.0
        @test TravelingWaveSimulations.slope_val(pair_detected_fronts[1]) > 0.0
        @test TravelingWaveSimulations.slope_loc(pair_detected_fronts[2]) > 2.0
        @test TravelingWaveSimulations.slope_val(pair_detected_fronts[2]) < 0.0
    end

end

@testset "Traveling fronts detection" begin

end