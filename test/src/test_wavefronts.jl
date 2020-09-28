
using Test, AxisIndices, Statistics
using Simulation73: population_repeat
using TravelingWaveSimulations: get_velocities, substantial_fronts, 
    detect_all_fronts, link_persistent_fronts, slope_loc, slope_val,
    ExecutionClassifications, get_prototype, execute_single_modification,
    Wavefront

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
        sfi_detected_fronts = detect_all_fronts(single_front_increasing)
        @test length(sfi_detected_fronts) == 1
        @test slope_loc(sfi_detected_fronts[1]) == 0.0
        @test slope_val(sfi_detected_fronts[1]) > 0.0

        sfd_detected_fronts = detect_all_fronts(single_front_decreasing)
        @test length(sfd_detected_fronts) == 1
        @test slope_loc(sfd_detected_fronts[1]) == 0.0
        @test slope_val(sfd_detected_fronts[1]) < 0.0
    
        pair_detected_fronts = detect_all_fronts(pair_of_wavefronts)
        @test length(pair_detected_fronts) == 2
        @test slope_loc(pair_detected_fronts[1]) < 0.0
        @test slope_val(pair_detected_fronts[1]) > 0.0
        @test slope_loc(pair_detected_fronts[2]) > 2.0
        @test slope_val(pair_detected_fronts[2]) < 0.0
    end

    @testset "Substantial front detection" begin
        sfi_detected_fronts = substantial_fronts(single_front_increasing)
        @test length(sfi_detected_fronts) == 1
        @test slope_loc(sfi_detected_fronts[1]) == 0.0
        @test slope_val(sfi_detected_fronts[1]) > 0.0

        sfd_detected_fronts = substantial_fronts(single_front_decreasing)
        @test length(sfd_detected_fronts) == 1
        @test slope_loc(sfd_detected_fronts[1]) == 0.0
        @test slope_val(sfd_detected_fronts[1]) < 0.0
    
        pair_detected_fronts = substantial_fronts(pair_of_wavefronts)
        @test length(pair_detected_fronts) == 2
        @test slope_loc(pair_detected_fronts[1]) < 0.0
        @test slope_val(pair_detected_fronts[1]) > 0.0
        @test slope_loc(pair_detected_fronts[2]) > 2.0
        @test slope_val(pair_detected_fronts[2]) < 0.0
    end

end


@testset "Traveling fronts detection" begin
    traveling_increasing_wavefront = [AxisArray(traveling_increasing_wavefront_fn.(xs, t), xs) for t in ts]
    traveling_decreasing_wavefront = [AxisArray(traveling_decreasing_wavefront_fn.(xs, t), xs) for t in ts]
    traveling_pair = [AxisArray(traveling_pair_fn.(xs, t), xs) for t in ts]

    motionless_decreasing_wavefront = [AxisArray(traveling_decreasing_wavefront_fn.(xs, 0.0), xs) for t in ts]
    motionless_pair = [AxisArray(traveling_pair_fn.(xs, 0.0), xs) for t in ts]

    traveling_increasing_persistent_fronts = link_persistent_fronts(substantial_fronts.(traveling_increasing_wavefront), ts)
    @test length(traveling_increasing_persistent_fronts) == 1
    @test mean(get_velocities(traveling_increasing_persistent_fronts[1])) ≈ 1.

    traveling_decreasing_persistent_fronts = link_persistent_fronts(substantial_fronts.(traveling_decreasing_wavefront), ts)
    @test length(traveling_decreasing_persistent_fronts) == 1
    @test mean(get_velocities(traveling_decreasing_persistent_fronts[1])) ≈ 1.

    traveling_pair_persistent_fronts = link_persistent_fronts(substantial_fronts.(traveling_pair), ts)
    @test length(traveling_pair_persistent_fronts) == 2
    @test mean(get_velocities(traveling_pair_persistent_fronts[1])) ≈ 1.

    motionless_decreasing_persistent_fronts = link_persistent_fronts(substantial_fronts.(motionless_decreasing_wavefront), ts)
    @test length(motionless_decreasing_persistent_fronts) == 1
    @test mean(get_velocities(motionless_decreasing_persistent_fronts[1])) ≈ 0.

    motionless_pair_persistent_fronts = link_persistent_fronts(substantial_fronts.(motionless_pair), ts)
    @test length(motionless_pair_persistent_fronts) == 2
    @test mean(get_velocities(motionless_pair_persistent_fronts[1])) ≈ 0.
end

@testset "Faux-execution classifications" begin
    traveling_increasing_wavefront = [population_repeat(AxisArray(traveling_increasing_wavefront_fn.(xs, t), xs), 2) for t in ts]
    traveling_decreasing_wavefront = [population_repeat(AxisArray(traveling_decreasing_wavefront_fn.(xs, t), xs), 2) for t in ts]
    traveling_pair = [population_repeat(AxisArray(traveling_pair_fn.(xs, t), xs), 2) for t in ts]

    motionless_decreasing_wavefront = [population_repeat(AxisArray(traveling_decreasing_wavefront_fn.(xs, 0.0), xs), 2) for t in ts]
    motionless_pair = [population_repeat(AxisArray(traveling_pair_fn.(xs, 0.0), xs), 2) for t in ts]

    traveling_increasing_wavefront_classification = ExecutionClassifications(substantial_fronts.(traveling_increasing_wavefront), 
                                                                            ts, xs, traveling_increasing_wavefront[end]; origin_radius=1.0)
    @test traveling_increasing_wavefront_classification.has_propagation == true
    
    traveling_decreasing_wavefront_classification = ExecutionClassifications(substantial_fronts.(traveling_decreasing_wavefront), 
                                                                            ts, xs, traveling_decreasing_wavefront[end]; origin_radius=1.0)
    @test traveling_decreasing_wavefront_classification.has_propagation == true
    
    traveling_pair_classification = ExecutionClassifications(substantial_fronts.(traveling_pair), 
                                                                            ts, xs, traveling_pair[end]; origin_radius=1.0)
    @test traveling_pair_classification.has_propagation == true
    
    motionless_decreasing_classification = ExecutionClassifications(substantial_fronts.(motionless_decreasing_wavefront), 
                                                                            ts, xs, motionless_decreasing_wavefront[end]; origin_radius=1.0)
    @test motionless_decreasing_classification.has_propagation == false
    
    motionless_pair_classification = ExecutionClassifications(substantial_fronts.(motionless_pair), 
                                                                            ts, xs, motionless_pair[end]; origin_radius=1.0)
    @test motionless_pair_classification.has_propagation == false
end


@testset "WCM execution classifications" begin

    # view with: Simulation73Plotting.animate_execution("$(tempname()).mp4", EXEC_VARIABLE)

    prototype_name = "ring_blocking"
    line_prototype = get_prototype(prototype_name)
    localized_mods = (Aee=40.0,Aei=200.0, Aie=73.0, 
        blocking_θI=25.0,
        θE=6.0,firing_θI=7.0, 
        n_lattice=256, x_lattice=600.0,
        other_opts=Dict())
    (_, localized_exec) = execute_single_modification(line_prototype, localized_mods)

    # @testset "WaveClassification internals on localized WCM example" begin
    #     let exec = localized_exec
    #         l_frame_fronts = exec.saved_values.saveval
    #         ts = exec.saved_values.t
    #         final_frame = exec.solution.u[end]
    #         xs = axes_keys(final_frame)[1]
    #         persistent_fronts = link_persistent_fronts(l_frame_fronts, ts)
    #         pf_measurements = SpatiotemporalWaveMeasurements.(persistent_fronts)
    #         @show pf_measurements .|> pfm -> pfm.velocities
    #         @test all(pf_measurements .|> pfm -> all(pfm.velocities .≈ 0.0))
    #     end
    # end

    localized_exec_cls = ExecutionClassifications(localized_exec)
    @test localized_exec_cls.has_propagation == false

    prototype_name = "ring_blocking"
    line_prototype = get_prototype(prototype_name)
    expanding_front_mods = (Aee=150.0,Aei=50.0, Aie=73.0, 
        blocking_θI=25.0,
        θE=6.0,firing_θI=7.0, 
        n_lattice=256, x_lattice=600.0,
        other_opts=Dict())
    (_, expanding_front_exec) = execute_single_modification(line_prototype, expanding_front_mods)   
    expanding_front_exec_cls = ExecutionClassifications(expanding_front_exec)
    @test expanding_front_exec_cls.has_propagation == true

    # # This case is ambiguous
    # prototype_name = "ring_blocking"
    # line_prototype = get_prototype(prototype_name)
    # truncated_expanding_front_mods = (Aee=150.0,Aei=50.0, Aie=73.0, 
    #     blocking_θI=25.0,
    #     θE=6.0,firing_θI=7.0, 
    #     n_lattice=256, x_lattice=200.0,
    #     other_opts=Dict())
    # (_, truncated_expanding_front_exec) = execute_single_modification(line_prototype, truncated_expanding_front_mods)
    # truncated_expanding_front_exec_cls = ExecutionClassifications(truncated_expanding_front_exec)
    # @test truncated_expanding_front_exec_cls.has_propagation == true
end