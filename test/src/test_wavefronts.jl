
using Test, AxisIndices, Statistics
using Simulation73: population_repeat
using TravelingWaveSimulations: get_velocities, substantial_fronts, 
    detect_all_fronts, link_persistent_fronts, slope_loc, slope_val,
    get_prototype, execute_single_modification,
    Wavefront, MinimalPropagationClassification,
    already_reduced_to_min_propagation_cls

# from NeuralModels/src/nonlinearity.jl
function _simple_sigmoid(x, a, theta)
    1.0 / (1 + exp(-a * (x - theta)))
end

function wavefront(x, center, steepness, height, increasing::Val{true})
    height * _simple_sigmoid(x, steepness, center)
end
function wavefront(x, center, steepness, height, increasing::Val{false})
    height * (1 - _simple_sigmoid(x, steepness, center))
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

traveling_parameters = (velocity_threshold=1e-4, n_traveling_frames_threshold=20)

xs = -20.0:0.1:20.0
ts = 0.0:0.1:20.0
front_steepness = 1.
detection_slope_min = 1e-4
traveling_increasing_wavefront_fn = make_traveling_wavefront(; velocity=1., x_0=0., steepness=front_steepness, height=1.)
traveling_decreasing_wavefront_fn = make_traveling_wavefront(; velocity=1., x_0=0., steepness=front_steepness, height=1., increasing=false)
traveling_pair_fn = make_traveling_wavefront_pair(; velocity=1., x_0=-3., width=5., steepness=front_steepness, height=1.)

@testset "Static front detection" begin
    single_front_increasing = AxisArray(traveling_increasing_wavefront_fn.(xs, 0.0), xs)
    single_front_decreasing = AxisArray(traveling_decreasing_wavefront_fn.(xs, 0.0), xs)
    trunc_front_loc = xs[end] - front_steepness
    single_front_truncated = AxisArray(traveling_decreasing_wavefront_fn.(xs, trunc_front_loc), xs)
    pair_of_wavefronts = AxisArray(traveling_pair_fn.(xs, 0.0), xs)

    @testset "Substantial front detection" begin
        sfi_detected_fronts = substantial_fronts(single_front_increasing, false)
        @test length(sfi_detected_fronts) == 1
        @test isapprox(slope_loc(sfi_detected_fronts[1]), 0.0; atol=sqrt(eps()))
        @test slope_val(sfi_detected_fronts[1]) > 0.0

        sfd_detected_fronts = substantial_fronts(single_front_decreasing, false)
        @test length(sfd_detected_fronts) == 1
        @test isapprox(slope_loc(sfd_detected_fronts[1]), 0.0; atol=sqrt(eps()))
        @test slope_val(sfd_detected_fronts[1]) < 0.0
    
        pair_detected_fronts = substantial_fronts(pair_of_wavefronts, false)
        @test length(pair_detected_fronts) == 2
        @test slope_loc(pair_detected_fronts[1]) < 0.0
        @test slope_val(pair_detected_fronts[1]) > 0.0
        @test slope_loc(pair_detected_fronts[2]) > 2.0
        @test slope_val(pair_detected_fronts[2]) < 0.0

        trunc_detected_fronts = substantial_fronts(single_front_truncated, false)
        @test length(trunc_detected_fronts) == 1
        @test isapprox(slope_loc(trunc_detected_fronts[1]), trunc_front_loc; atol=sqrt(eps()))
        @test slope_val(trunc_detected_fronts[1]) < 0.0
    end

end


@testset "Toy traveling fronts detection (persistent fronts method)" begin

    traveling_increasing_wavefront = [AxisArray(traveling_increasing_wavefront_fn.(xs, t), xs) for t in ts]
    traveling_decreasing_wavefront = [AxisArray(traveling_decreasing_wavefront_fn.(xs, t), xs) for t in ts]
    traveling_pair = [AxisArray(traveling_pair_fn.(xs, t), xs) for t in ts]

    motionless_decreasing_wavefront = [AxisArray(traveling_decreasing_wavefront_fn.(xs, 0.0), xs) for t in ts]
    motionless_pair = [AxisArray(traveling_pair_fn.(xs, 0.0), xs) for t in ts]

    @testset "Persistent fronts method" begin
        traveling_increasing_persistent_fronts = link_persistent_fronts(substantial_fronts.(traveling_increasing_wavefront, false, detection_slope_min), ts)
        @test length(traveling_increasing_persistent_fronts) == 1
        @show maximum(get_velocities(traveling_increasing_persistent_fronts[1]))
        @test mean(get_velocities(traveling_increasing_persistent_fronts[1])) ≈ 1.

        traveling_decreasing_persistent_fronts = link_persistent_fronts(substantial_fronts.(traveling_decreasing_wavefront, false, detection_slope_min), ts)
        @test length(traveling_decreasing_persistent_fronts) == 1
        @test mean(get_velocities(traveling_decreasing_persistent_fronts[1])) ≈ 1.

        traveling_pair_persistent_fronts = link_persistent_fronts(substantial_fronts.(traveling_pair, false, detection_slope_min), ts)
        @test length(traveling_pair_persistent_fronts) == 2
        @test mean(get_velocities(traveling_pair_persistent_fronts[1])) ≈ 1.

        motionless_decreasing_persistent_fronts = link_persistent_fronts(substantial_fronts.(motionless_decreasing_wavefront, false, detection_slope_min), ts)
        @test length(motionless_decreasing_persistent_fronts) == 1
        @test mean(get_velocities(motionless_decreasing_persistent_fronts[1])) ≈ 0.

        motionless_pair_persistent_fronts = link_persistent_fronts(substantial_fronts.(motionless_pair, false, detection_slope_min), ts)
        @test length(motionless_pair_persistent_fronts) == 2
        @test mean(get_velocities(motionless_pair_persistent_fronts[1])) ≈ 0.
    end
end

@testset "Faux-execution classifications" begin
    traveling_increasing_wavefront = [population_repeat(AxisArray(traveling_increasing_wavefront_fn.(xs, t), xs), 2) for t in ts]
    traveling_decreasing_wavefront = [population_repeat(AxisArray(traveling_decreasing_wavefront_fn.(xs, t), xs), 2) for t in ts]
    traveling_pair = [population_repeat(AxisArray(traveling_pair_fn.(xs, t), xs), 2) for t in ts]

    motionless_decreasing_wavefront = [population_repeat(AxisArray(traveling_decreasing_wavefront_fn.(xs, 0.0), xs), 2) for t in ts]
    motionless_pair = [population_repeat(AxisArray(traveling_pair_fn.(xs, 0.0), xs), 2) for t in ts]

    @testset "MinimalPropagationClassification" begin
        minprop(wavefront) = MinimalPropagationClassification(wavefront, ts,
                            xs, true;  min_dist_for_propagation = 10.).has_propagation
        
        @test minprop(traveling_increasing_wavefront) == true
        @test minprop(traveling_decreasing_wavefront) == true
        @test minprop(traveling_pair) == true
        @test minprop(motionless_decreasing_wavefront) == false
        @test minprop(motionless_pair) == false
    end

end


@testset "WCM execution classifications" begin
    wcm_n_lattice = 512
    wcm_x_lattice = 1800.0
    wcm_dx = wcm_x_lattice / wcm_n_lattice
    general_mods = (n_lattice=wcm_n_lattice, x_lattice=wcm_x_lattice, 
        step_reduction=nothing)
    wcm_dist_for_prop = wcm_x_lattice * 0.1
    min_prop_parameters = (min_dist_for_propagation = wcm_dist_for_prop, const_jitter = wcm_dx * 2,
            vel_jitter = 1.)

    # view with: Simulation73Plotting.animate_execution("$(tempname()).mp4", EXEC_VARIABLE)

    @testset "Localized parameterizations" begin

        prototype_name = "ring_blocking"
        line_prototype = get_prototype(prototype_name)
        localized_mods = (Aee=40.0,Aei=200.0, Aie=73.0, 
            blocking_θI=25.0,
            θE=6.0,firing_θI=7.0, 
            other_opts=Dict(), save_on=true)
        (_, localized_exec) = execute_single_modification(line_prototype, merge(general_mods, localized_mods))

        @testset "MinimalPropagationClassification" begin
            min_prop_cls = MinimalPropagationClassification(localized_exec; min_prop_parameters...)
            @test min_prop_cls.has_propagation == false
        end
    end


    @testset "Propagating parameterizations" begin
        prototype_name = "ring_blocking"
        line_prototype = get_prototype(prototype_name)
        expanding_front_mods = (Aee=150.0,Aei=50.0, Aie=73.0, 
            blocking_θI=25.0,
            θE=6.0,firing_θI=7.0, save_on=true,
            other_opts=Dict())
        (_, expanding_front_exec) = execute_single_modification(line_prototype, merge(general_mods, expanding_front_mods))

        @testset "MinimalPropagationClassification" begin
            expanding_front_min_prop_cls = MinimalPropagationClassification(expanding_front_exec; min_prop_parameters...)
            @test expanding_front_min_prop_cls.has_propagation == true
        end

        
    end

    @testset "MinimalPropagationClassification on-line defaults" begin
        prototype_name = "ring_blocking"
        line_prototype = get_prototype(prototype_name)
        localized_mods = (Aee=40.0,Aei=200.0, Aie=73.0, 
            blocking_θI=25.0,
            θE=6.0,firing_θI=7.0, 
            other_opts=Dict())
        (_, localized_exec) = execute_single_modification(line_prototype, merge(general_mods, localized_mods))

        localized_min_prop_cls = already_reduced_to_min_propagation_cls(localized_exec.solution)
        @test localized_min_prop_cls.propagation.has_propagation == false

        prototype_name = "ring_blocking"
        line_prototype = get_prototype(prototype_name)
        expanding_front_mods = (Aee=150.0,Aei=50.0, Aie=73.0, 
            blocking_θI=25.0,
            θE=6.0,firing_θI=7.0,
            other_opts=Dict())
        (_, expanding_front_exec) = execute_single_modification(line_prototype, merge(general_mods, expanding_front_mods))

        expanding_front_min_prop_cls = already_reduced_to_min_propagation_cls(expanding_front_exec.solution)
        @test expanding_front_min_prop_cls.propagation.has_propagation == true
    end

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