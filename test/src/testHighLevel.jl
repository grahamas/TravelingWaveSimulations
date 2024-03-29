using Test, TravelingWaveSimulations

@testset "Pre-made prototypes test" begin
    prototype_name = "full_dynamics_blocking"
    line_prototype = @test_nowarn get_prototype(prototype_name)
    these_mods = (Aee=40.0,Aei=200.0, Aie=73.0, 
        blocking_θI=25.0,
        θE=6.0,firing_θI=7.0, 
        n_lattice=128, x_lattice=300.0)
    (mod_name, exec) = @test_nowarn execute_single_modification(line_prototype, these_mods)
    # wp = @test_nowarn ExecutionClassifications(exec; velocity_threshold=1e-4, n_traveling_frames_threshold=50)
    # mp = @test_nowarn MinimalPropagationClassification(exec; min_dist_for_propagation = 30.)
end
