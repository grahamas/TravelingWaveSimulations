
@testset "Pre-made prototypes test" begin
    prototype_name = "dos_ring"
    line_prototype = @test_nowarn get_prototype(prototype_name)
    these_mods = (Aee=40.0,Aei=200.0, Aie=73.0, 
        blocking_θE=25.0,blocking_θI=25.0,
        firing_θE=6.0,firing_θI=7.0, 
        other_opts=Dict())
    (mod_name, exec) = @test_nowarn execute_single_modification(line_prototype, these_mods)
    wp = @test_nowarn ExecutionClassifications(exec)
end
