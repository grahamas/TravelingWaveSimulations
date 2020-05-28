
@testset "Pre-made examples test" begin
    example_name = "reduced_line_dos_effectively_sigmoid"
    line_example = @test_nowarn get_example(example_name)
    these_mods = (Aee=40.0,Aei=200.0, Aie=73.0, 
        blocking_θE=25.0,blocking_θI=25.0,
        firing_θE=6.0,firing_θI=7.0, 
        other_opts=Dict())
    (mod_name, exec) = @test_nowarn execute_single_modification(line_example, these_mods)
    wp = @test_nowarn get_wave_properties(exec)
end
