@testset "Sanity checks" begin
    @test_broken isempty(detect_unbound_args(Simulation73))
end
