@testset "Sanity checks" begin
    @test isempty(detect_unbound_args(WilsonCowanModel))
end
