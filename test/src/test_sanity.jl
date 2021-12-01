using Test, TravelingWaveSimulations

@testset "Sanity checks" begin
    @test isempty(detect_unbound_args(TravelingWaveSimulations))
end
