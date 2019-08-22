@testset "High-level tests" begin
    @test based_on_example(; example_name="sigmoid_normal", analyses=["radial_slice"], modifications=["iiS=1.7"]) != 0
end
