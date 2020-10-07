using SafeTestsets

@time begin
    # tested by wavefront #@time @safetestset "High level integration test" begin include("src/testHighLevel.jl") end
    @time @safetestset "Sanity checks" begin include("src/test_sanity.jl") end
    @time @safetestset "Wavefront detection algorithm tests" begin include("src/test_wavefronts.jl") end
end
# TODO: test parallel
