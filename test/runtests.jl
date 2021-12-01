using SafeTestsets

@time begin
    # tested by wavefront 
    @time @safetestset "High level integration test" begin include("src/testHighLevel.jl") end
    @time @safetestset "Sanity checks" begin include("src/test_sanity.jl") end
end
# TODO: test parallel
