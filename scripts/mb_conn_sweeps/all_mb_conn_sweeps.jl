
let sweep_dir = ""
    include(joinpath(@__DIR__, sweep_dir, "sweep_A_monotonic.jl"))
    include(joinpath(@__DIR__, sweep_dir, "sweep_S_monotonic.jl"))
    include(joinpath(@__DIR__, sweep_dir, "sweep_A_blocking.jl"))
    include(joinpath(@__DIR__, sweep_dir, "sweep_S_blocking.jl"))
end
