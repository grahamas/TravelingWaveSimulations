include(joinpath(@__DIR__, "init_proc.jl"))
using Distributed
@warn "Activating n_procs = $n_procs"
@assert nprocs() == 1
addprocs(n_procs)
@everywhere include(joinpath(@__DIR__, "init_proc.jl"))
