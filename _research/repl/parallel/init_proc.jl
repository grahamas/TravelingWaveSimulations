using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using DrWatson
include(projectdir("repl", "setup", "basic.jl"))