using Pkg

Pkg.add("DrWatson")
using DrWatson
quickactivate(@__DIR__, "WilsonCowanModel")

push!(LOAD_PATH, "./src/")
Pkg.develop("WCMExamples")
Pkg.develop("WCMAnalysis")
