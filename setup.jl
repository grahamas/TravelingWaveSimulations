using Pkg

Pkg.add("DrWatson")
using DrWatson
quickactivate(@__DIR__, "WilsonCowanModel")

Pkg.develop(PackageSpec(path="src/WCMExamples"))
Pkg.develop(PackageSpec(path="src/WCMAnalysis"))
