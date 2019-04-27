using DrWatson
quickactivate(@__DIR__)

using WilsonCowanModel

example_name = Symbol(ARGS[1])

@eval Examples.$(example_name)()
