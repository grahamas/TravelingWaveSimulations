# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,julia//jl
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.0
#   kernelspec:
#     display_name: Julia 1.4.0
#     language: julia
#     name: julia-1.4
# ---

# +
using DrWatson
@quickactivate "TravelingWaveSimulations"
using WilsonCowanModel, TravelingWaveSimulations, BenchmarkTools

example = get_example("line_dos_effectively_sigmoid")

# +
function (wcm::WilsonCowanModel.WCMSpatialAction{T,N,P})(dA,A,p,t) where {T,N,P}
    wcm.stimulus(dA, A, t)
    wcm.connectivity(dA, A, t)
    wcm.nonlinearity(dA, A, t)
    @views for i in 1:P
        dAi = population(dA,i); Ai = population(A,i)
        dAi .*= wcm.β[i] .* (1.0 .- Ai)
        dAi .+= -wcm.α[i] .* Ai
        dAi ./= wcm.τ[i]
    end
end

@btime execute(example())

# +
function (wcm::WilsonCowanModel.WCMSpatialAction{T,N,P})(dA,A,p,t) where {T,N,P}
    wcm.stimulus(dA, A, t)
    wcm.connectivity(dA, A, t)
    wcm.nonlinearity(dA, A, t)
    for i in 1:P
        dAi = population(dA,i); Ai = population(A,i)
        dAi .*= wcm.β[i] .* (1.0 .- Ai)
        dAi .+= -wcm.α[i] .* Ai
        dAi ./= wcm.τ[i]
    end
end

@btime execute(example())

# +
function (wcm::WilsonCowanModel.WCMSpatialAction{T,N,P})(dA,A,p,t) where {T,N,P}
    @views wcm.stimulus(dA, A, t)
    @views wcm.connectivity(dA, A, t)
    @views wcm.nonlinearity(dA, A, t)
    @views for i in 1:P
        dAi = population(dA,i); Ai = population(A,i)
        dAi .*= wcm.β[i] .* (1.0 .- Ai)
        dAi .+= -wcm.α[i] .* Ai
        dAi ./= wcm.τ[i]
    end
end

@btime execute(example())
# -

a = Dict(:a => 1:2); b = Dict(:b => 2.0:5.0); c = Dict(:Orange => 2.0:30)
abc = merge(a, b, c)
names = keys(abc) |> Tuple
@show names
vals = values(abc)
(NamedTuple{names})(zip(Iterators.product(vals...)...))

NamedTuple{(:a,:b,:O)}((1.0, 3.0, 3))
