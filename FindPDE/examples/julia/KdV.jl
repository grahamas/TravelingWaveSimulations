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
#     display_name: Julia 1.4.1
#     language: julia
#     name: julia-1.4
# ---

using Revise
using FindPDE
using MLJ
#using Plots
using Makie
using AbstractPlotting
AbstractPlotting.inline!(true)

# # Generate a few solutions to each equation
#
# For any $c$, the function below is a very simple closed form solution to the KdV and advection equations.  However, it always solves the KdV equation regardless of the choice of $c$ while solutions with different values of $c$ will each solve an advection equation but they won't be the same one.

# Soliton solution to KdV/advection equations
function soliton(x,t,c,a)
    return c/2* cosh(sqrt(c)/2*(x-c*t-a))^(-2)
end

# +
c1 = 5.0
c2 = 1.0
c3 = 7.0
n_x = 512
n_t = 50
x = range(-10, stop=12, length=n_x)
dt = 0.025
dx = x[2]-x[1]
t = range(dt, stop=n_t*dt, length=n_t)

T = Float64
U1 = zeros(T, n_x, n_t)
U2 = zeros(T, n_x, n_t)
U3 = zeros(T, n_x, n_t)

for i_x in 1:n_x
    for j_t in 1:n_t
        U1[i_x, j_t] = soliton(x[i_x], t[j_t], c1, -3)
        U2[i_x, j_t] = soliton(x[i_x], t[j_t], c2, -1)
        U3[i_x, j_t] = soliton(x[i_x], t[j_t], c3, -4)
    end
end
# -

scene1 = surface(x, t, U1)
scale!(scene1, 1.0, 10.0, 1.0)
scene1

scene2 = surface(x, t, U2)
scale!(scene2, 1.0, 10.0, 1.0)
scene2

scene3 = surface(x, t, U3)
scale!(scene3, 1.0, 10.0, 1.0)
scene3

U1t, X1, space_terms_desc = FindPDE.build_linear_system(U1, dt, dx, 3, 2, 2)
# lam = 10^-5
# d_tol = 5
# w1 = FindPDE.train_STRidge(X1, U1t, lam, d_tol)
# w1_el = FindPDE.train_ElasticNet(X1, U1t, space_terms_desc, lam, d_tol)
# print("STR: $(sum(w1 .* space_terms_desc))\n")
# print("ElasticNet: $(sum(fitted_params(w1_el) .* space_terms_desc))\n")

U2t, X2, space_terms_desc = build_linear_system(U2, dt, dx, 3, 2, 2)
w2 = FindPDE.train_STRidge(X2, U2t, lam, d_tol)
sum(w2 .* space_terms_desc) |> string

U3t, X3, space_terms_desc = build_linear_system(U3, dt, dx, 3, 2, 2)
w3 = FindPDE.train_STRidge(X3, U3t, lam, d_tol, 100)
sum(w3 .* space_terms_desc) |> string

Ut = [U1t; U2t; U3t]
X = [X1; X2; X3]
lam = 10^-5
d_tol = 0.5
w = FindPDE.train_STRidge(X,Ut,lam,d_tol)#,100, 25)
sum(w .* space_terms_desc) |> string

Ut = [U1t; U2t; U3t]
X = [X1; X2; X3]
lam = 10^-5
d_tol = 0.5
w = FindPDE.train_STRidge(X,Ut,lam,d_tol)#,100, 25)
sum(w .* space_terms_desc) |> string

idx=9
@show space_terms_desc[idx]
surf = surface(x,t,reshape(X[:,idx], (length(x), length(t))))#, title=string(space_terms_desc[idx]))
scale!(surf, 1.0, 10.0, 1.0)
surf

using BenchmarkTools

names(df)

A = Float64[1 2 3 4 5; 10 20 30 40 50; 100 200 300 400 500; 1000 2000 3000 4000 5000]

info("ElasticNetRegressor", pkg="MLJLinearModels")

# # WCM Traveling Waves
#
# Try fitting simulation output (two populations)

using TravelingWaveSimulations


