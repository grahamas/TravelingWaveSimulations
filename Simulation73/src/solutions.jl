
Base.minimum(solution::DESolution) = minimum(map(minimum, solution.u))
Base.maximum(solution::DESolution) = maximum(map(maximum, solution.u))

using DiffEqBase: AbstractTimeseriesSolution

struct MissingSolution <: AbstractTimeseriesSolution{Missing,Missing,Missing} end
function Base.show(io::IO, A::MissingSolution)
  print(io,"MissingSolution")
end
function Base.show(io::IO, m::MIME"text/plain", A::MissingSolution)
  print(io,"MissingSolution")
end


struct BareSolution{S,N,U<:Array{<:Array{<:S,N}},X,T} <: AbstractTimeseriesSolution{S,N,U}
    u::U
    x::X
    t::T
    tslocation::Int
end
BareSolution(; u, x, t) = BareSolution(u,x,t)
BareSolution(u::U, x::X, t::T) where {S,N,ELT<:Union{S,Array{<:S,N}},U<:Array{ELT},X,T} = BareSolution{S,N,U,X,T}(u,x,t,0)
BareSolution(u::Missing, x::Missing, t::Missing) = MissingSolution()
function Base.getproperty(bs::BareSolution, sym::Symbol)
    if sym ∈ [:dense, :prob]
        return false
    else
        return Base.getfield(bs, sym)
    end
end
function DiffEqBase.add_labels!(labels,x,dims,sol::BareSolution)
    return
end
timepoints(bs::BareSolution) = bs.t
coordinates(bs::BareSolution) = bs.x
function Base.show(io::IO, A::BareSolution)
  print(io,"t: ")
  show(io, A.t)
  println(io)
  print(io,"u: ")
  show(io, A.u)
end
function Base.show(io::IO, m::MIME"text/plain", A::BareSolution)
  print(io,"t: ")
  show(io,m,A.t)
  println(io)
  print(io,"u: ")
  show(io,m,A.u)
end
@generated function Simulation73.population_timepoint(solution::BareSolution{T,NP}, pop_dx::Int, time_dx::Int) where {T,NP}
    N = NP - 1 # N + pops
    colons = [:(:) for i in 1:N]
    :(solution[time_dx][$(colons...), pop_dx])
end
