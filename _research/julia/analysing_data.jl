# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all
#     notebook_metadata_filter: all,-language_info
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Julia 1.3.0
#     language: julia
#     name: julia-1.3
# ---

# %%
using Simulation73, TravelingWaveSimulations
using DiffEqBase: AbstractTimeseriesSolution
# abstract type AbstractSolution{S} end
# const SomeSolution{S} = Union{AbstractTimeseriesSolution{S}, AbstractSolution{S}}
struct BareSolution{S,N,U<:Array{<:Array{<:S,N}},X,T} <: AbstractTimeseriesSolution{S,N}
    u::U
    x::X
    t::T
end
BareSolution(; u::U, x::X, t::T) where {S,N,U<:Array{<:Array{<:S,N}},X,T} = BareSolution{S,N,U,X,T}(u,x,t)
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

# %%
using Lazy, JuliaDB

data_root = joinpath(homedir(), "sim_data")

function _load_data(sim_path)
    @warn "loading"
   JuliaDB.load(sim_path)
end
"Load most recent simulation"
function load_data(data_root, example_name)
    nsims = length(readdir(joinpath(data_root, example_name)))
    load_data(data_root, example_name, nsims)
end
"Load nth simulation, ordered by time"
function load_data(data_root, example_name, nth::Int)
    ex_path = joinpath(data_root, example_name)
    sims = readdir(ex_path)
    sorted_sims = sort(joinpath.(Ref(ex_path), sims), by=mtime)
    sim_path = sorted_sims[nth]
    return (get_example(example_name), @lazy map(_load_data, joinpath.(Ref(sim_path), readdir(sim_path))))
end

macro ifsomething(ex)
    quote
        result = $(esc(ex))
        result === nothing && return nothing
        result
    end
end
struct MultiDBRowIter
    dbs
end
function _get_row(dbs, db, dbs_state, row_state...)
    row_tuple = iterate(db, row_state...)
    while row_tuple === nothing
        db, dbs_state = @ifsomething iterate(dbs)
        row_tuple = iterate(db)
    end
    (row, row_state) = row_tuple
    return (row, (dbs_state, db, row_state))
end
function iterate(it::MultiDBRowIter)
    (db, dbs_state) = @ifsomething iterate(it.dbs)
    return _get_row(it.dbs, db, dbs_state)
end
function iterate(it::MultiDBRowIter, (db_state, db, row_state))
    return _get_row(it.dbs, db, dbs_state, row_state)
end
    
    
    



# %%
(example, dbs) = load_data(data_root, "sigmoid_normal_fft");

# %%
n_db = 5
keys1 = select(dbs[5], Keys())
vals1 = select(dbs[5], Not(Keys()))
first(dbs[1])

# %%
n=3000
mdl1 = example(; keys1[n]...)
anim1 = custom_animate(Execution(mdl1, BareSolution(; vals1[n]...)))

# %%
mp4(anim1, "tmp/tmp.mp4")
