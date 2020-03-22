
function mods_to_pkeys(modifications)::Array{Symbol,1}
    pkeys = keys(modifications[1]) |> collect
    disallowed_keys = [:algorithm, :u, :x, :t, :n, :n_points, :extent, :save_idxs]
    pkeys = filter((x) -> !(x in disallowed_keys), pkeys)
    pkeys = filter((x) -> !(modifications[1][x] === nothing), pkeys)
end

function getkeys(d, keys)
    [d[key] for key in keys]
end

function extract_data_namedtuple(execution::Execution)
    soln = execution.solution
    u = soln.u
    t = soln.t
    x = coordinates(space(execution)) |> collect
    return execution.simulation.global_reduction((u=u, t=t, x=x))
end


function extract_data_namedtuple(execution::ReducedExecution)
    execution.simulation.global_reduction(extract_data_namedtuple(execution.saved_values))
end

function extract_data_namedtuple(execution::AugmentedExecution)
    soln = execution.solution
    u = soln.u
    t = soln.t
    x = coordinates(space(execution))
    execution.simulation.global_reduction((u=u, t=t, x=x, extract_data_namedtuple(execution.saved_values)...))
end

function extract_data_namedtuple(sv::SavedValues{T,<:Array{<:Wavefront,1}}) where T
    return (wavefronts=sv.saveval, wavefronts_t=sv.t)
end

function extract_data_namedtuple(execution::Missing)
    return execution.simulation.global_reduction((u=missing, t=missing, x=missing))
end

function extract_params_tuple(modification, pkeys)
    return NamedTuple{Tuple(pkeys)}(getkeys(modification, pkeys))
end

function init_results_tuple(pkeys::Array{Symbol}, results::NamedTuple{NAMES,TYPES}) where {NAMES,TYPES}
    N = length(pkeys)
    ALL_NAMES = Tuple([pkeys...,NAMES...])
    ARR_TYPES = Tuple{[Array{Float64,1} for _ in 1:N]..., [Array{Union{Missing,TYPE},1} for TYPE in TYPES.parameters]...}
    empty_arrs = [[Float64[] for _ in 1:N]..., [Union{Missing,TYPE}[] for TYPE in TYPES.parameters]...]
    NamedTuple{ALL_NAMES,ARR_TYPES}(empty_arrs)
end

function init_failures_tuple(pkeys::Array{Symbol})
    N = length(pkeys)
    NamedTuple{Tuple(pkeys),NTuple{N,Array{Float64,1}}}([Float64[] for _ in 1:N])
end

Base.getindex(nt::NamedTuple, dx::Array{Symbol}) = getindex.(Ref(nt), dx)

function merge_ddb(::Nothing, tbl)
    return distribute(tbl, nworkers())
end
function merge_ddb(ddb::JuliaDB.DNDSparse, tbl::NDSparse)
    @error "Sparse unsupported"
    return merge(ddb, tbl)
end
function merge_ddb(ddb::JuliaDB.DIndexedTable, tbl::IndexedTable)
    return merge(ddb, tbl)
end

function push_namedtuple!(::Nothing, mods::NamedTuple{NAMES,TYPES}) where {NAMES,TYPES}
    arrd_TYPES = Tuple{[Array{T,1} for T in TYPES.parameters]...}
    NamedTuple{NAMES, arrd_TYPES}([[val] for val in values(mods)])
end
function push_namedtuple!(tup::NamedTuple, mods)
    for mod in pairs(mods)
        push!(tup[mod[1]], mod[2])
    end
end
push_namedtuple!(::Nothing, ::Nothing) = nothing
push_namedtuple!(nt::NamedTuple, ::Nothing) = nt