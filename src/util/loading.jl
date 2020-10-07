using JuliaDB

macro ifsomething(ex)
    quote
        result = $(esc(ex))
        result === nothing && return nothing
        result
    end
end

function _load_data(sim_path)
    @show sim_path
   JuliaDB.load(sim_path)
end
struct MultiDB
    fns::Array{String}
    function MultiDB(fns::AbstractArray)
        if isempty(fns)
            @warn "Empty MultiDB initialized."
        end
        new(fns)
    end
end
is_valid_db_path(::Nothing) = false
function is_valid_db_path(fn)
    splitted = splitext(fn)
    if length(splitted) == 1
        return false
    elseif splitted[2] != ".jdb"
        return false
    elseif basename(splitted[1]) == "failures"
        return false
    else
        return true
    end
end
function is_valid_mdb_path(path)
    any(map(is_valid_db_path, joinpath.(Ref(path), readdir(path))))
end
function get_sorted_db_subpaths(path)
    filter(is_valid_db_path, joinpath.(Ref(path), readdir(path)))
end
function get_sorted_mdb_subpaths(path)
    @show path
    mdb_subpaths = filter(is_valid_mdb_path, joinpath.(Ref(path), readdir(path)))
    sorted_mdb_subpaths = sort(mdb_subpaths, by=mtime)
    return sorted_mdb_subpaths
end
function MultiDB(path::AbstractString)
    return MultiDB(get_sorted_db_subpaths(path))
end
Base.length(mdb::MultiDB) = length(mdb.fns)

function Base.iterate(mdb::MultiDB, fns_state...)
    next_fn = nothing
    while !is_valid_db_path(next_fn)
        (next_fn, fns_state) = @ifsomething iterate(mdb.fns, fns_state...)
    end
    (_load_data(next_fn), fns_state)
end

function parse_mod_val(val)
    if val == "nothing"
        return nothing
    elseif val == "false"
        return false
    elseif val == "true"
        return true
    else
        parse(Float64, val)
    end
end
function parse_mod_val(val1, val2)
    parse(Float64, val1):parse(Float64, val2)
end
function parse_mod_val(val1, val2, val3)
    parse(Float64, val1):parse(Float64, val2):parse(Float64, val3)
end
function parse_mod(str)
    mod_name, mod_range = split(str, "=")
    return (Symbol(mod_name), parse_mod_val(split(mod_range,":")...))
end
function get_sim_name(mdb::MultiDB)
    splitpath(mdb.fns[1])[end-1]
end
function get_prototype_name(mdb::MultiDB)
    splitpath(mdb.fns[2])[end-2]
end
function get_mods(mdb::MultiDB)
    read_modifications_from_data_path(joinpath(splitpath(mdb.fns[begin])[begin:end-1]...)) 
end

function equals_str(key,val)
    "$key=$val"
end
function equals_strs(assocs)
    [equals_str(p...) for p in pairs(assocs)]
end
function mods_filename(; kwargs...)
    join(equals_strs(kwargs), "_")
end


# Data loading functions


struct DBRowIter
    db
end
Base.length(it::DBRowIter) = length(it.db)
function Base.iterate(it::DBRowIter, ((keydb, valdb), (keydb_state, valdb_state)))
    key, keydb_state = @ifsomething iterate(keydb, (keydb_state...,))
    val, valdb_state = iterate(valdb, (valdb_state...,))
    return ((key, val), ((keydb, valdb), (keydb_state, valdb_state)))
end
function Base.iterate(it::DBRowIter)
    keydb = JuliaDB.select(it.db, JuliaDB.Keys())
    valdb = JuliaDB.select(it.db, JuliaDB.Not(JuliaDB.Keys()))
    key, keydb_state = @ifsomething iterate(keydb)
    val, valdb_state = iterate(valdb)
    return ((key, val), ((keydb, valdb), (keydb_state, valdb_state)))
end

struct DBExecIter
    prototype
    db
    constant_mods
end
Base.length(it::DBExecIter) = length(it.db)
function Base.iterate(it::DBExecIter, (db_iter, db_state))
    (key, val), db_state = @ifsomething iterate(db_iter, db_state)
    model = it.prototype(; key..., it.constant_mods...)
    exec = Execution(model, BareSolution(; val...))
    return ((key, exec), (db_iter, db_state))
end
function Base.iterate(it::DBExecIter)
    db_iter = DBRowIter(it.db)
    (key, val), db_state = @ifsomething iterate(db_iter)
    model = it.prototype(; key..., it.constant_mods...)
    exec = Execution(model, BareSolution(; val...))
    return ((key, exec), (db_iter, db_state))
end

struct MultiDBRowIter
    mdb
end
function Base.iterate(it::MultiDBRowIter)
    db, mdb_state = @ifsomething iterate(it.mdb)
    db_exec_iter = DBRowIter(db)
    return iterate(it, (mdb_state, (db_exec_iter, ())))
end
function Base.iterate(it::MultiDBRowIter, (mdb_state, (db_exec_iter, wrapped_db_state)))
    row_tuple = iterate(db_exec_iter, wrapped_db_state...)
    while row_tuple === nothing
        db, mdb_state = @ifsomething iterate(it.mdb, mdb_state)
        db_exec_iter = DBRowIter(db)
        row_tuple = iterate(db_exec_iter)
    end
    row, db_state = row_tuple
    return (row, (mdb_state, (db_exec_iter, (db_state,))))
end

struct MultiDBExecIter
    prototype
    dbs::MultiDB
    constant_mods
end
function Base.iterate(it::MultiDBExecIter, (dbs_state, (db_iter, wrapped_db_state)))
    exec_tuple = iterate(db_iter, wrapped_db_state...)
    while exec_tuple === nothing
        db, dbs_state = @ifsomething iterate(it.dbs, dbs_state)
        db_iter = DBExecIter(it.prototype, db, it.constant_mods)
        exec_tuple = iterate(db_iter)
    end
    mod_exec, db_state = exec_tuple
    return (mod_exec, (dbs_state, (db_iter, (db_state,))))
end
function Base.iterate(it::MultiDBExecIter)
    (db, dbs_state) = @ifsomething iterate(it.dbs)
    db_iter = DBExecIter(it.prototype, db, it.constant_mods)
    return iterate(it, (dbs_state, (db_iter,())))
end
    
function mod_idx(particular_keys, particular_values, range_keys, range_values)
    idxs = Array{Int}(undef, length(particular_keys))
    for (pk, pv) in zip(particular_keys, particular_values)
        keydx = findfirst(range_keys .== pk)
        if keydx === nothing
            @show range_keys
            @show pk
        end
        idxs[keydx] = findfirst(range_values[keydx] .== pv)
    end
    return CartesianIndex(idxs...)
end

function get_recent_simulation_data_paths(prototype_path, max_n)
    sorted_mdb_subpaths = get_sorted_mdb_subpaths(prototype_path)
    return join(sorted_mdb_subpaths[end:-1:end-max_n+1], "\n")
end

function get_recent_simulation_data_path(prototype_path, nth::Int=0)
    sorted_mdb_subpaths = get_sorted_mdb_subpaths(prototype_path)
    if isempty(sorted_mdb_subpaths)
        error("no valid recent simulation data")
    end
    mdb_path = if nth > 0
        sorted_mdb_subpaths[nth]
    else
        sorted_mdb_subpaths[end+nth]
    end
    return mdb_path
end

function load_simulation_data(path)
    return MultiDB(path)
end

"Load nth simulation, ordered by time"
function load_simulation_data_recent(data_root, prototype_name, nth::Int=0)
    prototype_path = joinpath(data_root, prototype_name)
    sim_path = get_recent_simulation_data_path(prototype_path, nth) 
    return load_simulation_data(sim_path)
end


function load_ExecutionClassifications_recent(type::Type, prototype_name, offset_from_current=0; data_root = datadir(), kwargs...)
    recent_path = get_recent_simulation_data_path(joinpath(data_root, prototype_name), offset_from_current)
    #sim_name = splitpath(recent_path)[end] # SIM NAME UNHELPFUL FIXME
    return (load_ExecutionClassifications(type, recent_path; kwargs...), nothing) 
end

_is_varying(x::Nothing) = false
_is_varying(x::AbstractArray) = true
_is_varying(x::Number) = false

# FIXME type can be either <:AbstractArray or DataFrame -- should impact return type
# not sure why it was necessary... never actually implemented
function load_ExecutionClassifications(type::Type, data_path)
    mdb = load_simulation_data(data_path)
    mods = get_mods(mdb)

    fixed_names = [name for name in keys(mods) if !_is_varying(mods[name])]
    varying_names = [name for name in keys(mods) if _is_varying(mods[name])]
    fixed_mods = NamedTuple{Tuple(fixed_names)}([mods[key] for key in fixed_names])
    varying_mods = NamedTuple{Tuple(varying_names)}([mods[key] for key in varying_names])

    init_mod_array(T) = NamedAxisArray{keys(varying_mods)}(Array{Union{Bool,Missing}}(undef, length.(values(varying_mods))...), values(varying_mods))
    #mod_dict = Dict(name => val for (name, val) in zip(mod_names, mod_values))
    
    
    classification_names = fieldnames(ExecutionClassifications)
    classifications_A = Dict(name => init_mod_array(Bool) for name in classification_names) 

    for (this_mod, this_result) in MultiDBRowIter(mdb)
        exec_classification = this_result[:wave_properties]
        #A_idx = this_mod[varying_names]
        A_idx = NamedTuple{Tuple(varying_names)}(this_mod[varying_names])
        if exec_classification === missing
            for name in classification_names
                setindex!(classifications_A[name], missing; A_idx...)
            end
        else
            for name in classification_names
                setindex!(classifications_A[name], getproperty(exec_classification, name); A_idx...)
            end
        end
    end
    return classifications_A
end
