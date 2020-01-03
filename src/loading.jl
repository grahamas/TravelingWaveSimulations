using JuliaDB

macro ifsomething(ex)
    quote
        result = $(esc(ex))
        result === nothing && return nothing
        result
    end
end

function _load_data(sim_path)
    @warn "loading"
   JuliaDB.load(sim_path)
end
struct MultiDB
    fns
end
Base.length(mdb::MultiDB) = length(mdb.fns)

function Base.iterate(mdb::MultiDB, fns_state...)
    (next_fn, new_state) = @ifsomething iterate(mdb.fns, fns_state...)
    (_load_data(next_fn), new_state)
end

function parse_mod_val(val)
    parse(Float64, val)
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
function get_sim_name(fn::String)
    splitpath(fn)[end-1]
end
function get_example_name(fn::String)
    splitpath(fn)[end-2]
end
function get_mods(fn::String)
    sim_name = get_sim_name(fn)
    @warn "Working on $sim_name..."
    mods_str = split(sim_name, "_20")[1]
    mods_str_arr = split(mods_str, ";")
    mod_tuples = map(parse_mod, mods_str_arr)
    mod_dict = Dict(map((x) -> Pair(x...), mod_tuples)...)   
end
function get_mods(mdb::MultiDB)
    @assert all(map((path) -> splitpath(path)[end-1], mdb.fns) .== splitpath(mdb.fns[1])[end-1])
    get_mods(mdb.fns[1])
end



# Data loading functions

export DBExecIter, MultiDBExecIter, MultiDB, load_data

struct DBExecIter
    example
    db
    constant_mods
end
Base.length(it::DBExecIter) = length(it.db)
function Base.iterate(it::DBExecIter, ((keydb, valdb), (keydb_state, valdb_state)))
    key, keydb_state = @ifsomething iterate(keydb, (keydb_state...,))
    val, valdb_state = iterate(valdb, (valdb_state...,))
    model = it.example(; key..., it.constant_mods...)
    exec = Execution(model, BareSolution(; val...))
    return ((key, exec), ((keydb, valdb), (keydb_state, valdb_state)))
end
function Base.iterate(it::DBExecIter)
    keydb = select(it.db, Keys())
    valdb = select(it.db, Not(Keys()))
    key, keydb_state = @ifsomething iterate(keydb)
    val, valdb_state = iterate(valdb)
    model = it.example(; key..., it.constant_mods...)
    exec = Execution(model, BareSolution(; val...))
    return ((key, exec), ((keydb, valdb), (keydb_state, valdb_state)))
end
    

struct MultiDBExecIter
    example
    dbs::MultiDB
    constant_mods
end
function Base.iterate(it::MultiDBExecIter, (dbs_state, (db_iter, wrapped_db_state)))
    exec_tuple = iterate(db_iter, wrapped_db_state...)
    while exec_tuple === nothing
        db, dbs_state = @ifsomething iterate(it.dbs, dbs_state)
        db_iter = DBExecIter(it.example, db, it.constant_mods)
        exec_tuple = iterate(db_iter)
    end
    exec, db_state = exec_tuple
    return (exec, (dbs_state, (db_iter, (db_state,))))
end
function Base.iterate(it::MultiDBExecIter)
    (db, dbs_state) = @ifsomething iterate(it.dbs)
    db_iter = DBExecIter(it.example, db, it.constant_mods)
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
    sim_path = if nth > 0
        sorted_sims[nth]
    else
        sorted_sims[end+nth]
    end
    return (get_example(example_name), MultiDB(joinpath.(Ref(sim_path), readdir(sim_path))))
end
