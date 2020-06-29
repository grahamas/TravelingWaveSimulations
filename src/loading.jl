using JuliaDB

macro ifsomething(ex)
    quote
        result = $(esc(ex))
        result === nothing && return nothing
        result
    end
end

function _load_data(sim_path)
   JuliaDB.load(sim_path)
end
struct MultiDB
    fns
end
Base.length(mdb::MultiDB) = length(mdb.fns)

is_valid_fn(::Nothing) = false
function is_valid_fn(fn)
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

function Base.iterate(mdb::MultiDB, fns_state...)
    next_fn = nothing
    while !is_valid_fn(next_fn)
        (next_fn, fns_state) = @ifsomething iterate(mdb.fns, fns_state...)
    end
    (_load_data(next_fn), fns_state)
end

function parse_mod_val(val)
    if val == "nothing"
        return nothing
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
function get_sim_name(fn::String)
    splitpath(fn)[end-1]
end
function get_example_name(fn::String)
    splitpath(fn)[end-2]
end
function get_mods(fn::String)
    sim_name = get_sim_name(fn)
    @warn "Working on $sim_name..."
    if sim_name[1:4] ∈ ["2019", "2020"] 
        # No mods
        return []
    end
    mods_str = split(sim_name, "_20")[1]
    mods_str_arr = split(mods_str, ";")
    mod_tuples = map(parse_mod, mods_str_arr)
    mod_dict = Dict(map((x) -> Pair(x...), mod_tuples)...)   
    return mod_dict
end
function get_mods(mdb::MultiDB)
    @assert all(map((path) -> splitpath(path)[end-1], mdb.fns) .== splitpath(mdb.fns[1])[end-1])
    get_mods(mdb.fns[1])
end



# Data loading functions

export DBRowIter, MultiDBRowIter, DBExecIter, MultiDBExecIter, MultiDB, load_data

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
    example
    db
    constant_mods
end
Base.length(it::DBExecIter) = length(it.db)
function Base.iterate(it::DBExecIter, (db_iter, db_state))
    (key, val), db_state = @ifsomething iterate(db_iter, db_state)
    model = it.example(; key..., it.constant_mods...)
    exec = Execution(model, BareSolution(; val...))
    return ((key, exec), (db_iter, db_state))
end
function Base.iterate(it::DBExecIter)
    db_iter = DBRowIter(it.db)
    (key, val), db_state = @ifsomething iterate(db_iter)
    model = it.example(; key..., it.constant_mods...)
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
    mod_exec, db_state = exec_tuple
    return (mod_exec, (dbs_state, (db_iter, (db_state,))))
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

export load_ExecutionClassifications_recent
function load_ExecutionClassifications_recent(example_name, data_root = datadir(), offset_from_current=0)
    (example, mdb) = load_data(data_root, example_name, offset_from_current)
    example_name = get_example_name(mdb.fns[1])
    sim_name = get_sim_name(mdb.fns[1])

    mods = get_mods(mdb)

    is_range(x::Nothing) = false
    is_range(x::AbstractRange) = true
    is_range(x::Number) = false
    
    all_mod_names = keys(mods) |> collect
    all_mod_values = values(mods) |> collect
    varied_mods = is_range.(all_mod_values)
    fixed_mods = .!varied_mods
    fixed_mods_dict = Dict(name => mods[name] for name in all_mod_names[fixed_mods])
    
    mod_names = all_mod_names[varied_mods] |> sort
    mod_values = [mods[name] for name in mod_names]
    mod_dict = Dict(name => val for (name, val) in zip(mod_names, mod_values))
    
    mod_array(T) = NamedAxisArray{Tuple(mod_names)}(Array{Union{Bool,Missing}}(undef, length.(mod_values)...), Tuple(mod_values))
    
    first_result = MultiDBRowIter(mdb) |> first
    classification_names = fieldnames(ExecutionClassifications)
    classifications_A = Dict(name => mod_array(Bool) for name in classification_names) 

    for (this_mod, this_result) in MultiDBRowIter(mdb)
        exec_classification = this_result[:wave_properties]
        classifications_A_idx = this_mod[mod_names]
        if exec_classification === missing
            for name in classification_names
                classifications_A[name][classifications_A_idx...] = missing
            end
        else
            for name in classification_names
                classifications_A[name][classifications_A_idx...] = getproperty(exec_classification, name)
            end
        end
    end
    return classifications_A
end
