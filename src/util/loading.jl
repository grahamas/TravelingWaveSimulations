
macro ifsomething(ex)
    quote
        result = $(esc(ex))
        result === nothing && return nothing
        result
    end
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
# function get_sim_name(mdb::MultiDB)
#     splitpath(mdb.fns[1])[end-1]
# end
# function get_prototype_name(mdb::MultiDB)
#     splitpath(mdb.fns[2])[end-2]
# end

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

is_valid_bson_path(::Nothing) = false
function is_valid_bson_path(fn)
    @show fn
    splitted = splitext(fn)
    if length(splitted) == 1
        return false
    elseif splitted[2] != ".bson"
        return false
    elseif basename(splitted[1]) == "failures"
        return false
    else
        return true
    end
end
function is_valid_mbson_path(path)
    any(map(is_valid_bson_path, joinpath.(Ref(path), readdir(path))))
end
function get_sorted_bson_subpaths(path)
    filter(is_valid_bson_path, joinpath.(Ref(path), readdir(path)))
end
function get_sorted_mbson_subpaths(path)
    @show path
    mbson_subpaths = filter(is_valid_mbson_path, joinpath.(Ref(path), readdir(path)))
    sorted_mbson_subpaths = sort(mbson_subpaths, by=mtime)
    return sorted_mbson_subpaths
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
    sorted_mdb_subpaths = get_sorted_bson_subpaths(prototype_path)
    return join(sorted_mdb_subpaths[end:-1:end-max_n+1], "\n")
end

function get_recent_simulation_data_path(prototype_path, nth::Int=0)
    sorted_mbson_subpaths = get_sorted_mbson_subpaths(prototype_path)
    if isempty(sorted_mbson_subpaths)
        error("no valid recent simulation data")
    end
    mdb_path = if nth > 0
        sorted_mbson_subpaths[nth]
    else
        sorted_mbson_subpaths[end+nth]
    end
    return mdb_path
end

function load_simulation_data(path)
    # `only` because currently only support single bson file
    bson_path = get_sorted_bson_subpaths(path) |> only
    BSON.@load bson_path table
    mods = read_modifications_file(path)
    return (table, mods)
end

"Load nth simulation, ordered by time"
function load_simulation_data_recent(data_root, prototype_name, nth::Int=0)
    prototype_path = joinpath(data_root, prototype_name)
    sim_path = get_recent_simulation_data_path(prototype_path, nth)
    return load_simulation_data(sim_path)
end


function load_classifications_recent(prototype_name, offset_from_current=0; data_root = datadir(), kwargs...)
    recent_path = get_recent_simulation_data_path(joinpath(data_root, prototype_name), offset_from_current)
    #sim_name = splitpath(recent_path)[end] # SIM NAME UNHELPFUL FIXME
    return (load_classifications(recent_path; kwargs...), nothing) 
end

_is_varying(x::Nothing) = false
_is_varying(x::AbstractArray) = true
_is_varying(x::Number) = false

function findfirsts_nt(nt_arrs::NamedTuple{NAMES}, nt_vals::NamedTuple{NAMES}) where NAMES
    NamedTuple{NAMES}(
        [findfirst(isapprox.(arr, val, atol=1e-8)) for (arr, val) ∈ zip(nt_arrs, nt_vals)]
    )
end
# FIXME type can be either <:AbstractArray or DataFrame -- should impact return type
# not sure why it was necessary... never actually implemented
function load_classifications(data_path; classification_name=:propagation)
    data_table, mods = load_simulation_data(data_path)

    fixed_names = [name for name in keys(mods) if !_is_varying(mods[name])]
    varying_names = [name for name in keys(mods) if _is_varying(mods[name])]
    fixed_mods = NamedTuple{Tuple(fixed_names)}([mods[key] for key in fixed_names])
    varying_mods = NamedTuple{Tuple(varying_names)}([mods[key] for key in varying_names])

    init_mod_array(T) = NamedDimsArray{keys(varying_mods)}(Array{Union{Bool,Missing}}(undef, length.(values(varying_mods))...))
    #mod_dict = Dict(name => val for (name, val) in zip(mod_names, mod_values))
    
    first_row::NamedTuple{NAMES} = first(rows(data_table))
    first_result = first_row[classification_name]
    mods_from_row = (row) -> NamedTuple{NAMES[NAMES != classification_name]}(row[NAMES != classification_name])
    first_mods = mods_from_row(first_row)
    cls_sym, cls_type = get_cls_type(first_result)
    classification_names = fieldnames(cls_type)
    classifications_A = Dict(name => init_mod_array(Bool) for name in classification_names) 

    for row in rows(data_table)
        exec_classification = row[cls_sym]
        this_mod = mods_from_row(row)
        #A_idx = this_mod[varying_names]
        A_idx = findfirsts_nt(varying_mods, this_mod)
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

# function get_cls_type(example_result::NamedTuple{NAMES}) where NAMES
#     if :wave_properties ∈ NAMES
#         return (:wave_properties, ExecutionClassifications)
#     elseif :propagation ∈ NAMES
#         return (:propagation, MinimalPropagationClassification)
#     else
#         error("Unknown NamedTuple structure for execution classification")
#     end
# end
