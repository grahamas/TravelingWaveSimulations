const modifications_prefix_filename = "modifications_prefix.txt"

function init_data_path(modifications; data_root, prototype_name, 
                                       experiment_name="",
                                       unique_id="$(Dates.now())_$(gitdescribe())")
    data_path = joinpath(data_root, prototype_name, experiment_name,
                               unique_id)
    mkpath(data_path)
    open(joinpath(data_path, modifications_prefix_filename), "w") do io
        println.(Ref(io), modifications)
    end
    return data_path
end

function read_modifications_from_data_dir(data_dir)
    modifications = Dict{Symbol,Any}()
    open(joinpath(data_dir, modifications_prefix_filename), "r") do io
        for line in readlines(io)
            parse_to_dict!(modifications, line)
        end
    end
    return modifications
end


function mods_to_pkeys(modifications)::Array{Symbol,1}
    pkeys = keys(modifications[1]) |> collect
    disallowed_keys = [:algorithm, :u, :x, :t, :n, :n_points, :extent, :save_idxs]
    pkeys = filter((x) -> !(x in disallowed_keys), pkeys)
    pkeys = filter((x) -> !(modifications[1][x] === nothing), pkeys)
end

function getkeys(d, keys)
    [d[key] for key in keys]
end

function reduce_to_namedtuple(exec::AbstractExecution)
    exec.simulation.global_reduction(exec)
end

#function extract_data_namedtuple(execution::Execution)
#    soln = execution.solution
#    u = soln.u
#    t = soln.t
#    x = coordinates(reduced_space(execution)) |> collect
#    return execution.simulation.global_reduction((u=u, t=t, x=x))
#end
#
#function extract_data_namedtuple(execution::ReducedExecution)
#    # xs is necessary for the algorithm
#    nt = (xs=frame_xs(execution), final_frame=execution.solution.u[end],
#          extract_data_namedtuple(execution.saved_values)...)
#    execution.simulation.global_reduction(nt)
#end
#
#function extract_data_namedtuple(execution::AugmentedExecution)
#    soln = execution.solution
#    u = soln.u
#    t = soln.t
#    x = coordinates(reduced_space(execution))
#    execution.simulation.global_reduction((us=u, ts=t, xs=x, extract_data_namedtuple(execution.saved_values)...))
#end
#
#function extract_data_namedtuple(sv::SavedValues{T,<:Array{<:Wavefront,1}}) where T
#    return (wavefronts=sv.saveval, wavefronts_t=sv.t)
#end

function extract_params_tuple(modification, pkeys)
    return NamedTuple{Tuple(pkeys)}(getkeys(modification, pkeys))
end

function init_results_tuple(pkeys::Array{Symbol}, results::NamedTuple{NAMES,TYPES}, init_length=0) where {NAMES,TYPES}
    N = length(pkeys)
    ALL_NAMES = Tuple([pkeys...,NAMES...])
    ARR_TYPES = Tuple{[Array{Float64,1} for _ in 1:N]..., [Array{Union{Missing,TYPE},1} for TYPE in TYPES.parameters]...}
    empty_arrs = [[Vector{Float64}(undef,init_length) for _ in 1:N]..., [Vector{Union{Missing,TYPE}}(undef, init_length) for TYPE in TYPES.parameters]...]
    NamedTuple{ALL_NAMES,ARR_TYPES}(empty_arrs)
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
function push_namedtuple!(tup::NamedTuple, mods::Union{Dict,NamedTuple})
    for mod in pairs(mods)
        push!(tup[mod[1]], mod[2])
    end
    return tup
end
push_namedtuple!(::Nothing, ::Nothing) = nothing
push_namedtuple!(nt::NamedTuple, ::Nothing) = nt

function setindex_namedtuple!(tup_arrs::NamedTuple, tup::NamedTuple, idx)
    for mod in pairs(tup)
        tup_arrs[mod[1]][idx] = mod[2]
    end
    return tup_arrs
end
function setindex_namedtuple!(tup_arrs::NamedTuple, arr::Vector{<:NamedTuple}, idx::AbstractArray)
    for key in keys(tup_arrs)
        tup_arr = tup_arrs[key]
        for i in 1:length(idx)
            tup_arr[i] = arr[i][key]
        end
    end
    return tup_arrs
end

function append_namedtuple_arr!(target_tup::NamedTuple, arr::Array{<:NamedTuple, 1})
    for key in keys(arr[1])
        append!(target_tup[key], [nt[key] for nt in arr])
    end
    return target_tup
end

function init_missing_data(sample_data::NamedTuple{NAMES,TYPES}) where {NAMES,TYPES}
     NamedTuple{NAMES}([missing for _ in NAMES])
end

function make_path_windows_safe(path)
	@assert ';' ∉ path "Can't have ';' in $path; Reserved for replacing ':'"
	# Colons are not allowed, but I need them, so I'll replace them with ;
	path = replace(path, ":" => "#")
    windows_bad_chars = Vector{Char}("""<>":|?*""")
	@assert all(windows_bad_chars .∉ path) "Your path is not Windows safe: $path"
    return path
end

function full_name(name; path="", prefix="", sep="_")
	if prefix != ""
		name = join([prefix, name], sep)
	end
	full_path = joinpath(path, name)
	@assert length(full_path) < 4096
	return full_path
end
