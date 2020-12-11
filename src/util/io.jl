
# FIXME move exports
# FIXME move imports
# FIXME finish refactoring io, saving, and loading

const modifications_prefix_filename = "modifications_prefix.txt"

export read_params_from_data_path, make_path_windows_safe

function read_params_from_data_path(fpath)
    sim_params = read_modifications_from_data_path(fpath)
    path_components = splitpath(fpath)
    data_dx = findfirst(path_components .== "data")
    prototype_name = path_components[data_dx + 1]
    return prototype_name, sim_params
end

function read_modifications_from_data_path(data_path)
    mod_txt_fn = joinpath(data_path, modifications_prefix_filename) 
    if isfile(mod_txt_fn)
        @info "Found modifications text file"
        modifications_dict = Dict{Symbol,Any}()
        open(mod_txt_fn, "r") do io
            for modification_line in readlines(io)
                parse_modification_to_dict!(modifications_dict, modification_line)
            end
        end
        return modifications_dict
    else
        @warn "DEPRECATED: modifications file not found; reading from filename"
        return parse_modifications_filename(basename(data_path))
    end
end

function write_modifications!(dir, mods::AbstractVector{<:AbstractString})
    open(joinpath(dir, modifications_prefix_filename), "w") do io
        println.(Ref(io), mods)
    end
end
write_modifications!(dir, mods::Union{Dict,NamedTuple}, args...) = write_modifications!(dir, [string(pair) for pair in pairs(mods)], args...)


function init_data_path(modifications; data_root, 
                prototype_name, 
                experiment_name="",
                unique_id="$(Dates.now())_$(gitdescribe())")
    data_path = joinpath(data_root, prototype_name, experiment_name, unique_id)
    mkpath(data_path)
    write_modifications!(data_path, modifications)
    return data_path
end

function make_path_windows_safe(path)
	@assert ';' ∉ path "Can't have ';' in $path; Reserved for replacing ':'"
	# Colons are not allowed, but I need them, so I'll replace them with ;
	path = replace(path, ":" => "#")
    windows_bad_chars = Vector{Char}("""<>":|?*""")
	@assert all(windows_bad_chars .∉ path) "Your path is not Windows safe: $path"
    return path
end

function parse_prototype_name_from_mdb_path(path)
    prototype_dir = if isdirpath(path)
        splitpath(path)[end-2]
    else
        splitpath(path)[end-3]
    end
    return prototype_dir
end