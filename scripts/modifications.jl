
modifications_dict = Dict()
modifications_path = joinpath(@__DIR__, "modifications")
top_level = walkdir(modifications_path) |> first
modification_basenames = filter(top_level[3]) do name
    name[1] != '.'
end
modification_names = [splitext(basename)[1] for basename in modification_basenames]
modification_paths = joinpath.(modifications_path, modification_basenames)
for (name, path) in zip(modification_names, modification_paths)
    include(path)
    modifications_dict[name] = modifications
end

function get_modification(modification_name)
    return modifications_dict[modification_name]
end
