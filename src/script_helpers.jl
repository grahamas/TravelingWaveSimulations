script_helpers_path = joinpath(srcdir(), "script_helpers")
top_level = walkdir(script_helpers_path) |> first
script_helper_basenames = filter(top_level[3]) do name
    name[1] != '.'
end
script_helper_paths = joinpath.(script_helpers_path, script_helper_basenames)
for path in script_helper_paths
    include(path)
end
