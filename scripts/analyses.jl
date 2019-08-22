
analyses_dict = Dict()
analyses_path = joinpath(scriptsdir(), "analyses")
top_level = walkdir(analyses_path) |> first
analysis_basenames = filter(top_level[3]) do name
    name[1] != '.'
end
analysis_names = [splitext(basename)[1] for basename in analysis_basenames]
analysis_paths = joinpath.(analyses_path, analysis_basenames)
for (name, path) in zip(analysis_names, analysis_paths)
    include(path)
    analyses_dict[name] = analyses
end

function get_analysis(analysis_name)
    return analyses_dict[analysis_name]
end
