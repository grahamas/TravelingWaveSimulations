
plot_specs_dict = Dict()
path_to_here = "scripts"
plot_specs_path = "plot_specs"
top_level = walkdir(joinpath(path_to_here, plot_specs_path)) |> first
plot_spec_basenames = filter(top_level[3]) do name
    name[1] != '.'
end
plot_spec_names = [splitext(basename)[1] for basename in plot_spec_basenames]
plot_spec_paths = joinpath.(plot_specs_path, plot_spec_basenames)
for (name, path) in zip(plot_spec_names, plot_spec_paths)
    include(path)
    @show plot_spec
    plot_specs_dict[name] = plot_spec
end

function get_plot_spec(plot_spec_name)
    return plot_specs_dict[plot_spec_name]
end
