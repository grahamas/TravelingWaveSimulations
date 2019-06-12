
include("example_helpers.jl")

examples_dict = Dict()
path_to_here = "src"
examples_path = "examples"
top_level = walkdir(joinpath(path_to_here, examples_path)) |> first
example_basenames = filter(top_level[3]) do name
    name[1] != '.'
end
example_names = [splitext(basename)[1] for basename in example_basenames]
example_paths = joinpath.(examples_path, example_basenames)
for (name, path) in zip(example_names, example_paths)
    include(path)
    @show example
    examples_dict[name] = example
end

function get_example(example_name)
    return examples_dict[example_name]
end
