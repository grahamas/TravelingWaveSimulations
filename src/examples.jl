# TODO: Autogenerate (macro) this mapping
function get_example(example_name)
    examples = Dict(
                    "neuman_line" => neuman_line,
                    "neuman_square" => neuman_square,
                    "meijer_line" => meijer_line,
                    "meijer_square" => meijer_square,
                    "neuman_square_noise" => neuman_square_noise
                    )
    return examples[example_name]
end

include("example_helpers.jl")

path_to_here = "src"
examples_path = "examples"
top_level = walkdir(joinpath(path_to_here, examples_path)) |> first
example_basenames = filter(top_level[3]) do name
    name[1] != '.'
end
examples = joinpath.(examples_path, example_basenames)
include.(examples)
