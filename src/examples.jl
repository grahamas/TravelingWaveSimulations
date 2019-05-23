function get_example(example_name)
    examples = Dict(
                    "neuman_line" => neuman_line,
                    "neuman_square" => neuman_square,
                    "meijer_line" => meijer_line
                    )
    return examples[example_name]
end

include("example_helpers.jl")

path_to_here = "src"
examples_path = "examples"
top_level = walkdir(joinpath(path_to_here, examples_path)) |> first
examples = joinpath.(examples_path, top_level[3])
include.(examples)
