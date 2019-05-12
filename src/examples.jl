function get_example(example_name)
    examples = Dict(
                    "neuman_line" => neuman_line,
                    "neuman_square" => neuman_square
                    )
    return examples[example_name]
end

include("example_helpers.jl")

examples_path = "examples"
include(joinpath(examples_path, "neuman_line.jl"))
include(joinpath(examples_path, "neuman_square.jl"))
