function get_example(example_name)
    examples = Dict(
                    "neuman_line" => neuman_line,
                    "neuman_square" => neuman_square,
                    "meijer_line" => meijer_line
                    )
    return examples[example_name]
end

include("example_helpers.jl")

examples_path = "examples"
include(joinpath(examples_path, "neuman_line.jl"))
include(joinpath(examples_path, "neuman_square.jl"))
include(joinpath(examples_path, "meijer_line.jl"))
