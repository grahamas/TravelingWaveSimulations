export parse_modifications_argument, parse_plot_specs_argument

parse_range(start, stop) = parse(Float64, start):parse(Float64, stop)
parse_range(start, step, stop) = parse(Float64, start):parse(Float64, step):parse(Float64, stop)

parse_array(first::T, args...) where T = T[first, args...]
parse_array(array_strs::AbstractArray{<:AbstractString}) = parse_array(parse.(Float64, array_strs)...)
parse_array(array_str::AbstractString) = @> array_str strip(['[', ']']) split(',') parse_array

function parse_modification(str::AbstractString)
    if occursin("=", str)
        name_str, value_str = split(str, "=")
        if value_str[1] == '['
            @assert value_str[end] == ']'
            return [Dict(Symbol(name_str) => val) for val in parse_array(value_str)]
        end
        if occursin(":", value_str)
            return [Dict(Symbol(name_str) => val) for val in parse_range(split(value_str,":")...)]
        else
            return Dict(Symbol(name_str) => parse(Float64, value_str))
        end
    else
        return get_modification(str)
    end
end

must_be_list(x::AbstractArray) = x
must_be_list(x) = [x]

function make_prefix(strs, path_prefix="")
    prefix = "$(join(sort(strs), ';'))_"
    if length(prefix) >= 255
        return make_prefix(strs[2:end], path_prefix=joinpath(path_prefix, strs[1]))
    else
        return joinpath(path_prefix,prefix)
    end
end

function parse_modifications_array(modification_strs::AbstractArray)
    parsed_modifications = @> modification_strs begin
        parse_modification.()
        must_be_list.()
    end
    modification_cases = Iterators.product(parsed_modifications...)
    modification_cases = map((case) -> merge(case...), modification_cases)
    return modification_cases
end

function parse_modifications_argument(modification_strs)
    if modification_strs != []
        parsed_modifications = parse_modifications_array(modification_strs)
        modifications_prefix = make_prefix(modification_strs)
    else
        parsed_modifications = [Dict()]
        modifications_prefix = ""
    end
    return parsed_modifications, modifications_prefix
end

function parse_plot_specs_array(plot_spec_strs::AbstractArray)
    parsed_plot_specs = @> plot_spec_strs begin
        parse_plot_spec.()
        must_be_list.()
    end
    return cat(parsed_plot_specs..., dims=1)
end

parse_plot_spec(plot_spec_str) = get_plot_spec(plot_spec_str)

function parse_plot_specs_argument(plot_spec_strs)
    if length(plot_spec_strs) == 0
        return [], ""
    end
    parsed_plot_specs = parse_plot_specs_array(plot_spec_strs)
    plot_specs_prefix = make_prefix(plot_spec_strs)
    return parsed_plot_specs, plot_specs_prefix
end
