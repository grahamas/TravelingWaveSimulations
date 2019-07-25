export parse_modifications_argument, parse_plot_specs_argument

MOD_SEP = "+"

parse_num(str) = try
        parse(Int, str)
    catch
        parse(Float64, str)
    end

parse_range(start, stop) = parse_num(start):parse_num(stop)
parse_range(start, step, stop) = parse_num(start):parse_num(step):parse_num(stop)

parse_array(first::T, args...) where T = T[first, args...]
parse_array(array_strs::AbstractArray{<:AbstractString}) = parse_array(parse_num.(array_strs)...)
parse_array(array_str::AbstractString) = @> array_str strip(['[', ']']) split(',') parse_array

function parse_modification(str::AbstractString)
    if occursin("=", str)
        name_str, value_str = split(str, "=")
        if value_str[1] == '['
            @assert value_str[end] == ']'
            return [Dict(Symbol(name_str) => val) for val in parse_array(value_str)]
        elseif occursin(":", value_str)
            return [Dict(Symbol(name_str) => val) for val in parse_range(split(value_str,":")...)]
        else
            return Dict(Symbol(name_str) => parse_num(value_str))
        end
    else
        return get_modification(str)
    end
end

must_be_list(x::AbstractArray) = x
must_be_list(x) = [x]

function parse_modifications_array(modification_strs::AbstractArray)
    parsed_modifications = @> modification_strs begin
        parse_modification.()
        must_be_list.()
    end
    modification_cases = Iterators.product(parsed_modifications...)
    modification_cases = map(modification_cases) do cases
        any_dict = Dict{Symbol,Any}()
        merge!(any_dict, cases...)
        return any_dict
    end
    return modification_cases
end

function parse_modifications_argument(modification_strs)
    if modification_strs != []
        parsed_modifications = parse_modifications_array(modification_strs)
        modifications_prefix = """$(join(sort(modification_strs), MOD_SEP))_"""
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
    plot_specs_prefix = """$(join(sort(plot_spec_strs), MOD_SEP))"""
    return parsed_plot_specs, plot_specs_prefix
end
