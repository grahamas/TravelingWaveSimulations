MOD_SEP = "+"

parse_num(str) = try
        parse(Int, str)
    catch e
        parse(Float64, str)
    end

parse_val(val_str) = eval(Meta.parse(val_str))

parse_range(start, stop) = parse_num(start):parse_num(stop)
parse_range(start, step, stop) = parse_num(start):parse_num(step):parse_num(stop)
parse_range(str::AbstractString) = parse_range(split(str, ':')...)

parse_array(first::T, args...) where T = T[first, args...]
parse_array(array_strs::AbstractArray{<:AbstractString}) = parse_array(parse_val.(array_strs)...)
parse_array(array_str::AbstractString) = @> array_str strip(['[', ']']) split(',') parse_array

function parse_rhs(value_str)
    if value_str[1] == '['
        @assert value_str[end] == '['
        return parse_array(value_str)
    elseif occursin(":", value_str)
        return parse_range(value_str)
    else
        return parse_val(value_str)
    end
end

function parse_modification(str::AbstractString)
    name_str, value_str = split(str, "=")
    return Dict(Symbol(name_str) => parse_rhs(value_str))
end

function dict_array_from_array_dict(dct)
    key = only(keys(dct))
    arr = only(values(dct))
    return [Dict(key => val) for val in arr]
end

function parse_modification_to_dict!(dict, modification_str)
    merge!(dict, parse_modification(modification_str))
end

function make_modifications_prefix(strs, path_prefix="")
    prefix = "$(join(sort(strs), ';'))"
end

function parse_modifications_array(modification_strs::AbstractArray)
    parsed_modifications = parse_modification.(modification_strs)
    arrays_of_cases = dict_array_from_array_dict.(parsed_modifications)
    modification_cases = Iterators.product(arrays_of_cases...)
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
    else
        parsed_modifications = [Dict()]
    end
    return parsed_modifications, modification_strs
end
