export parse_modifications_argument, parse_analyses_argument

MOD_SEP = "+"

parse_num(str) = try
        parse(Int, str)
    catch
        parse(Float64, str)
    end

parse_val(val_str) = eval(Meta.parse(val_str))

parse_range(start, stop) = parse_num(start):parse_num(stop)
parse_range(start, step, stop) = parse_num(start):parse_num(step):parse_num(stop)

parse_array(first::T, args...) where T = T[first, args...]
parse_array(array_strs::AbstractArray{<:AbstractString}) = parse_array(parse_val.(array_strs)...)
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
            return Dict(Symbol(name_str) => parse_val(value_str))
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
        modifications_prefix = make_prefix(modification_strs)
    else
        parsed_modifications = [Dict()]
        modifications_prefix = ""
    end
    return parsed_modifications, modifications_prefix
end

function parse_analyses_array(analysis_strs::AbstractArray)
    parsed_analyses = @> analysis_strs begin
        parse_analysis.()
        must_be_list.()
    end
    return cat(parsed_analyses..., dims=1)
end

parse_analysis(analysis_str) = get_analysis(analysis_str)

function parse_analyses_argument(analysis_strs)
    if length(analysis_strs) == 0
        return [], ""
    end
    parsed_analyses = parse_analyses_array(analysis_strs)
    analyses_prefix = make_prefix(analysis_strs)
    return parsed_analyses, analyses_prefix
end
