struct Value{T_LOC,T_VAL}
    loc::T_LOC
    val::T_VAL
end
Value(loc::Int, val::T) where T = (@assert loc > 0; Value{Int,T}(loc,val))
Base.:(==)(a::Value, b::Value) = a.val == b.val
Base.:(==)(a::Value, b::Number) = a.val == b
Base.isless(a::Value, b::Value) = a.val < b.val
Base.isless(a::Value, b::Number) = a.val < b
Base.isless(a::Number, b::Value) = a < b.val

# function Base.minimum(arr::Array{<:Value})
#     val, idx = findmin([elt.val for elt in arr]) 
#     return arr[idx]
# end

function linear_interpolate((x1,x2), (y1,y2), x)
    @assert x1 <= x <= x2 || x2 <= x <= x1
    if x1 == x2
        return y1 + (y2 - y1) / 2 # FIXME: if flat, just split the difference
    end
    y1 + (x - x1) * ((y2 - y1) / (x2 - x1)) 
end


function linear_interpolate_loc((left,right)::Tuple{Value,Value}, val)
    loc = linear_interpolate((left.val,right.val), (left.loc,right.loc), val)
    Value(loc, val)
end
function linear_interpolate_val((left,right)::Tuple{Value,Value}, loc)
    val = linear_interpolate((left.loc,right.loc), (left.val,right.val), loc)
    Value(loc, val)
end

all_but_last(itr, n=1) = Iterators.take(itr, length(itr)-n)

struct ValuedSpace{T,C,N,AT<:AbstractArray{T,N},AC<:AbstractArray{C,N}} <: AbstractArray{T,N}
    values::AT
    coordinates::AC
    function ValuedSpace(values::AT, coordinates::AC) where {T,C,N,AT<:AbstractArray{T,N},AC<:AbstractArray{C,N}}
        if !(size(values) == size(coordinates))
            error("Values and coordinates must match")
        end
        return new{T,C,N,AT,AC}(values, coordinates)
    end
end


Base.BroadcastStyle(::Type{<:ValuedSpace}) = Broadcast.ArrayStyle{ValuedSpace}()

Base.size(vs::ValuedSpace) = size(vs.values)
Base.getindex(vs::ValuedSpace, i::Int) = vs.values[i]
Base.getindex(vs::ValuedSpace, inds::Vararg{Int,N}) where N = @view vs.values[inds...]#ValuedSpace(vs.values[I], vs.coordinates[I])
getvalue(vs::ValuedSpace, idx::Union{CartesianIndex,Int}) = Value(vs.coordinates[idx], vs.values[idx])
#Base.getindex(vs::ValuedSpace, idx::Union{CartesianIndex,Int}) = ValuedSpace([vs.values[idx]], [vs.coordinates[idx]])
function Base.getindex(vs::ValuedSpace, fidx::AbstractFloat)
    lidx = findlast(vs.coordinates .< fidx)
    if lidx == length(vs)
        return getvalue(vs,lidx)
    end
    val = linear_interpolate_val((getvalue(vs,lidx), getvalue(vs,lidx+1)), fidx)
    return val
end
function getslice(vs::ValuedSpace, (left,right)::Tuple{Value,Value})
    return getslice(vs, (left.loc, right.loc))
end
function getslice(vs::ValuedSpace, (left_loc,right_loc)::Tuple{AbstractFloat,AbstractFloat})
    xs = vs.coordinates
    maybe_left = findlast(xs .< left_loc)
    left_idx = if maybe_left === nothing
        findfirst(xs .>= left_loc)
    else
        maybe_left + 1
    end
    maybe_right = findfirst(xs .> right_loc)
    right_idx = if maybe_right === nothing
        findlast(xs .<= right_loc)
    else
        maybe_right - 1
    end
    return @views ValuedSpace(vs.values[left_idx:right_idx], vs.coordinates[left_idx:right_idx])
end
Base.setindex!(vs::ValuedSpace, val, inds::Vararg{Int,N}) where N = vs.values[inds...] = val
Base.showarg(io::IO, vs::ValuedSpace, toplevel) = print(io, typeof(vs), ":\n\tValues: $(vs.values)\n\tCoordinates: $(vs.coordinates)")

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{ValuedSpace}}, ::Type{ElType}) where {ElType}
    # Scan the inputs for the ValuedSpace:
    vs = find_valuedspace(bc)
    # Use the char field of A to create the output
    ValuedSpace(similar(Array{ElType}, axes(bc)), vs.coordinates)
end

"`A = find_valuedspace(bc)` returns the first ValuedSpace among the arguments."
find_valuedspace(bc::Base.Broadcast.Broadcasted) = find_valuedspace(bc.args)
find_valuedspace(args::Tuple) = find_valuedspace(find_valuedspace(args[1]), Base.tail(args))
find_valuedspace(x) = x
find_valuedspace(a::ValuedSpace, rest) = a
find_valuedspace(::Any, rest) = find_valuedspace(rest)


function zero_crossing(A::ValuedSpace, i, j, atol=1e-4)
    if (A[i] <= 0 && A[j] >= 0) || (A[i] >= 0 && A[j] <= 0)
        return linear_interpolate_loc((getvalue(A,i),getvalue(A,j)), 0.0)
    elseif (A[i] < -atol && A[j] >= -atol) || (A[i] > atol && A[j] <= atol)
        return getvalue(A,j)
    else
        return nothing
    end
end

Base.diff(vs::ValuedSpace) = ValuedSpace(diff(vs.values) ./ diff(vs.coordinates), collect(all_but_last(vs.coordinates)) .+ (0.5 .* diff(vs.coordinates)))

@noinline function translate(val::Value, vs::ValuedSpace)
    vs[val.loc]
end
translate(::Nothing, args...) = nothing
