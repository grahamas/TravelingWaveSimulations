@doc """
    abs_difference(edge)

Return the difference between two points in euclidean space, given an edge between those points.

# Example
```jldoctest
julia> abs_difference( (5,1) )
4

julia> abs_difference( ((2,2), (5,-5)) )
(3, 7)
```
"""
abs_difference(edge::Tuple{T,T}) where T<:Number = abs(edge[1] - edge[2])
# FIXME should really be L2 norm
abs_difference(edge::Tuple{Tup,Tup}) where {Nminusone,T,Tup<:Tuple{T,Vararg{T,Nminusone}}} = abs.(edge[1] .- edge[2])



"""
    abs_difference_periodic(edge, period)

Return the difference between two points in euclidean space as in abs_difference, but let the space wrap with period.

# Example
```jldoctest
julia> abs_difference_periodic( (5,1), 3 )
-1

julia> abs_difference_periodic( ((2,2), (5,-5)), (3,4) )
(0, -3)
```
"""
function abs_difference_periodic(edge::Tuple{T,T}, period::T) where T<:Number
    diff = abs_difference(edge)
    if diff > period / 2
        return period - diff
    else
        return diff
    end
end
function abs_difference_periodic(edge::Tuple{Tup,Tup}, periods::Tup) where {Nminusone,T,Tup<:Tuple{T,Vararg{T,Nminusone}}}
    diffs = abs_difference(edge)
    diffs = map(zip(diffs, periods)) do (diff, period)
        if diff > period / 2
            return period - diff
        else
            return diff
        end
    end
    return Tup(diffs)
end
