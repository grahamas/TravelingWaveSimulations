
#region AbstractModel Deconstruction & Reconstruction
#########################################################
######### MODEL DECONSTRUCTION & RECONSTRUCTION #########
# Deconstruct a model (which is naturally hierarchical)
# into an array representation.
# Alter array representation.
# Reconstruct hierarchical model.
#########################################################

abstract type AbstractDeconstruction{OBJ} end

struct Deconstruction{OBJ} <: AbstractDeconstruction{OBJ}
    source::OBJ
    deconstruction::Tuple{Type,Array}
end


"""
    set_deep_dx!(model_deconstruction, dxs, val)

Within an object deconstructed into an array, set one element to `val`,
where `dxs` indicate that element by successive indexing into nested arrays.
"""
function set_deep_dx!(model_deconstruction::Tuple{Type,Array}, dxs, val)
    tmp = model_deconstruction
    for dx in dxs[1:end-1]
        tmp = tmp[2][dx]
    end
    tmp[2][dxs[end]] = val
end


abstract type AbstractDeconstructor end
struct Deconstructor <: AbstractDeconstructor end


function deconstruct(val::Union{AbstractString,Number}, ad::AbstractDeconstructor=Deconstructor())
    return (typeof(val), val)
end

function deconstruct(arr::AA, ad::AbstractDeconstructor=Deconstructor()) where {AA<:AbstractArray}
    deconstruction = map(ad, arr)
    return (base_type(typeof(arr)), [deconstruction...]) # remade incase static
end

function deconstruct(tup::Tuple, ad::AbstractDeconstructor=Deconstructor())
    deconstruction = map(ad, tup)
    return (base_type(typeof(tup)), [deconstruction...])
end

function base_type(T::Int)
    T
end


function base_type(::Type{T}) where T <: Real
    return T
end

function base_type(::Type{Array{T,N}}) where {T, N}
    BT = base_type(T)
    return Array{BT,N}
end

function base_type(::Type{Tuple{T,S}}) where {T,S}
    BT = base_type(T)
    BS = base_type(S)
    return Tuple{BT,BS}
end

# function base_type(::Type{SA}) where {N,M,T,TUP,SA<:Union{SArray{TUP,T,N,M},MArray{TUP,T,N,M},SizedArray{TUP,T,N,M}}}
#     BT = base_type(T)
#     return SArray{TUP,BT,N,M}
# end


function reconstruct(tup::Tuple{Type,<:Union{Number,AbstractString}})
    return tup[2]
end

function reconstruct(tup::Tuple{Type,Array})
    typ, arr = tup
    base_typ = base_type(typ)
    if typ <: Union{AbstractArray, Tuple}
        return base_typ(reconstruct.(arr))
    else
        return base_typ(reconstruct.(arr)...)
    end
end

############################################
##endregion
