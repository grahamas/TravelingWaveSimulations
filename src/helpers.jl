
# from MikeInnes
# Speed up anonymous functions 10x
# @dotimed 10^8 (x -> x^2)(rand()) # 3.9 s
# @dotimed 10^8 (@fn x -> x^2)(rand()) # 0.36 s
macro fn(expr::Expr)
  @assert expr.head in (:function, :->)
  name = gensym()
  args = expr.args[1]
  args = typeof(args) == Symbol ? [args] : args.args
  body = expr.args[2]
  @eval $name($(args...)) = $body
  name
end

using MacroTools: splitdef, combinedef, splitarg
macro memoize(args...)
    if length(args) == 1
        dicttype = :(IdDict)
        ex = args[1]
    elseif length(args) == 2
        (dicttype, ex) = args
    else
        error("Memoize accepts at most two arguments")
    end

    def_dict = try
        splitdef(ex)
    catch
        error("@memoize must be applied to a method definition")
    end
    # a return type declaration of Any is a No-op because everything is <: Any
    rettype = get(def_dict, :rtype, Any)
    f = def_dict[:name]
    def_dict_unmemoized = copy(def_dict)
    def_dict_unmemoized[:name] = u = Symbol("##",f,"_unmemoized")

    args = def_dict[:args]
    kws = def_dict[:kwargs]
    # Set up arguments for tuple
    tup = [splitarg(arg)[1] for arg in vcat(args, kws)]

    # Set up identity arguments to pass to unmemoized function
    identargs = map(args) do arg
        arg_name, typ, slurp, default = splitarg(arg)
        if slurp
            Expr(:..., arg_name)
        else
            arg_name
        end
    end
    identkws = map(kws) do kw
        arg_name, typ, slurp, default = splitarg(kw)
        if slurp
            Expr(:..., arg_name)
        else
            Expr(:kw, arg_name, arg_name)
        end
    end

    fcachename = Symbol("##",f,"_memoized_cache")
    fcache = eval(:(@isdefined $fcachename)) ?
             eval(:($fcachename)) :
             eval(:(const $fcachename = ($dicttype)()))

    if length(kws) == 0
        lookup = :($fcache[($(tup...),)]::Core.Compiler.return_type($u, typeof(($(identargs...),))))
    else
        lookup = :($fcache[($(tup...),)])
    end

    def_dict[:body] = quote
        haskey($fcache, ($(tup...),)) ? $lookup :
        ($fcache[($(tup...),)] = $u($(identargs...),; $(identkws...)))
    end
    esc(quote
        $(combinedef(def_dict_unmemoized))
        empty!($fcache)
        Base.@__doc__($(combinedef(def_dict)))
    end)

end

# Thanks to Tamas Papp
function view_slice_last(arr::AbstractArray{T,N}, dx::Int) where {T,N}
    view(arr, ntuple(_ -> Colon(), N - 1)..., dx)
end

function view_slice_last(arr::AbstractArray{T,N}, dx::CartesianIndex{DX}) where {T,N,DX}
    view(arr, ntuple(_ -> Colon(), N - DX)..., dx)
end
