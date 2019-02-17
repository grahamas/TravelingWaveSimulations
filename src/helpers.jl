
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
