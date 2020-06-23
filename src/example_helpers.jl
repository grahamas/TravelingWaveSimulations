"""
Parse a function with vector kwargs into function taking corresponding EI scalar arguments
"""
function EI(fn_def)
  fn_dct = MacroTools.splitdef(fn_def)
  new_kwargs = map(fn_dct[:kwargs]) do kwarg
    EI_parsekwarg(MacroTools.splitarg(kwarg)...) .|> (x) -> MacroTools.combinearg(x...)
  end
  fn_dct[:kwargs] = cat(new_kwargs..., dims=1)
  (fn_dct[:name], MacroTools.combinedef(fn_dct))
end

function parse_EIarr_arg(arr, name)
  if @capture(arr, [E_, I_])
    replaced = :([$(Symbol(name,:E)), $(Symbol(name,:I))])
    return (replaced, [
      (Symbol(name,:E), :Any, false, E),
      (Symbol(name,:I), :Any, false, I)
    ])
  elseif @capture(arr, [EE_ EI_; IE_ II_])
    replaced = :([$(Symbol(name,:EE)) $(Symbol(name,:EI)); $(Symbol(name,:IE)) $(Symbol(name,:II))])
    return (replaced, [
      (Symbol(name,:EE), :Any, false, EE),
      (Symbol(name,:IE), :Any, false, IE),
      (Symbol(name,:EI), :Any, false, EI),
      (Symbol(name,:II), :Any, false, II)
    ])
  else
    return (nothing, [])
  end
end


function EI_parsekwarg(name::Symbol, type::Symbol, splat::Bool, default::Expr)
  kwarg_list = []
  MacroTools.prewalk(default) do x
    replacement, entries = parse_EIarr_arg(x, name)
    if length(entries) > 0
      append!(kwarg_list, entries)
      append!(kwarg_list, [(name, type, splat, replacement)])
      return nothing
    elseif @capture(x, op_(arr_, other_))
      replacement, entries = parse_EIarr_arg(arr, name)
      if length(entries) > 0
        append!(kwarg_list, entries)
        append!(kwarg_list, [(name, type, splat, :($op($replacement, $other)))])
        return nothing
      else
        return x
      end
    else
      return x
    end
  end
  if length(kwarg_list) == 0
    kwarg_list = [(name, type, splat, default)]
  end
  return kwarg_list
end

function EI_parsekwarg(name::Symbol, type::Symbol, splat::Bool, default::Any)
  return [(name, type, splat, default)]
end

"""Makes a top-level kwarg for every kwarg in entire function body."""
function kw_example(fn_expr)
  fn_dct = MacroTools.splitdef(fn_expr)
  fn_dct[:name] = gensym(fn_dct[:name])
  # postwalk processes leaves first, so they can be used to set defaults in larger kwargs
  duplicated_kwarg_names = []
  seen_kwarg_names = []
  # find all duplicated kwarg names to exclude them from the function kwargs signature
  MacroTools.postwalk(fn_dct[:body]) do x
    if isa(x, Expr) && (x.head == :kw)
      kwarg_name = x.args[1]
      if kwarg_name ∈ seen_kwarg_names
        @warn "Duplicated: $kwarg_name"
        push!(duplicated_kwarg_names, kwarg_name)
      else
        push!(seen_kwarg_names, kwarg_name)
      end
    end
    x
  end
  fn_dct[:body] = MacroTools.postwalk(fn_dct[:body]) do x
    if isa(x, Expr) && (x.head == :kw) && x.args[1] ∉ duplicated_kwarg_names
      push!(fn_dct[:kwargs], Expr(:kw, x.args[1], x.args[2]))
      Expr(:kw, x.args[1], x.args[1]) # INTENTIONALLY set default to same symbol
    else
      x
    end
  end
  sym = MacroTools.combinedef(fn_dct)
  @show sym
end

macro EI_kw_example(fn_expr)
  kw_fn = fn_expr |> kw_example #|> EI
  fn_dct = MacroTools.splitdef(kw_fn)
  quote
    $kw_fn
    example = $(fn_dct[:name])
  end |> esc
end

# From DrWatson
function increment_backup_num(filepath)
    path, filename = splitdir(filepath)
    fname, suffix = splitext(filename)
    m = match(r"^(.*)_#([0-9]+)$", fname)
    if m == nothing
        return joinpath(path, "$(fname)_#1$(suffix)")
    end
    newnum = string(parse(Int, m.captures[2]) +1)
    return joinpath(path, "$(m.captures[1])_#$newnum$(suffix)")
end
function recursively_clear_path(cur_path)
    isfile(cur_path) || return
    new_path=increment_backup_num(cur_path)
    if isfile(new_path)
        recursively_clear_path(new_path)
    end
    mv(cur_path, new_path)
end
