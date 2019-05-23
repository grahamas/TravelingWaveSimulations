"""
Parse a function with vector kwargs into function taking corresponding EI scalar arguments
"""
function EI(fn_def)
  fn_dct = MacroTools.splitdef(fn_def)
  new_kwargs = map(fn_dct[:kwargs]) do kwarg
    EI_parsekwarg(MacroTools.splitarg(kwarg)...) .|> (x) -> MacroTools.combinearg(x...)
  end
  fn_dct[:kwargs] = cat(new_kwargs..., dims=1)
  MacroTools.combinedef(fn_dct)
end

function EI_parsekwarg(name::Symbol, type::Symbol, splat::Bool, default::Expr)
  if @capture(default, [E_, I_])
    [
      (Symbol(name,:E), :Any, false, E),
      (Symbol(name,:I), :Any, false, I),
      (name, type, splat, :([$(Symbol(name,:E)), $(Symbol(name,:I))]))
    ]
  elseif @capture(default, [EE_ EI_; IE_ II_])
    [
      (Symbol(name,:EE), :Any, false, EE),
      (Symbol(name,:IE), :Any, false, IE),
      (Symbol(name,:EI), :Any, false, EI),
      (Symbol(name,:II), :Any, false, II),
      (name, type, splat, :([$(Symbol(name,:EE)) $(Symbol(name,:EI)); $(Symbol(name,:IE)) $(Symbol(name,:II))]))
    ]
  else
    [(name, type, splat, default)]
  end
end

function EI_parsekwarg(name::Symbol, type::Symbol, splat::Bool, default::Any)
  return [(name, type, splat, default)]
end

"""Makes a top-level kwarg for every kwarg in entire function body."""
function kw_example(fn_expr)
  fn_dct = MacroTools.splitdef(fn_expr)
  # postwalk processes leaves first, so they can be used to set defaults in larger kwargs
  fn_dct[:body] = MacroTools.postwalk(fn_dct[:body]) do x
    if isa(x, Expr) && (x.head == :kw)
      push!(fn_dct[:kwargs], Expr(:kw, x.args[1], x.args[2]))
      Expr(:kw, x.args[1], x.args[1]) # INTENTIONALLY set default to same symbol
    else
      x
    end
  end
  MacroTools.combinedef(fn_dct)
end

macro EI_kw_example(fn_expr)
  fn_expr |> kw_example |> EI |> esc
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
