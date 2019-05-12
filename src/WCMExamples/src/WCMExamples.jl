module WCMExamples

using WilsonCowanModel
using Simulation73
using MacroTools
import DifferentialEquations: Euler

function get_example(example_name)
    examples = Dict(
                    "neuman_line" => neuman_line,
                    "neuman_square" => neuman_square
                    )
    return examples[example_name]
end

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

@EI_kw_example function neuman_line(N=1, P=2)
    simulation = Simulation(;
      model = WCMSpatial{Float64,N,P}(;
        pop_names = ["E", "I"],
        α = [1.1, 1.0],
        β = [1.1, 1.1],
        τ = [0.1, 0.18],
        space = Pops{P}(Segment{Float64}(; n_points=301, extent=100.0)),
        nonlinearity = pops(SigmoidNonlinearity{Float64};
          a = [1.2, 1.0],
          θ = [2.6, 8.0]),
        stimulus = pops(NoisyStimulus{Float64,N};
          strength=[1.2, 1.2],
          time_window=[(0.0, 0.55), (0.0, 0.55)],
          width=[2.81, 2.81],
          SNR=[80.0, 80.0],
          stim_type=[SharpBumpStimulus{Float64,N}, SharpBumpStimulus{Float64,N}]),
        connectivity = pops(ShollConnectivity{Float64};
          amplitude = [16.0 -18.2;
                       27.0 -4.0],
          spread = [2.5 2.7;
                    2.7 2.5])
        ),
      solver = Solver{Float64}(;
        stop_time = 1.8,
        dt = 0.01,
        space_save_every=1,
        time_save_every=1,
        #stiffness=:stiff
        algorithm=Euler()
      )
    )
    return simulation
end

@EI_kw_example function neuman_square(N=2,P=2)
  simulation = Simulation(;
    model = WCMSpatial{Float64,N,P}(;
      pop_names = ["E", "I"],
      α = [1.1, 1.0],
      β = [1.1, 1.1],
      τ = [0.1, 0.18],
      space = Pops{P}(Grid{Float64}(; n_points=(51,51),
        extent=(50.0,50.0))),
      nonlinearity = pops(SigmoidNonlinearity{Float64};
        a = [1.2, 1.0],
        θ = [2.6, 8.0]),
        # a=Float64[BV(1.2, (0.1, 2.0)), BV(1.0, (0.1, 2.0))],
        # θ=Float64[BV(8.0, (2.0, 9.0)), BV(2.6, (2.0,9.0))]),
      stimulus = pops(NoisyStimulus{Float64,N};
          strength = [1.2, 1.2],
          width = [2.81, 2.81],
          SNR = [80.0, 80.0],
          time_window = [(0.0, 0.55), (0.0, 0.55)],
          stim_type=[SharpBumpStimulus{Float64,N}, SharpBumpStimulus{Float64,N}]),
      connectivity = pops(GaussianConnectivity{Float64};
          amplitude = [16.0 -18.2;
                       27.0 -4.0],
          spread = [(2.5,2.5) (2.7,2.7);
                    (2.7,2.7) (2.5,2.5)])
      ),
    solver = Solver{Float64}(;
        stop_time = 1.8,
        dt = 0.01,
        space_save_every=1,
        time_save_every=1,
        algorithm=Euler()
      )
  )
end

end
