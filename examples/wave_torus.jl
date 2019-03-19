# ENV["GKSwstype"] = "100" # For headless plotting (on server)
# ENV["MPLBACKEND"]="Agg"

using Simulation73
using WilsonCowanModel
using DifferentialEquations: Euler
using JLD2
using Plots; pyplot()
using Random

Random.seed!(9999)

if !@isdefined(v)
  const v = Float64
  const BV = BoundedVariable
end
N=2
P=2
simulation = Simulation(;
  model = WCMSpatial{v,N,P}(;
    pop_names = ["E", "I"],
    α = v[1.1, 1.0],
    β = v[1.1, 1.1],
    τ = v[0.1, 0.18],
    space = Pops{P}(Torus{v}(; n_points=(51,51), extent=(50.0,50.0))),
    nonlinearity = pops(SigmoidNonlinearity{v};
      a = v[1.2, 1.0],
      θ = v[2.6, 8.0]),
      # a=v[BV(1.2, (0.1, 2.0)), BV(1.0, (0.1, 2.0))],
      # θ=v[BV(8.0, (2.0, 9.0)), BV(2.6, (2.0,9.0))]),
    # stimulus = pops(NoisyStimulus{v};
    #   strength=v[1.2, 1.2],
    #   window=Tuple{v,v}[(0.0,0.55), (0.0,0.55)],
    #   width=v[0.2, 0.2],
    #   SNR=v[80.0, 80.0],
    #   stim_type=[Sech2BumpStimulus{v}, Sech2BumpStimulus{v}]),
    stimulus = [NoisyStimulus{v,N}(;
      strength=1.2,
      time_window=Tuple{v,v}((0.0,0.55)),
      width=2.81,
      SNR= 80.0,
      stim_type=SharpBumpStimulus{v,N}),
      NoisyStimulus{v,N}(;
        strength=1.2,
        time_window=Tuple{v,v}((0.0,0.55)),
        width=2.81,
        SNR= 80.0,
        stim_type=SharpBumpStimulus{v,N})
    ],
    connectivity = pops(GaussianConnectivity{v};
      amplitude = v[16.0 -18.2
                    27.0 -4.0],
      # spread = v[BV(2.5, (1.0, 4.0)) BV(2.7, (1.0, 4.0));
      #            BV(2.7, (1.0, 4.0)) BV(2.5, (1.0, 4.0))])
      spread = Tuple{v,v}[(2.5,2.5) (2.7,2.7);
                 (2.7,2.7) (2.5,2.5)])
    ),
  solver = Solver{v}(;
    stop_time = 2.0,
    dt = 0.01,
    space_save_every=1,
    time_save_every=1,
    #stiffness=:stiff
    algorithm=Euler()
    )
  )

analyses = [
Animate(;
  fps = 20
  ),
NonlinearityPlot(;
  fn_bounds = (-1,15)
  )
# SpaceTimePlot(),
# SubsampledPlot(
#   plot_type=WaveStatsPlot,
#   time_subsampler=Subsampler(
#     Δ = 0.01,
#     window = (1.2, 1.8)
#   ),
#   space_subsampler=Subsampler(
#       window = ((0.0, 0.0),(Inf,Inf))
#     )
#   )
]

output = SingleOutput(;
  root = joinpath(homedir(), "simulation-results"),
  simulation_name = "neuman/torus"
  )

@save "parameters.jld2" simulation

execution = Execution(simulation)
analyse.(analyses, Ref(execution), Ref(output))

output(((name, obj) -> @save name obj), "parameters.jld2", simulation)
