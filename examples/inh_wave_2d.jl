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
    α = v[1.0, 1.0],
    β = v[1.0, 1.0],
    τ = v[10.0, 10.0],
    space = Pops{P}(Grid{v}(; n_points=(51,51), extent=(500.0,500.0))),
    nonlinearity = pops(GaussianNonlinearity{v};
      sd = v[6.7, 3.2],
      θ = v[18.0, 10.0]),
    # nonlinearity = pops(SigmoidNonlinearity{v};
    #   a = v[2.0, 0.95],
    #   θ = v[12.41, 7.33]),
    stimulus = [NoisyStimulus{v,N}(;
      strength=10.0,
      time_window=Tuple{v,v}((0.0,10.0)),
      width=100.0,
      mean = 1.0,
      SNR = 80.0,
      stim_type=SharpBumpStimulus{v,N}),
      NoStimulus{v,N}()
    ],
    connectivity = pops(GaussianConnectivity{v};
    amplitude = v[280.0 -297.0;
                  270.0 -1.4],
      spread = Tuple{v,v}[(70.0,70.0) (90.0,90.0);
                 (90.0,90.0) (70.0,70.0)])
    ),
  solver = Solver{v}(;
    stop_time = 100.0,
    dt = 1.0,
    space_save_every=1,
    time_save_every=1,
    algorithm=Euler()
    )
  )

analyses = [
Animate(;
  fps = 15
  ),
# NonlinearityPlot(;
#   fn_bounds = (-1,15)
#   ),
# SpaceTimePlot(),
# SubsampledPlot(
#   plot_type=WaveStatsPlot,
#   time_subsampling=Dict(
#     :Δsubsampled => 0.01,
#     :scalar_window => (1.2, 1.8)
#     ),
#   space_subsampling=Dict(
#     :scalar_window => (5.0,Inf)
#     )
#   )
]

output = SingleOutput(;
  root = joinpath(homedir(), "simulation-results"),
  simulation_name = "inh_wave_2d"
  )

execution = Execution(simulation)
analyse.(analyses, Ref(execution), Ref(output))

output(((name, obj) -> @save name obj), "parameters.jld2", simulation)
