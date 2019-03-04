# ENV["GKSwstype"] = "100" # For headless plotting (on server)
# ENV["MPLBACKEND"]="Agg"

using Simulation73
using WilsonCowanModel
using DifferentialEquations: Euler
using JLD2

if !@isdefined(v)
  const v = Float64
  const BV = BoundedVariable
end
N=1
P=2
simulation = Simulation(;
  model = WCMSpatial{v,N,P}(;
    pop_names = ["E", "I"],
    α = v[1.1, 1.0],
    β = v[1.1, 1.1],
    τ = v[0.1, 0.18],
    space = PopSegment{v,P}(; n_points=301, extent=100.0),
    nonlinearity = pops(Sech2Nonlinearity{v}; 
      a = v[1.2, 1.0],
      θ = v[2.6, 8.0]),
      # a=v[BV(1.2, (0.1, 2.0)), BV(1.0, (0.1, 2.0))],
      # θ=v[BV(8.0, (2.0, 9.0)), BV(2.6, (2.0,9.0))]),
    stimulus = pops(NoisyStimulus{v}; 
      strength=v[1.2, 1.2],
      window=Tuple{v,v}[(0.0,0.55), (0.0,0.55)],
      width=v[2.81, 2.81],
      SNR=v[80.0, 80.0],
      stim_type=[Sech2BumpStimulus{v}, Sech2BumpStimulus{v}]),
    connectivity = pops(ShollConnectivity{v};
      amplitude = v[16.0 -18.2
                    27.0 -4.0],
      # spread = v[BV(2.5, (1.0, 4.0)) BV(2.7, (1.0, 4.0));
      #            BV(2.7, (1.0, 4.0)) BV(2.5, (1.0, 4.0))])
      spread = v[2.5 2.7;
                 2.7 2.5])
    ),
  solver = Solver(;
    stop_time = 1.8,
    dt = 0.005,
    space_save_every=1,
    time_save_every=1,
    algorithm=Euler()
    )
  )

analyses = [
# Animate(;
#   fps = 20
#   ),
# NonlinearityPlot(;
#   fn_bounds = (-1,15)
#   ),
# SpaceTimePlot(),
SubsampledPlot(
  plot_type=WaveStatsPlot,
  time_subsampling=Dict(
    :Δsubsampled => 0.01,
    :scalar_window => (1.2, 1.8)
    ),
  space_subsampling=Dict(
    :scalar_window => (5.0,Inf)
    )
  )
]

output = SingleOutput(;
  root = "/home/grahams/Dropbox/simulation-73/results/",
  simulation_name = "neuman/sech2"
  )

@save "parameters.jld2" simulation

analyse.(analyses, Ref(simulation), Ref(output))