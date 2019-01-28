ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"

using Simulation73
using WilsonCowanModel
using DifferentialEquations: Euler
using JLD2
using Plots; pyplot()

if !@isdefined(v)
  const v = MaybeVariable{Float64}
  const BV = BoundedVariable
end
N=1
P=2
p_search = ParameterSearch(;
  varying_model = WCMSpatial1D{v,N,P}(;
    pop_names = ["E", "I"],
    α = v[1.1, 1.0],
    β = v[1.1, 1.1],
    τ = v[BV(0.1, (0.05, 0.2)), BV(0.18, (0.05,0.2))],#v[0.1, 0.18], #REVERSED
    space = PopSegment{v,P}(; n_points=301, extent=100.0),
    nonlinearity = pops(Sech2Nonlinearity{v}; 
      a = v[1.2, 1.0],
      θ = v[BV(2.6, (2.0,5.0)), BV(8.0, (5.0,9.0))]#v[BV(2.6,(1.0,5.0)), BV(8.0,(5.0,9.0))]
      # a=v[BV(1.2, (0.1, 2.0)), BV(1.0, (0.1, 2.0))],
      # θ=v[BV(8.0, (2.0, 9.0)), BV(2.6, (2.0,9.0))]
      ),
    stimulus = pops(NoisySharpBumpStimulus{v}; 
      strength=v[1.2, 1.2],
      window=Tuple{v,v}[(0.0,0.55), (0.0,0.55)],
      width=v[2.81, 2.81],
      SNR=v[80.0, 80.0]),
    connectivity = pops(ShollConnectivity{v};
      amplitude = v[16.0 -18.2;
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
    ),
  target=MatchExample(
    file_name="data/wave_example.jld2"
    ),
  MaxSteps=10000
  )

analyses = [
Animate(;
  fps = 20
  ),
# NonlinearityPlot(;
#   fn_bounds = (-1,15)
#   ),
# SpaceTimePlot(),
SubsampledPlot(;
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
  simulation_name = "search/sech2"
  )

@save "parameters.jld2" p_search

analyse.(analyses, Ref(result_simulation(p_search)), Ref(output))