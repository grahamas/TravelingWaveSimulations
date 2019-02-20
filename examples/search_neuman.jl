ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"

using Simulation73
using WilsonCowanModel
using DifferentialEquations: Euler
using JLD2
using Plots; pyplot()

if !@isdefined(v)
  const v = Float64
  const MV = MaybeVariable{v}
  const BV = BoundedVariable{v}
end
N=1
P=2
p_search = ParameterSearch(;
  varying_model = WCMSpatial1D{MV,N,P}(;
    pop_names = ["E", "I"],
    α = MV[1.1, 1.0],
    β = MV[1.1, 1.1],
    τ = MV[BV(0.15, (0.05, 0.2)), BV(0.09, (0.05,0.2))],#MV[0.1, 0.18], #REVERSED
    space = PopSegment{MV,P}(; n_points=301, extent=100.0),
    nonlinearity = pops(SigmoidNonlinearity{MV};
      a = MV[1.2, 1.0],
      θ = MV[2.6, 8.0]#MV[BV(2.6,(1.0,5.0)), BV(8.0,(5.0,9.0))]
      # a=MV[BV(1.2, (0.1, 2.0)), BV(1.0, (0.1, 2.0))],
      # θ=MV[BV(8.0, (2.0, 9.0)), BV(2.6, (2.0,9.0))]
      ),
      stimulus = [NoisyStimulus{MV}(;
        strength=1.2,
        window=Tuple{MV,MV}((0.0,0.55)),
        width=2.81,
        SNR= 80.0,
        stim_type=SharpBumpStimulus{MV}),
        NoisyStimulus{MV}(;
          strength=1.2,
          window=Tuple{MV,MV}((0.0,0.55)),
          width=2.81,
          SNR= 80.0,
          stim_type=SharpBumpStimulus{MV})
      ],
    connectivity = pops(ShollConnectivity{MV};
      amplitude = MV[BV(10.0, (5.0, 20.0)) BV(-11.0, (-30.0, -10.0));
                    BV(15.0, (10.0, 30.0)) BV(-1.0,(-10.0, 0.0))],
      # spread = MV[BV(2.5, (1.0, 4.0)) BV(2.7, (1.0, 4.0));
      #            BV(2.7, (1.0, 4.0)) BV(2.5, (1.0, 4.0))])
      spread = MV[2.5 2.7;
                 2.7 2.5])
    ),
  solver = Solver{v}(;
    stop_time = 1.8,
    dt = 0.005,
    space_save_every=1,
    time_save_every=1,
    algorithm=Euler()
    ),
  target=MatchData(
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
SubsampledPlot(
  plot_type=WaveStatsPlot,
  time_subsampler=Subsampler(
    Δ = 0.01,
    window = (1.2, 1.8)
  ),
  space_subsampler=Subsampler(
      window = (5.0,Inf)
    )
  )
]

output = SingleOutput(;
  root = joinpath(homedir(), "simulation-results"),
  simulation_name = "search/neuman"
  )

@save "parameters.jld2" p_search

analyse.(analyses, Ref(result_simulation(p_search)), Ref(output))

output(((name, obj) -> @save name obj), "parameters.jld2", p_search)
