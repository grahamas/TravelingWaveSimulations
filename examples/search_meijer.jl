# ENV["GKSwstype"] = "100" # For headless plotting (on server)
# ENV["MPLBACKEND"]="Agg"

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
  varying_model = WCMSpatial{MV,N,P}(;
    pop_names = ["E", "I"],
    α = MV[1.0, 1.0],
    β = MV[1.0, 1.0],
    τ = MV[10.0, 10.0],#MV[BV(10.0, (9.5, 10.1)), BV(10.0, (9.5, 10.1))],
    space = PopSegment{MV,P}(; n_points=301, extent=1000.0),
    nonlinearity = pops(GaussianNonlinearity{MV};
      sd = MV[6.7, 3.2],
      θ = MV[BV(10.0, (5.0, 20.0)), BV(15.0, (5.0, 20.0))]),
    # nonlinearity = pops(SigmoidNonlinearity{MV};
    #   a = MV[2.0, 0.95],
    #   θ = [BV(12.41, (5.0, 15.0)), BV(7.33, (5.0, 15.0))]
    # ),
      stimulus = [NoisyStimulus{MV}(;
        strength=BV(10.0, (1.0, 11.0)),
        window=Tuple{MV,MV}((0.0,10.0)),
        width=100.0,
        mean = 1.0,
        SNR = 80.0,
        stim_type=SharpBumpStimulus{MV}),
        NoStimulus{MV}()
      ],
      connectivity = pops(MeijerConnectivity{MV};
        amplitude = MV[BV(1.0, (0.0, 3.0)) BV(-1.0, (-3.0, 0.0));
                      BV(1.0, (0.0, 3.0)) BV(-1.0, (-3.0, 0.0))],
        spread = MV[70.0 90.0;
                   90.0 70.0])
      ),
    solver = Solver{v}(;
      stop_time = 100.0,
      dt = 1.0,
      space_save_every=1,
      time_save_every=1,
      algorithm=Euler()
      ),
  target=MatchData(
    file_name="data/meijer_wave_example.jld2"
    ),
  MaxSteps=10000
  )

analyses = [
Animate(;
  fps = 15
  ),
NonlinearityPlot(;
  fn_bounds = (-1,20)
  ),
# SpaceTimePlot(),
SubsampledPlot(
  plot_type=WaveStatsPlot,
  time_subsampler=Subsampler(
    Δ = 1.0,
    window = (30.0, 80.0)
  ),
  space_subsampler=Subsampler(
      window = (50.0,Inf)
    )
  )
]

output = SingleOutput(;
  root = joinpath(homedir(), "simulation-results"),
  simulation_name = "search/meijer/sigmoid"
  )

output(((name, obj) -> @save name obj), "parameters.jld2", p_search)

analyse.(analyses, Ref(p_search.result_simulation), Ref(output))
