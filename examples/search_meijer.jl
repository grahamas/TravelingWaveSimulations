# ENV["GKSwstype"] = "100" # For headless plotting (on server)
# ENV["MPLBACKEND"]="Agg"

using Simulation73
using WilsonCowanModel
using DifferentialEquations: Euler
using JLD2
using Plots; pyplot()

if !@isdefined(MV)
  const MV = MaybeVariable{Float64}
  const BV = BoundedVariable
end
N=1
P=2
p_search = ParameterSearch(;
  varying_model = WCMSpatial1D{MV,N,P}(;
    pop_names = ["E", "I"],
    α = MV[1.0, 1.0],
    β = MV[1.0, 1.0],
    τ = MV[BV(10.0, (9.5, 10.1)), BV(10.0, (9.5, 10.1))],#v[0.1, 0.18], #REVERSED
    space = PopSegment{MV,P}(; n_points=301, extent=1000.0),
    nonlinearity = pops(GaussianNonlinearity{MV};
      sd = MV[6.7, 3.2],
      θ = MV[18.0, 10.0]),
      stimulus = [NoisyStimulus{MV}(;
        strength=10.0,
        window=Tuple{MV,MV}((0.0,10.0)),
        width=100.0,
        mean = 1.0,
        SNR = 80.0,
        stim_type=SharpBumpStimulus{MV}),
        NoStimulus{MV}()
      ],
      connectivity = pops(MeijerConnectivity{MV};
        amplitude = MV[2.0 -1.65;
                      1.5 -0.01],
        spread = MV[70.0 90.0;
                   90.0 70.0])
      ),
    solver = Solver(;
      stop_time = 100.0,
      dt = 1.0,
      space_save_every=1,
      time_save_every=1,
      algorithm=Euler()
      ),
  target=MatchExample(
    file_name="data/meijer_wave_example.jld2"
    ),
  MaxSteps=1000
  )

analyses = [
Animate(;
  fps = 15
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
  root = joinpath(homedir(), "simulation-results"),
  simulation_name = "search/meijer"
  )

output(((name, obj) -> @save name obj), "parameters.jld2", p_search)

analyse.(analyses, Ref(p_search.result_simulation), Ref(output))
