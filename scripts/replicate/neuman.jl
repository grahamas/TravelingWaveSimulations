# ENV["GKSwstype"] = "100" # For headless plotting (on server)
# ENV["MPLBACKEND"]="Agg"
using DrWatson
quickactivate(@__DIR__, "WilsonCowanModel")
using Simulation73
using WilsonCowanModel
using DifferentialEquations: Euler
using Random
using Dates

Random.seed!(9999)

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
    space = Pops{P}(Segment{v}(; n_points=301, extent=100.0)),
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
    connectivity = pops(ShollConnectivity{v};
      amplitude = v[16.0 -18.2
                    27.0 -4.0],
      # spread = v[BV(2.5, (1.0, 4.0)) BV(2.7, (1.0, 4.0));
      #            BV(2.7, (1.0, 4.0)) BV(2.5, (1.0, 4.0))])
      spread = v[2.5 2.7;
                 2.7 2.5])
    ),
  solver = Solver{v}(;
    stop_time = 1.8,
    dt = 0.01,
    space_save_every=1,
    time_save_every=1,
    #stiffness=:stiff
    algorithm=Euler()
  )
)

execution = execute(simulation);

replication_directory = joinpath(datadir(), "sim", "replicate", "neuman")
this_time_commit_filename = join([Dates.now(),current_commit()*".bson"], "_")
mkpath(replication_directory)
tagsave(joinpath(replication_directory, this_time_commit_filename), @dict execution; safe=true)



# analyses = [
# Animate(;
#   fps = 20
#   ),
# NonlinearityPlot(;
#   fn_bounds = (-1,15)
#   ),
# # SpaceTimePlot(),
# SubsampledPlot(
#   plot_type=WaveStatsPlot,
#   time_subsampler=Subsampler(
#     Δ = 0.01,
#     window = (1.2, 1.8)
#   ),
#   space_subsampler=Subsampler(
#       window = (5.0,Inf)
#     )
#   )
# ]

# analyse.(analyses, Ref(simulation), Ref(output))
