@EI_kw_example function neuman_square(N=2,P=2)
  simulation = Simulation(;
    model = WCMSpatial{Float64,N,P}(;
      pop_names = ["E", "I"],
      α = [1.1, 1.0],
      β = [1.1, 1.1],
      τ = [0.1, 0.18],
      space = Lattice{Float64,N}(; n_points=(51,51), extent=(50.0,50.0)),
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
      connectivity = pops(ExpSumSqDecayingConnectivity{Float64,N};
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
