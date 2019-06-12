@EI_kw_example function example(N=2,P=2; spread_scale=1.0, amplitude_scale=1.0)
  simulation = Simulation(;
    model = WCMSpatial{Float64,N,P}(;
      pop_names = ["E", "I"],
      α = [1.1, 1.0],
      β = [1.1, 1.1],
      τ = [10.0, 18.0],
      space = Lattice{Float64,N}(; n_points=(51,51), extent=(500.0,500.0)),
      nonlinearity = pops(SigmoidNonlinearity{Float64};
        a = [1.2, 1.0],
        θ = [2.6, 8.0]),
      stimulus = pops(NoisyStimulus{Float64,N};
          strength = [1.2, 1.2],
          width = [28.1, 28.1],
          SNR = [80.0, 80.0],
          time_windows = [[(0.0, 55.0)], [(0.0, 55.0)]],
          stim_type=[SharpBumpStimulus{Float64,N}, SharpBumpStimulus{Float64,N}]),
      connectivity = pops(ExpSumSqDecayingConnectivity{Float64,N};
          amplitude = [16.0 -18.2;
                       27.0 -4.0] .* amplitude_scale,
          #spread = [(70.0,70.0) (90.0,90.0);
          #          (90.0,90.0) (70.0,70.0)] .|> (tup) -> map((x) -> x * spread_scale, tup)
          spread = [(25.0,25.0) (27.0,27.0);
                    (27.0,27.0) (25.0,25.0)] .|> (tup) -> map((x) -> x * spread_scale, tup)
          )
      ),
    solver = Solver{Float64}(;
        stop_time = 180.0,
        dt = 1.0,
        space_save_every=1,
        time_save_every=1,
        algorithm=Euler()
      )
  )
end
