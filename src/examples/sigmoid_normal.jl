@EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; spread_scale=1.0, amplitude_scale=1.0, 
                                                     SNR_scale=80.0, auto=25.0, 
                                                     ca_ratio=(27.0/25.0), cross=(ca_ratio*auto),
                                                     inh_scale=1.0)
  simulation = Simulation(;
    model = WCMSpatial{Float64,N_ARR,N_CDT,P}(;
      pop_names = ["E", "I"],
      α = [1.0, 1.0],
      β = [1.0, 1.0],
      τ = [10.0, 10.0],
      space = CompactLattice{Float64,N_ARR}(; n_points=(35,35), extent=(350.0,350.0)),
      nonlinearity = pops(SigmoidNonlinearity{Float64};
        a = [1.2, 1.0],
        θ = [2.6, 8.0]),
      stimulus = pops(NoisyStimulus{Float64,N_CDT};
          strength = [1.2, 1.2],
          width = [28.1, 28.1],
          SNR = [1.0, 1.0] .* SNR_scale,
          time_windows = [[(0.0, 45.0)], [(0.0, 45.0)]],
          stim_type=[SharpBumpStimulus{Float64,N_CDT}, SharpBumpStimulus{Float64,N_CDT}]),
      connectivity = pops(ExpSumSqDecayingConnectivity{Float64,N_CDT};
          amplitude = [16.0 (-18.2 * inh_scale);
                       27.0 (-4.0 * inh_scale)] .* amplitude_scale,
          #spread = [(70.0,70.0) (90.0,90.0);
          #          (90.0,90.0) (70.0,70.0)] .|> (tup) -> map((x) -> x * spread_scale, tup)
          spread = [(auto,auto) (cross,cross);
                    (cross,cross) (auto,auto)] .|> (tup) -> map((x) -> x * spread_scale, tup)
          )
      ),
    solver = Solver{Float64}(;
        stop_time = 100.0,
        dt = 1.0,
        space_save_every=1,
        time_save_every=1,
        algorithm=Euler()
      )
  )
end
