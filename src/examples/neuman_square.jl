@EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; spread_scale=1.0, amplitude_scale=1.0, SNR_scale=80.0, stop_time=180.0,
                                                     auto=25.0, ca_ratio=(27.0/25.0), cross=(ca_ratio*auto), n=51, x=500.0)
  simulation = Simulation(
   WCMSpatial{Float64,N_CDT,P}(;
      pop_names = ["E", "I"],
      α = [1.1, 1.0],
      β = [1.1, 1.1],
      τ = [10.0, 18.0],
      nonlinearity = pops(SigmoidNonlinearity{Float64};
        a = [1.2, 1.0],
        θ = [2.6, 8.0]),
      stimulus = pops(NoisyStimulus{Float64,N_CDT};
          strength = [1.2, 1.2],
          width = [28.1, 28.1],
          SNR = [1.0, 1.0] .* SNR_scale,
          time_windows = [[(0.0, 55.0)], [(0.0, 55.0)]],
          stim_type=[SharpBumpStimulus{Float64,N_CDT}, SharpBumpStimulus{Float64,N_CDT}]),
      connectivity = pops(ExpSumSqDecayingConnectivity{Float64,N_CDT};
          amplitude = [16.0 -18.2;
                       27.0 -4.0] .* amplitude_scale,
          spread = [(auto,auto) (cross,cross);
                    (cross,cross) (auto,auto)] .|> (tup) -> map((x) -> x * spread_scale, tup)
          )
      );
    space = CompactLattice{Float64,N_ARR}(; n_points=(n,n), extent=(x,x)),
    tspan=(0.0, stop_time),
    dt = 1.0,
    algorithm=Euler()
  )
end
