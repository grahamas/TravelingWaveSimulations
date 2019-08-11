@EI_kw_example function example(N_ARR=2, N_CDT=3, P=2; circle_spread=0.4,
      circle_auto_cross_ratio = 7.0/9.0,
      auto_spread=70.0, cross_spread=90.0,
      circle_auto_spread=(circle_spread * circle_auto_cross_ratio),
      circle_cross_spread=circle_spread,
      grid_axis_n_points=51,
      grid_axis_extent=500.0,
      circle_n_points=10,
      circle_extent=2π)
    simulation = Simulation(;
      model = WCMSpatial{Float64,N_ARR,N_CDT,P}(;
        pop_names = ["E", "I"],
        α = [1.1, 1.0],
        β = [1.1, 1.1],
        τ = [10.0, 10.0], # in ms
        space = RandomlyEmbeddedLattice(;
          lattice=CompactLattice{Float64,2}(; n_points=round.(Ref(Int),(grid_axis_n_points,grid_axis_n_points)),
            extent=(grid_axis_extent,grid_axis_extent)),
          embedded_lattice=PeriodicLattice(; n_points=(round(Int, circle_n_points),),
            extent=(circle_extent,))
        ),
        nonlinearity = pops(GaussianNonlinearity{Float64};
          sd = [6.7, sqrt(3.2)],
          θ = [18.0, 10.0]),
        stimulus = [NoisyStimulus{Float64,N_CDT}(;
          strength=10.0,
          time_windows=[(0.0,10.0)],
          width=100.0,
          SNR=80.0,
          mean=1.0,
          stim_type=SharpBumpStimulus{Float64,N_CDT}),
          NoStimulus{Float64,N_CDT}()],
        connectivity = pops(ExpSumSqDecayingAnisotropicConnectivity{Float64,N_CDT};
          amplitude = [280.0 -297.0;
                       270.0 -1.4],
          spread = [(auto_spread, auto_spread, circle_auto_spread) (cross_spread, cross_spread, circle_cross_spread);
                    (cross_spread, cross_spread, circle_cross_spread) (auto_spread, auto_spread, circle_auto_spread)])
        ),
      solver = Solver{Float64}(;
        stop_time = 100.0,
        dt = 1.0,
        space_save_every=1,
        time_save_every=1,
        #stiffness=:stiff
        algorithm=Euler()
      )
    )
    return simulation
end
