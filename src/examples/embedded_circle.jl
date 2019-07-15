@EI_kw_example function example(N_ARR=2, N_CDT=3, P=2)
    simulation = Simulation(;
      model = WCMSpatial{Float64,N_ARR,N_CDT,P}(;
        pop_names = ["E", "I"],
        α = [1.1, 1.0],
        β = [1.1, 1.1],
        τ = [10.0, 10.0], # In ms
        space = RandomlyEmbeddedLattice(;
          lattice=Grid{Float64}(; n_points=(51,51), extent=(500.0,500.0)),
          embedded_lattice=Circle{Float64}(; n_points=(10,), extent=(2π,))
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
        connectivity = pops(ExpSumSqDecayingConnectivity{Float64,N_CDT};
          amplitude = [280.0 -297.0;
                       270.0 -1.4],
          spread = [(70.0, 70.0, 2π/10.0) (90.0, 90.0, 2π/10.0);
                    (90.0, 90.0, 2π/10.0) (70.0, 70.0, 2π/10.0)])
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
