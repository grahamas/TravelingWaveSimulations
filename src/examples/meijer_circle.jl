@EI_kw_example function example(N=1, P=2)
    simulation = Simulation(;
      model = WCMSpatial{Float64,N,P}(;
        pop_names = ["E", "I"],
        α = [1.1, 1.0],
        β = [1.1, 1.1],
        τ = [10.0, 10.0], # In ms
        space = Circle{Float64}(; n_points=(301,), extent=(1000.0,)),
        nonlinearity = pops(GaussianNonlinearity{Float64};
          sd = [6.7, 3.2],
          θ = [18.0, 10.0]),
        stimulus = [NoisyStimulus{Float64,N}(;
          strength=10.0,
          time_windows=[(0.0,10.0)],
          width=100.0,
          SNR=80.0,
          mean=1.0,
          stim_type=SharpBumpStimulus{Float64,N}),
          NoStimulus{Float64,N}()],
        connectivity = pops(ExpSumSqDecayingConnectivity{Float64,N};
          amplitude = [280.0 -297.0;
                       270.0 -1.4],
          spread = [(70.0,) (90.0,);
                    (90.0,) (70.0,)]
          )
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
