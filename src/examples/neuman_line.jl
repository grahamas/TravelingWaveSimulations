@EI_kw_example function neuman_line(N=1, P=2)
    simulation = Simulation(;
      model = WCMSpatial{Float64,N,P}(;
        pop_names = ["E", "I"],
        α = [1.1, 1.0],
        β = [1.1, 1.1],
        τ = [0.1, 0.18],
        space = Pops{P}(Segment{Float64}(; n_points=301, extent=100.0)),
        nonlinearity = pops(SigmoidNonlinearity{Float64};
          a = [1.2, 1.0],
          θ = [2.6, 8.0]),
        stimulus = pops(NoisyStimulus{Float64,N};
          strength=[1.2, 1.2],
          time_window=[(0.0, 0.55), (0.0, 0.55)],
          width=[2.81, 2.81],
          SNR=[80.0, 80.0],
          stim_type=[SharpBumpStimulus{Float64,N}, SharpBumpStimulus{Float64,N}]),
        connectivity = pops(ShollConnectivity{Float64};
          amplitude = [16.0 -18.2;
                       27.0 -4.0],
          spread = [2.5 2.7;
                    2.7 2.5])
        ),
      solver = Solver{Float64}(;
        stop_time = 1.8,
        dt = 0.01,
        space_save_every=1,
        time_save_every=1,
        #stiffness=:stiff
        algorithm=Euler()
      )
    )
    return simulation
end
