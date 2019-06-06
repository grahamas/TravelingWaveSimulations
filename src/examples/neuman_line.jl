@EI_kw_example function neuman_line(N=1, P=2)
    simulation = Simulation(;
      model = WCMSpatial{Float64,N,P}(;
        pop_names = ["E", "I"],
        α = [1.1, 1.0],
        β = [1.1, 1.1],
        τ = [10.0, 18.0],
        space = Segment{Float64}(; n_points=(301,), extent=(1000.0,)),
        nonlinearity = pops(SigmoidNonlinearity{Float64};
          a = [1.2, 1.0],
          θ = [2.6, 8.0]
        ),
        stimulus = pops(NoisyStimulus{Float64,N};
            strength = [1.2, 1.2],
            width = [28.1, 28.1],
            SNR = [80.0, 80.0],
            time_window = [(0.0, 55.0), (0.0, 55.0)],
            stim_type=[SharpBumpStimulus{Float64,N}, SharpBumpStimulus{Float64,N}]
        ),
        connectivity = pops(ExpSumAbsDecayingConnectivity{Float64,N};
          amplitude = [16.0 -18.2;
                       27.0 -4.0],
          spread = [70.0 90.0;
                    90.0 70.0] .|> (x) -> (x,))
        ),
      solver = Solver{Float64}(;
        stop_time = 180.0,
        dt = 1.0,
        space_save_every=1,
        time_save_every=1,
        #stiffness=:stiff
        algorithm=Euler()
      )
    )
    return simulation
end
