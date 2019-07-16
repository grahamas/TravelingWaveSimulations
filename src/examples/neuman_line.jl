@EI_kw_example function example(N_ARR=1,N_CDT=1, P=2; SNR_scale=80.0)
    simulation = Simulation(;
      model = WCMSpatial{Float64,N_ARR,N_CDT,P}(;
        pop_names = ["E", "I"],
        α = [1.1, 1.0],
        β = [1.1, 1.1],
        τ = [10.0, 18.0],
        space = Segment{Float64}(; n_points=(301,), extent=(1000.0,)),
        nonlinearity = pops(SigmoidNonlinearity{Float64};
          a = [1.2, 1.0],
          θ = [2.6, 8.0]
        ),
        stimulus = pops(NoisyStimulus{Float64,N_CDT};
            strength = [1.2, 1.2],
            width = [28.1, 28.1],
            SNR = [1.0, 1.0] .* SNR_scale,
            time_windows = [[(0.0, 55.0)], [(0.0, 55.0)]],
            stim_type=[SharpBumpStimulus{Float64,N_CDT}, SharpBumpStimulus{Float64,N_CDT}]
        ),
        connectivity = pops(ExpSumAbsDecayingConnectivity{Float64,N_CDT};
          amplitude = [16.0 -18.2;
                       27.0 -4.0],
          # spread = [70.0 90.0;
          #           90.0 70.0] .|> (x) -> (x,))
          spread = [(70.0,) (90.0,);
                    (90.0,) (70.0,)]
        )
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
