@EI_kw_example function example(N_ARR=2,N_CDT=2, P=2; sep=70.0)
    simulation = Simulation(;
      model = WCMSpatial{Float64,N_ARR,N_CDT,P}(;
        pop_names = ["E", "I"],
        α = [1.1, 1.0],
        β = [1.1, 1.1],
        τ = [10.0, 10.0], # In ms
        space = CompactLattice{Float64,N_ARR}(; n_points=(51,51), extent=(500.0,500.0)),
        nonlinearity = pops(GaussianNonlinearity{Float64};
          sd = [6.7, 3.2],
          θ = [18.0, 10.0]),
        stimulus = [MultipleDifferentStimuli{Float64,2}([
            GaussianNoiseStimulus{Float64,2}(
                SNR=80.0,
                mean=1.0
            ),
            MultipleSameStimuli{Float64,N_CDT,SharpBumpStimulus{Float64,N_CDT},2}(;
                strength=10.0,
                time_windows=[(0.0,10.0)],
                width=100.0,
                center = ((sep,sep), (-sep,-sep))
            )
            ]),
          NoStimulus{Float64,N_CDT}()],
        connectivity = pops(ExpSumSqDecayingConnectivity{Float64,N_CDT};
          amplitude = [280.0 -297.0;
                       270.0 -1.4],
          spread = [(70.0, 70.0) (90.0, 90.0);
                    (90.0, 90.0) (70.0, 70.0)])
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
