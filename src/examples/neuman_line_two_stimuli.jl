@EI_kw_example function example(N=1, P=2; SNR_scale = 80.0, sep=100.0)
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
        stimulus = MultipleDifferentStimuli{Float64,N}([
            pops(GaussianNoiseStimulus{Float64,N}; SNR = [1.0, 1.0] .* SNR_scale),
            pops(MultipleSameStimuli{Float64,N,SharpBumpStimulus{Float64,N},2};
              strength=[1.2, 1.2],
              width=[28.1, 28.1],
              time_windows = [[(0.0, 55.0)], [(0.0, 55.0)]],
              center = [((sep,), (-sep,)), ((sep,), (-sep,))]
            )
          ]
        ),
        connectivity = pops(ExpSumAbsDecayingConnectivity{Float64,N};
          amplitude = [16.0 -18.2;
                       27.0 -4.0],
          spread = [70.0 90.0;
                    90.0 70.0] .|> (x) -> (x,)
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
