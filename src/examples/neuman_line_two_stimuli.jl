@EI_kw_example function neuman_line_two_stimuli(N=1, P=2)
    T = Float64
    simulation = Simulation(;
      model = WCMSpatial{T,N,P}(;
        pop_names = ["E", "I"],
        α = [1.1, 1.0],
        β = [1.1, 1.1],
        τ = [10.0, 18.0],
        space = Segment{T}(; n_points=(301,), extent=(1000.0,)),
        nonlinearity = pops(SigmoidNonlinearity{T};
          a = [1.2, 1.0],
          θ = [2.6, 8.0]
        ),
        stimulus = MultipleDifferentStimuli{T,N}([
            pops(GaussianNoiseStimulus{T,N}; SNR = [80.0, 80.0]),
            pops(MultipleSameStimuli{SharpBumpStimulus{T,N},2};
              strength=[1.2, 1.2],
              width=[28.1, 28.1],
              time_windows = [[(0.0, 55.0)], [(0.0, 55.0)]]
              center = [((100.0,), (-100.0)), ((100.0,), (-100.0))]
            )
          ]
        )
        connectivity = pops(ExpSumAbsDecayingConnectivity{T,N};
          amplitude = [16.0 -18.2;
                       27.0 -4.0],
          spread = [70.0 90.0;
                    90.0 70.0] .|> (x) -> (x,)
        )
      ),
      solver = Solver{T}(;
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
