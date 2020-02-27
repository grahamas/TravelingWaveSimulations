ABS_STOP=300.0
TravelingWaveSimulations.@EI_kw_example function example(N_ARR=1,N_CDT=1,P=2; SNR_scale=80.0, stop_time=ABS_STOP,
                                                     Aee=70.0, See=25.0,
                                                     Aii=2.0, Sii=27.0,
                                                     Aie=35.0, Sie=25.0,
                                                     Aei=70.0, Sei=27.0,
                                                     n=256, x=700.0, 
                                                     stim_strength=6.0,
                                                     stim_width=28.1,
                                                     stim_duration=7.0,
                                                     save_idxs_arg=[IndexSubsampler((5,)), RightCutFromValue((0.0,))])
  simulation = Simulation(
    WCMSpatial(;
      pop_names = ("E", "I"),
      α = (1.0, 1.0),
      β = (1.0, 1.0),
      τ = (3.0, 3.0),
      nonlinearity = pops(DifferenceOfSigmoids;
        firing_θ = [6.0, 11.4],
        firing_a = [1.2, 1.0],
        blocking_θ = [30.0, 30.0],
        blocking_a = [1.2, 1.0]),
      stimulus = pops(SharpBumpStimulusParameter;
          strength = [stim_strength, stim_strength],
          width = [stim_width, stim_width],
          time_windows = [[(0.0, stim_duration)], [(0.0, stim_duration)]],
          baseline=[0.0, 0.0]),
      connectivity = FFTParameter(pops(GaussianConnectivityParameter;
          amplitude = [Aee -Aei;
                       Aie -Aii],
          spread = [(See,) (Sei,);
                    (Sie,) (Sii,)]
         ))
      );
      space = PeriodicLattice{Float64,N_ARR}(; n_points=(n,), extent=(x,)),
      tspan = (0.0,stop_time),
      dt = 0.1,
      algorithm=Tsit5(),
      save_idxs = save_idxs_arg,
      callback=DiscreteCallback(if !(save_idxs_arg === nothing)
        (u,t,integrator) -> begin
                    sub_u = u[integrator.opts.save_idxs];
                    t > 5 && ((all(isapprox.(sub_u, 0.0, atol=1e-4)) || (sub_u[end] > 0.005)))
                end
    else
        (u,t,integrator) -> begin
                    pop = population(u,1)
                    t > 5 && ((all(isapprox.(u, 0.0, atol=1e-4)) || (pop[end] > 0.005)))
            end
    end, terminate!)
  )
end

