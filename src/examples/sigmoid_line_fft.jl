ABS_STOP = 300.0
@EI_kw_example function example(N_ARR=1,N_CDT=1,P=2; stop_time=ABS_STOP,
                                                     Aee=24.0, See=25.0,
                                                     Aii=4.0, Sii=27.0,
                                                     Aie=27.0, Sie=25.0,
                                                     Aei=18.2, Sei=27.0,
                                                     n=128, x=700.0, stim_strength=1.2,
                                                     stim_width=28.1)
  simulation = Simulation(
    WCMSpatial(;
      pop_names = ("E", "I"),
      α = (1.0, 1.0),
      β = (1.0, 1.0),
      τ = (3.0, 3.0),
      nonlinearity = pops(SigmoidNonlinearity;
        a = [1.2, 1.0],
        θ = [2.6, 8.0]),
      stimulus =  pops(SharpBumpStimulusParameter;
          strength = [stim_strength,stim_strength],
          width = [stim_width, stim_width],
          time_windows = [[(0.0, 5.0)], [(0.0, 5.0)]]),
      connectivity = FFTParameter(pops(GaussianConnectivityParameter;
          amplitude = [Aee -Aei;
                       Aie -Aii],
          spread = [(See,) (Sei,);
                    (Sie,) (Sii,)]
          ))
      );
      space = PeriodicLattice{Float64,N_ARR}(; n_points=(n,), extent=(x,)),
      save_idxs=nothing,
      tspan = (0.0,stop_time),
      dt = 1.0,
      algorithm=Euler(),
      callback=DiscreteCallback(if !(save_idxs === nothing)
        (u,t,integrator) -> begin
                    sub_u = u[integrator.opts.save_idxs];
                    (all(sub_u .≈ 0.0) || (sub_u[end] > 0.01)) && t > 5
                end
    else
        (u,t,integrator) -> begin
                    pop = population(u,1)
                    (all(u .≈ 0.0) || (sum(pop[:,end]) / size(pop,1) > 0.01)) && t > 5
            end
    end, terminate!)
  )
end
