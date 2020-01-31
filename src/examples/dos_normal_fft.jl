# Aee: E->E amplitude
# See: E->E spread

# caS: cross-auto spread ratio
# caA: cross-auto amplitude ratio

# ioA: inhibitory output amplitude scale
# ioS: inhibitory output spread scale
# iiA: inhibitory input amplitude scale
# iiS: inhibitory input spread scale
ABS_STOP=300.0
dos_example = TravelingWaveSimulations.@EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; SNR_scale=80.0, stop_time=ABS_STOP,
                                                     Aee=250.0, See=25.0,
                                                     Aii=10.0, Sii=25.0,
                                                     Aie=50.0, Sie=27.0,
                                                     Aei=75.0, Sei=27.0,
                                                     n=256, x=700.0)
  simulation = Simulation(
    WCMSpatial(;
      pop_names = ("E", "I"),
      α = (1.0, 1.0),
      β = (1.0, 1.0),
      τ = (10.0, 10.0),
      nonlinearity = pops(DifferenceOfSigmoids;
        #sd = [6.7, sqrt(3.2)],
        #θ = [18.0, 10.0]),
        firing_θ = [10.0, 5.0],
        firing_a = [1.2, 1.0],
        blocking_θ = [18.0, 13.0],
        blocking_a = [1.2, 1.0]),
      stimulus = pops(SharpBumpStimulusParameter;
          strength = [10.0, 0.0],
          width = [28.1, 28.1],
          time_windows = [[(0.0, 10.0)], [(0.0, 10.0)]],
          baseline=[0.0, 0.0]),
      connectivity = FFTParameter(pops(GaussianConnectivityParameter;
          amplitude = [Aee -Aei;
                       Aie -Aii],
          spread = [(See,See) (Sei,Sei);
                    (Sie,Sie) (Sii,Sii)]
         ))
      );
      space = PeriodicLattice{Float64,N_ARR}(; n_points=(n,n), extent=(x,x)),
      save_idxs = RadialSlice(),
      tspan = (0.0,stop_time),
      dt = 1.0,
      algorithm=Tsit5(),
      callback=DiscreteCallback(if !(save_idxs === nothing)
        (u,t,integrator) -> begin
                    sub_u = u[integrator.opts.save_idxs];
                    (all(isapprox.(sub_u, 0.0, atol=1e-4)) || (sub_u[end] > 0.01)) && t > 5
                end
    else
        (u,t,integrator) -> begin
                    pop = population(u,1)
                    (all(isapprox.(u, 0.0, atol=1e-4)) || (sum(pop[:,end]) / size(pop,1) > 0.01)) && t > 5
            end
    end, terminate!)
  )
end

