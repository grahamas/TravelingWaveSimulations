# Aee: E->E amplitude
# See: E->E spread

# caS: cross-auto spread ratio
# caA: cross-auto amplitude ratio

# ioA: inhibitory output amplitude scale
# ioS: inhibitory output spread scale
# iiA: inhibitory input amplitude scale
# iiS: inhibitory input spread scale

@EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; SNR_scale=80.0, stop_time=30.0,
                                                     Aee=24.0, See=25.0,
                                                     Aii=4.0, Sii=27.0,
                                                     Aie=27.0, Sie=25.0,
                                                     Aei=18.2, Sei=27.0,
                                                     n=128, x=700.0, stim_strength=1.2)
  simulation = Simulation(
    WCMSpatial(;
      pop_names = ("E", "I"),
      α = (1.0, 1.0),
      β = (1.0, 1.0),
      τ = (3.0, 3.0),
      nonlinearity = pops(SigmoidNonlinearity;
        a = [1.2, 1.0],
        θ = [2.6, 8.0]),
      stimulus = [pops(GaussianNoiseStimulusParameter; SNR=[1.0, 1.0] .* SNR_scale), pops(SharpBumpStimulusParameter;
          strength = [stim_strength,stim_strength],
          width = [28.1, 28.1],
          time_windows = [[(0.0, 5.0)], [(0.0, 5.0)]])],
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
      algorithm=Euler()
  )
end
