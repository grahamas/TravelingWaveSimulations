# Aee: E->E amplitude
# See: E->E spread

# caS: cross-auto spread ratio
# caA: cross-auto amplitude ratio

# ioA: inhibitory output amplitude scale
# ioS: inhibitory output spread scale
# iiA: inhibitory input amplitude scale
# iiS: inhibitory input spread scale

@EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; SNR_scale=80.0, stop_time=100.0, 
                                                     Aee=280.0, See=70.0,
                                                     iiA=sqrt((270.0/1.4)/297.0), iiS=1.0,
                                                     ioA=1.0/(1.4iiA), ioS=1.0,
                                                     caA=270.0/(280.0iiA), caS=(90.0/70.0),
                                                     Aii=Aee*iiA*ioA, Sii=See*iiS*ioS,
                                                     Aie=Aee*iiA*caA, Sie=See*iiS*caS,
                                                     Aei=Aee*ioA*caA, Sei=See*ioS*caS,
                                                     n=35, x=350.0)
  simulation = Simulation(
    WCMSpatial{Float64,N_CDT,P}(;
      pop_names = ["E", "I"],
      α = [1.0, 1.0],
      β = [1.1, 1.1],
      τ = [10.0, 10.0],
      nonlinearity = pops(GaussianNonlinearity{Float64};
        sd = [6.7, sqrt(3.2)],
        θ = [18.0, 10.0]),
      stimulus = [NoisyStimulus{Float64,N_CDT}(;
          strength = 10.0,
          width = 100.0,
          SNR = SNR_scale,
          mean=1.0,
          time_windows = [(0.0, 10.0)],
          stim_type=SharpBumpStimulus{Float64,N_CDT}),
          NoStimulus{Float64,N_CDT}()],
      connectivity = pops(ExpSumSqDecayingConnectivity{Float64,N_CDT};
          amplitude = [Aee -Aei;
                       Aie -Aii],
          spread = [(See,See) (Sei,Sei);
                    (Sie,Sie) (Sii,Sii)]
          )
      );
      space = CompactLattice{Float64,N_ARR}(; n_points=(n,n), extent=(x,x)),
      tspan = (0.0,stop_time),
      dt = 1.0,
      algorithm=Euler()
  )
end
