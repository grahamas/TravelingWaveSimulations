# Aee: E->E amplitude
# See: E->E spread

# caS: cross-auto spread ratio
# caA: cross-auto amplitude ratio

# ioA: inhibitory output amplitude scale
# ioS: inhibitory output spread scale
# iiA: inhibitory input amplitude scale
# iiS: inhibitory input spread scale

@EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; SNR_scale=80.0, 
                                                     Aee=16.0, See=25.0,
                                                     iiA=sqrt((27.0/4.0)/18.2), iiS=1.0,
                                                     ioA=1.0/(4iiA), ioS=1.0,
                                                     caA=27.0/(16iiA), caS=(27.0/25.0),
                                                     Aii=Aee*iiA*ioA, Sii=See*iiS*ioS,
                                                     Aie=Aee*iiA*caA, Sie=See*iiS*caS,
                                                     Aei=Aee*ioA*caA, Sei=See*ioS*caS)
  simulation = Simulation(;
    model = WCMSpatial{Float64,N_ARR,N_CDT,P}(;
      pop_names = ["E", "I"],
      α = [1.0, 1.0],
      β = [1.0, 1.0],
      τ = [10.0, 10.0],
      space = CompactLattice{Float64,N_ARR}(; n_points=(35,35), extent=(350.0,350.0)),
      nonlinearity = pops(SigmoidNonlinearity{Float64};
        a = [1.2, 1.0],
        θ = [2.6, 8.0]),
      stimulus = pops(NoisyStimulus{Float64,N_CDT};
          strength = [1.2, 1.2],
          width = [28.1, 28.1],
          SNR = [1.0, 1.0] .* SNR_scale,
          time_windows = [[(0.0, 45.0)], [(0.0, 45.0)]],
          stim_type=[SharpBumpStimulus{Float64,N_CDT}, SharpBumpStimulus{Float64,N_CDT}]),
      connectivity = pops(ExpSumSqDecayingConnectivity{Float64,N_CDT};
          amplitude = [Aee -Aei;
                       Aie -Aii],
          spread = [(See,See) (Sei,Sei);
                    (Sie,Sie) (Sii,Sii)]
          )
      ),
    solver = Solver{Float64}(;
        stop_time = 100.0,
        dt = 1.0,
        space_save_every=1,
        time_save_every=1,
        algorithm=Euler()
      )
  )
end
