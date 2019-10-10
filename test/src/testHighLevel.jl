TravelingWaveSimulations.@EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; SNR_scale=80.0, stop_time=100.0,
                                                     Aee=16.0, See=25.0,
                                                     iiA=sqrt((27.0/4.0)/18.2), iiS=1.0,
                                                     ioA=1.0/(4iiA), ioS=1.0,
                                                     caA=27.0/(16iiA), caS=(27.0/25.0),
                                                     Aii=Aee*iiA*ioA, Sii=See*iiS*ioS,
                                                     Aie=Aee*iiA*caA, Sie=See*iiS*caS,
                                                     Aei=Aee*ioA*caA, Sei=See*ioS*caS,
                                                     n=36, x=350.0, stim_strength=1.2)
  simulation = Simulation(
    WCMSpatial(;
      pop_names = ("E", "I"),
      α = (1.0, 1.0),
      β = (1.0, 1.0),
      τ = (10.0, 10.0),
      nonlinearity = pops(SigmoidNonlinearity;
        a = [1.2, 1.0],
        θ = [2.6, 8.0]),
      stimulus = [pops(GaussianNoiseStimulusParameter; SNR=[1.0, 1.0] .* SNR_scale), pops(SharpBumpStimulusParameter;
          strength = [stim_strength,stim_strength],
          width = [28.1, 28.1],
          time_windows = [[(0.0, 45.0)], [(0.0, 45.0)]])],
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
      algorithm=Tsit5()
  )
end

@testset "High-level tests" begin
    #@test based_on_example(; example_name="sigmoid_normal_fft", analyses=["radial_slice"], modifications=["iiS=1.7"]) != 0
	@test example() != 0
	ex = example()
	@test execute(ex) != 0
end
