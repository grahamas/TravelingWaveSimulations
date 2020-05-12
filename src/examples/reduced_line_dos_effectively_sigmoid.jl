# -*- coding: utf-8 -*-
ABS_STOP=300.0
TravelingWaveSimulations.@EI_kw_example function example(N_ARR=1,N_CDT=1,P=2; SNR_scale=80.0, stop_time=ABS_STOP,
                                                     Aee=70.0, See=25.0,
                                                     Aii=2.0, Sii=27.0,
                                                     Aie=35.0, Sie=25.0,
                                                     Aei=70.0, Sei=27.0,
                                                     n=256, x=700.0, 
                                                     stim_strength=6.0,
                                                     stim_width=28.1,
                                                     stim_duration=7.0,other_opts=Dict(:saveat=>[0.0]))
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
      save_idxs=[IndexSubsampler((2,)), RightCutFromValue((0.0,))],
      step_reduction = ((u, t, integrator) -> begin
        # Need to get x from integrator (or simulation)
        sub_idx = save_idxs === nothing ? CartesianIndices(u) : integrator.opts.save_idxs
        sub_u = u[sub_idx]
        sub_x = [x[1] for x in space.arr[population(sub_idx,1)]]
        fronts = TravelingWaveSimulations.substantial_fronts(sub_u, sub_x)
        return fronts
      end, Array{Wavefront{Float64,Float64,Value{Float64,Float64}},1}),
      global_reduction = (data_named_tuple) -> (wave_properties=get_wave_properties(data_named_tuple),),
      callback=DiscreteCallback(if !(save_idxs === nothing)
            cb_E_is_fully_propagated_with_save_idxs
        else
            cb_E_is_fully_propagated
        end, terminate!), 
      other_opts...
  )
end

