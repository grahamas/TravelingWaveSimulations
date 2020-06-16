include("example_helpers.jl")

examples_dict = Dict()

ABS_STOP=300.0
X_PROP = 0.9
examples_dict["reduced_line_dos_effectively_sigmoid"] = @EI_kw_example function example(N_ARR=1,N_CDT=1,P=2; SNR_scale=80.0, stop_time=ABS_STOP,
                                                     Aee=70.0, See=25.0,
                                                     Aii=2.0, Sii=27.0,
                                                     Aie=35.0, Sie=25.0,
                                                     Aei=70.0, Sei=27.0,
                                                     n=256, x=700.0, 
                                                     stim_strength=6.0,
                                                     stim_width=28.1,
                                                     stim_duration=7.0,other_opts=Dict())
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
      save_idxs=[IndexSubsampler((2,)), RightCutProportionFromValue((0.0,),(X_PROP,))],
      step_reduction = nothing,#(reduce_to_fronts(save_idxs, space), front_array_type),
      global_reduction = reduce_to_wave_properties,
      callback=terminate_when_E_fully_propagates(save_idxs, proportion_full=X_PROP, min_time=0.1), 
      other_opts...
  )
end

function get_example(example_name)
    return examples_dict[example_name]
end
