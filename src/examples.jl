include("example_helpers.jl")

examples_dict = Dict()

const ABS_STOP=300.0
const X_PROP = 0.9
examples_dict["reduced_line_dos_effectively_sigmoid"] = (
                  N_ARR=1,N_CDT=1,P=2; 
                  SNR_scale=80.0, stop_time=ABS_STOP,
                  Aee=70.0, See=25.0,
                  Aii=2.0, Sii=27.0,
                  Aie=35.0, Sie=25.0,
                  Aei=70.0, Sei=27.0,
                  n_lattice=512, x_lattice=1400.0, 
                  firing_θE=6.0,
                  firing_θI=11.4,
                  blocking_θE=30.0,
                  blocking_θI=30.0,
                  stim_strength=6.0,
                  stim_radius=14.0,
                  stim_duration=7.0,
                  pop_names = ("E", "I"),
                  α = (1.0, 1.0),
                  β = (1.0, 1.0),
                  τ = (3.0, 3.0),
                  nonlinearity = pops(DifferenceOfSigmoids;
                      firing_θ = [firing_θE, firing_θI],
                      firing_a = [1.2, 1.0],
                      blocking_θ = [blocking_θE, blocking_θI],
                      blocking_a = [1.2, 1.0]
                  ),
                  stimulus = pops(CircleStimulusParameter;
                      strength = [stim_strength, stim_strength],
                      radius = [stim_radius, stim_radius],
                      time_windows = [[(0.0, stim_duration)], [(0.0, stim_duration)]],
                      baseline=[0.0, 0.0]
                  ),
                  connectivity = FFTParameter(pops(GaussianConnectivityParameter;
                      amplitude = [Aee -Aei;
                                   Aie -Aii],
                      spread = [(See,) (Sei,);
                                (Sie,) (Sii,)]
                     )
                  ),
                  space = PeriodicLattice{Float64,N_ARR}(; n_points=(n_lattice,), 
                                                           extent=(x_lattice,)),
                  tspan = (0.0,stop_time),
                  dt = 0.1,
                  algorithm=Tsit5(),
                  save_idxs=nothing,
                  step_reduction=(reduce_to_fronts(save_idxs, space), front_array_type),
                  global_reduction = reduce_to_wave_properties,
                  callback=terminate_when_E_fully_propagates(save_idxs, 
                                                             proportion_full=X_PROP, 
                                                             min_time=0.1), 
                  other_opts...
            ) -> begin
    Simulation(
        WCMSpatial{N_CDT,P}(;
            pop_names = pop_names,
            α = α,
            β = β,
            τ = τ,
            nonlinearity = nonlinearity,
            stimulus = stimulus,
            connectivity = connectivity
           );
        space = space,
        tspan = tspan,
        dt = dt,
        algorithm = algorithm,
        save_idxs = save_idxs,
        step_reduction = step_reduction,
        global_reduction = global_reduction,
        callback = callback,
        other_opts...
    )
end

examples_dict["simple_square"] = @EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; stop_time=ABS_STOP,
                                                     Aee=16.0, See=2.5,
                                                     Aii=4.0, Sii=2.7,
                                                     Aie=27.0, Sie=2.5,
                                                     Aei=18.2, Sei=2.7,
                                                     n=100, x=100,
                                                     stim_strength=1.2,
                                                     common_baseline = 0.1,
                                                     stim_radius=14.0,
                                                     other_opts=Dict()
                                                    )
  simulation = Simulation(
                          WCMSpatial{N_CDT,P}(;
      pop_names = ("E", "I"),
      α = (1.5, 1.0),
      β = (1.1, 1.1),
      τ = (10.0, 18.0),
      nonlinearity = PopulationActionsParameters(NeuralModels.ZeroedSigmoidNonlinearity(a = 1.2, θ=2.6), SigmoidNonlinearity(a = 1.0, θ = 8.0)),
      stimulus =  pops(CircleStimulusParameter;
          strength = [stim_strength,stim_strength],
          radius = [stim_radius, stim_radius],
          time_windows = [[(0.0, 5.0)], [(0.0, 5.0)]],
          baseline = [common_baseline, common_baseline]),
      connectivity = FFTParameter(#pops([
          pops(ExpAbsSumDecayingConnectivityParameter;
            amplitude = [Aee -Aei;
                       Aie -Aii],
            spread = [(See, See) (Sei, Sei);
                    (Sie, Sie) (Sii, Sii)]
          )#,  
          #pops([NoOpParameter{typeof(Atheta)}() ExpAbsSumDecayingConnectivityParameter(; amplitude=Atheta, spread=(Stheta_x, Stheta_theta));
           #     NoOpParameter{typeof(Atheta)}() NoOpParameter{typeof(Atheta)}()])
           # ])
        )
      );
      space = PeriodicLattice{Float64,N_ARR}(; n_points=(n,n), extent=(x,x)),
      save_idxs=nothing,
      algorithm=nothing,
      tspan = (0.0,stop_time),
      other_opts...
     )
end

examples_dict["orientation"] = @EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; stop_time=ABS_STOP,
                                                     Aee=16.0, See=2.5,
                                                     Aii=4.0, Sii=2.7,
                                                     Aie=27.0, Sie=2.5,
                                                     Aei=18.2, Sei=2.7,
                                                     Atheta=16.0,
                                                     Stheta_x=25.0, Stheta_theta=pi/12,
                                                     n_space=100, x_space=100.0, 
                                                     n_theta=8, x_theta=2pi,
                                                     stim_strength=1.2,
                                                     common_baseline = 0.0,
                                                     stim_x=10.0, stim_theta=0.1,
                                                     other_opts=Dict()
                                                    )
  simulation = Simulation(
                          WCMSpatial{N_CDT,P}(;
      pop_names = ("E", "I"),
      α = (1.5, 1.0),
      β = (1.1, 1.1),
      τ = (10.0, 18.0),
      nonlinearity = PopulationActionsParameters(NeuralModels.ZeroedSigmoidNonlinearity(a = 1.2, θ=2.6), SigmoidNonlinearity(a = 1.0, θ = 8.0)),
      stimulus =  pops(RectangleStimulusParameter;
          strength = [stim_strength,stim_strength],
          widths = [(stim_x, stim_theta), (stim_x, stim_theta)],
          time_windows = [[(0.0, 5.0)], [(0.0, 5.0)]],
          baseline = [common_baseline, common_baseline]),
      connectivity = FFTParameter(pops([
          pops(ExpAbsSumDecayingConnectivityParameter;
            amplitude = [Aee -Aei;
                       Aie -Aii],
            spread = [(See, Stheta_x) (Sei, Stheta_x);
                    (Sie, Stheta_x) (Sii, Stheta_x)]
          ),  
          pops([NoOpParameter{typeof(Atheta)}() ExpAbsSumDecayingConnectivityParameter(; amplitude=Atheta, spread=(Stheta_x, Stheta_theta));
                NoOpParameter{typeof(Atheta)}() NoOpParameter{typeof(Atheta)}()])
         ])
                                 )
      );
      space = PeriodicLattice{Float64,N_ARR}(; n_points=(n_space,n_theta), extent=(x_space,x_theta)),
      save_idxs=nothing,
      tspan = (0.0,stop_time),
      algorithm=nothing,
      other_opts...
     )
end

examples_dict["oscillating_pulse"] = @EI_kw_example function example(N_ARR=1,N_CDT=1,P=2; stop_time=ABS_STOP,
                                                     Aee=16.0, See=2.5,
                                                     Aii=4.0, Sii=2.7,
                                                     Aie=27.0, Sie=2.5,
                                                     Aei=18.2, Sei=2.7,
                                                     n=100, x=100.0, 
                                                     stim_strength=1.2,
                                                     common_baseline = 0.1,
                                                     stim_radius=14.0,
                                                     other_opts=Dict()
                                                    )
  simulation = Simulation(
                          WCMSpatial{N_CDT,P}(;
      pop_names = ("E", "I"),
      α = (1.5, 1.0),
      β = (1.1, 1.1),
      τ = (10.0, 18.0),
      nonlinearity = PopulationActionsParameters(NeuralModels.ZeroedSigmoidNonlinearity(a = 1.2, θ=2.6), SigmoidNonlinearity(a = 1.0, θ = 8.0)),
      stimulus =  pops(CircleStimulusParameter;
          strength = [stim_strength,stim_strength],
          radius = [stim_radius, stim_radius],
          time_windows = [[(0.0, 5.0)], [(0.0, 5.0)]],
          baseline = [common_baseline, common_baseline]),
      connectivity = FFTParameter(pops(ExpAbsSumDecayingConnectivityParameter;
          amplitude = [Aee -Aei;
                       Aie -Aii],
          spread = [(See,) (Sei,);
                    (Sie,) (Sii,)]
          ))
      );
      space = PeriodicLattice{Float64,N_ARR}(; n_points=(n,), extent=(x,)),
      save_idxs=nothing,
      tspan = (0.0,stop_time),
      other_opts...
  )
end



#### REPLICATIONS ####

examples_dict["FAULTY_neuman_fft"] = @EI_kw_example function example(N_ARR=1,N_CDT=1,P=2; stop_time=ABS_STOP,
                                                     Aee=16.0, See=2.5,
                                                     Aii=4.0, Sii=2.7,
                                                     Aie=27.0, Sie=2.5,
                                                     Aei=18.2, Sei=2.7,
                                                     n=100, x=100.0, 
                                                     stim_strength=1.2,
                                                     common_baseline = 0.1,
                                                     stim_radius=14.0,
                                                     other_opts=Dict()
                                                    )
  simulation = Simulation(
                          WCMSpatial{N_CDT,P}(;
      pop_names = ("E", "I"),
      α = (1.5, 1.0),
      β = (1.1, 1.1),
      τ = (10.0, 18.0),
      nonlinearity = PopulationActionsParameters(NeuralModels.ZeroedSigmoidNonlinearity(a = 1.2, θ=2.6), SigmoidNonlinearity(a = 1.0, θ = 8.0)),
      stimulus =  pops(CircleStimulusParameter;
          strength = [stim_strength,stim_strength],
          radius = [stim_radius, stim_radius],
          time_windows = [[(0.0, 5.0)], [(0.0, 5.0)]],
          baseline = [common_baseline, common_baseline]),
      connectivity = FFTParameter(pops(ExpAbsSumDecayingConnectivityParameter;
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
      other_opts...
  )
end

examples_dict["orientation_NOFFT"] = @EI_kw_example function example(N_ARR=2,N_CDT=2,P=2; stop_time=ABS_STOP,
                                                     Aee=16.0, See=2.5,
                                                     Aii=4.0, Sii=2.7,
                                                     Aie=27.0, Sie=2.5,
                                                     Aei=18.2, Sei=2.7,
                                                     Atheta=16.0,
                                                     Stheta_x=25.0, Stheta_theta=pi/12,
                                                     n_space=50, x_space=100.0, 
                                                     n_theta=8, x_theta=2pi,
                                                     stim_strength=1.2,
                                                     common_baseline = 0.0,
                                                     stim_x=10.0, stim_theta=0.1,
                                                     other_opts=Dict()
                                                    )
  simulation = Simulation(
                          WCMSpatial{N_CDT,P}(;
      pop_names = ("E", "I"),
      α = (1.5, 1.0),
      β = (1.1, 1.1),
      τ = (10.0, 18.0),
      nonlinearity = PopulationActionsParameters(NeuralModels.ZeroedSigmoidNonlinearity(a = 1.2, θ=2.6), SigmoidNonlinearity(a = 1.0, θ = 8.0)),
      stimulus =  pops(RectangleStimulusParameter;
          strength = [stim_strength,stim_strength],
          widths = [(stim_x, stim_theta), (stim_x, stim_theta)],
          time_windows = [[(0.0, 5.0)], [(0.0, 5.0)]],
          baseline = [common_baseline, common_baseline]),
      connectivity = pops([
          pops(ExpAbsSumDecayingConnectivityParameter;
            amplitude = [Aee -Aei;
                       Aie -Aii],
            spread = [(See, Stheta_x) (Sei, Stheta_x);
                    (Sie, Stheta_x) (Sii, Stheta_x)]
          ),  
          pops([NoOpParameter{typeof(Atheta)}() ExpAbsSumDecayingConnectivityParameter(; amplitude=Atheta, spread=(Stheta_x, Stheta_theta));
                NoOpParameter{typeof(Atheta)}() NoOpParameter{typeof(Atheta)}()])
         ])
      );
      space = PeriodicLattice{Float64,N_ARR}(; n_points=(n_space,n_theta), extent=(x_space,x_theta)),
      save_idxs=nothing,
      tspan = (0.0,stop_time),
      algorithm=nothing,
      other_opts...
     )
end

function get_example(example_name)
    return examples_dict[example_name]
end
