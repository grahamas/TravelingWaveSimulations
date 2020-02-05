# -*- coding: utf-8 -*-
ABS_STOP = 100.0
@EI_kw_example function example(N_ARR=1,N_CDT=1,P=2; stop_time=ABS_STOP,
                                                     Aee=16.0, See=2.5,
                                                     Aii=4.0, Sii=2.7,
                                                     Aie=27.0, Sie=2.5,
                                                     Aei=18.2, Sei=2.7,
                                                     n=101, x=100.0, 
                                                     stim_strength=1.2,
                                                     common_baseline = 0.1,
                                                     stim_width=28.1)
  simulation = Simulation(
    WCMSpatial(;
      pop_names = ("E", "I"),
      α = (1.5, 1.0),
      β = (1.1, 1.1),
      τ = (10.0, 18.0),
      nonlinearity = PopulationActionsParameters(SigmoidNonlinearity(a = 1.2, θ=2.6), NeuralModels.UnzeroedSigmoidNonlinearity(a = 1.0, θ = 8.0)),
      stimulus =  pops(SharpBumpStimulusParameter;
          strength = [stim_strength,stim_strength],
          width = [stim_width, stim_width],
          time_windows = [[(0.0, 5.0)], [(0.0, 5.0)]],
          baseline = [common_baseline, common_baseline]),
      connectivity = pops(ExpSumAbsDecayingConnectivityParameter;
          amplitude = [Aee -Aei;
                       Aie -Aii],
          spread = [(See,) (Sei,);
                    (Sie,) (Sii,)]
          )
      );
      space = PeriodicLattice{Float64,N_ARR}(; n_points=(n,), extent=(x,)),
      save_idxs=nothing,
      tspan = (0.0,stop_time),
      dt = 1.0,
      algorithm=Euler(),
#       callback=DiscreteCallback(if !(save_idxs === nothing)
#         (u,t,integrator) -> begin
#                     sub_u = u[integrator.opts.save_idxs];
#                     (all(sub_u .≈ 0.0) || (sub_u[end] > 0.01)) && t > 5
#                 end
#     else
#         (u,t,integrator) -> begin
#                     pop = population(u,1)
#                     (all(u .≈ 0.0) || (sum(pop[:,end]) / size(pop,1) > 0.01)) && t > 5
#             end
#     end, terminate!)
  )
end
