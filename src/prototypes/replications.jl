replications_dict = Dict()

replications_dict["faulty_neuman_fft"] = (
                  N_ARR=1,N_CDT=1,P=2; 
                  SNR_scale=80.0, stop_time=ABS_STOP,
                  Aee=16.0, See=2.5,
                  Aii=4.0, Sii=2.7,
                  Aie=27.0, Sie=2.5,
                  Aei=18.2, Sei=2.7,
                  n_lattice=100, x_lattice=100.0, 
                  θE=2.6,
                  θI=8.0,
                  aE=1.2,
                  aI=1.0,
                  stim_strength=1.2,
                  stim_radius=14.0,
                  stim_duration=5.0,
                  common_baseline=0.1,
                  pop_names = ("E", "I"),
                  α = (1.5, 1.0),
                  β = (1.1, 1.1),
                  τ = (10.0, 18.0),
                  nonlinearity = PopulationActionsParameters(
                                     ZeroedSigmoidNonlinearity(a=aE, θ=θE), 
                                     SigmoidNonlinearity(a=aI, θ=θI)
                  ),
                  stimulus = pops(CircleStimulusParameter;
                      strength = [stim_strength, stim_strength],
                      radius = [stim_radius, stim_radius],
                      time_windows = [[(0.0, stim_duration)], [(0.0, stim_duration)]],
                      baseline=[common_baseline, common_baseline]
                  ),
                  connectivity = FFTParameter(pops(ExpAbsSumDecayingConnectivityParameter;
                      amplitude = [Aee -Aei;
                                   Aie -Aii],
                      spread = [(See,) (Sei,);
                                (Sie,) (Sii,)]
                     )
                  ),
                  space = PeriodicLattice{Float64,N_ARR}(; n_points=(n_lattice,), 
                                                           extent=(x_lattice,)),
                  tspan = (0.0,stop_time),
                  algorithm=Euler(),
                  dt=1.0,
                  save_idxs=nothing,
                  step_reduction=nothing,
                  global_reduction=identity,
                  callback=nothing,
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
        algorithm = algorithm,
        dt = dt,
        save_idxs = save_idxs,
        step_reduction = step_reduction,
        global_reduction = global_reduction,
        callback = callback,
        other_opts...
    )
end
replications_dict["faulty_neuman_fft"] = (
                  N_ARR=1,N_CDT=1,P=2; 
                  SNR_scale=80.0, stop_time=ABS_STOP,
                  Aee=16.0, See=2.5,
                  Aii=4.0, Sii=2.7,
                  Aie=27.0, Sie=2.5,
                  Aei=18.2, Sei=2.7,
                  n_lattice=100, x_lattice=100.0, 
                  θE=2.6,
                  θI=8.0,
                  aE=1.2,
                  aI=1.0,
                  stim_strength=1.2,
                  stim_radius=14.0,
                  stim_duration=5.0,
                  pop_names = ("E", "I"),
                  α = (1.5, 1.0),
                  β = (1.1, 1.1),
                  τ = (10.0, 18.0),
                  nonlinearity = PopulationActionsParameters(
                                     ZeroedSigmoidNonlinearity(a=aE, θ=θE), 
                                     SigmoidNonlinearity(a=aI, θ=θI)
                  ),
                  stimulus = pops(CircleStimulusParameter;
                      strength = [stim_strength, stim_strength],
                      radius = [stim_radius, stim_radius],
                      time_windows = [[(0.0, stim_duration)], [(0.0, stim_duration)]],
                      baseline=[common_baseline, common_baseline]
                  ),
                  connectivity = FFTParameter(pops(ExpAbsSumDecayingConnectivityParameter;
                      amplitude = [Aee -Aei;
                                   Aie -Aii],
                      spread = [(See,) (Sei,);
                                (Sie,) (Sii,)]
                     )
                  ),
                  space = PeriodicLattice{Float64,N_ARR}(; n_points=(n_lattice,), 
                                                           extent=(x_lattice,)),
                  tspan = (0.0,stop_time),
                  algorithm=Euler(),
                  dt=1.0,
                  save_idxs=nothing,
                  step_reduction=nothing,
                  global_reduction=nothing,
                  callback=nothing,
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
        algorithm = algorithm,
        dt = dt,
        save_idxs = save_idxs,
        step_reduction = step_reduction,
        global_reduction = global_reduction,
        callback = callback,
        other_opts...
    )
end
replications_dict["faulty_neuman_fft"] = (
                  N_ARR=1,N_CDT=1,P=2; 
                  SNR_scale=80.0, stop_time=ABS_STOP,
                  Aee=16.0, See=2.5,
                  Aii=4.0, Sii=2.7,
                  Aie=27.0, Sie=2.5,
                  Aei=18.2, Sei=2.7,
                  n_lattice=100, x_lattice=100.0, 
                  θE=2.6,
                  θI=8.0,
                  aE=1.2,
                  aI=1.0,
                  stim_strength=1.2,
                  stim_radius=14.0,
                  stim_duration=5.0,
                  pop_names = ("E", "I"),
                  α = (1.5, 1.0),
                  β = (1.1, 1.1),
                  τ = (10.0, 18.0),
                  nonlinearity = PopulationActionsParameters(
                                     ZeroedSigmoidNonlinearity(a=aE, θ=θE), 
                                     SigmoidNonlinearity(a=aI, θ=θI)
                  ),
                  stimulus = pops(CircleStimulusParameter;
                      strength = [stim_strength, stim_strength],
                      radius = [stim_radius, stim_radius],
                      time_windows = [[(0.0, stim_duration)], [(0.0, stim_duration)]],
                      baseline=[common_baseline, common_baseline]
                  ),
                  connectivity = FFTParameter(pops(ExpAbsSumDecayingConnectivityParameter;
                      amplitude = [Aee -Aei;
                                   Aie -Aii],
                      spread = [(See,) (Sei,);
                                (Sie,) (Sii,)]
                     )
                  ),
                  space = PeriodicLattice{Float64,N_ARR}(; n_points=(n_lattice,), 
                                                           extent=(x_lattice,)),
                  tspan = (0.0,stop_time),
                  algorithm=Euler(),
                  dt=1.0,
                  save_idxs=nothing,
                  step_reduction=nothing,
                  global_reduction=nothing,
                  callback=nothing,
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
        algorithm = algorithm,
        dt = dt,
        save_idxs = save_idxs,
        step_reduction = step_reduction,
        global_reduction = global_reduction,
        callback = callback,
        other_opts...
    )
end
replications_dict["faulty_neuman_fft"] = (
                  N_ARR=1,N_CDT=1,P=2; 
                  SNR_scale=80.0, stop_time=ABS_STOP,
                  Aee=16.0, See=2.5,
                  Aii=4.0, Sii=2.7,
                  Aie=27.0, Sie=2.5,
                  Aei=18.2, Sei=2.7,
                  n_lattice=100, x_lattice=100.0, 
                  θE=2.6,
                  θI=8.0,
                  aE=1.2,
                  aI=1.0,
                  stim_strength=1.2,
                  stim_radius=14.0,
                  stim_duration=5.0,
                  pop_names = ("E", "I"),
                  α = (1.5, 1.0),
                  β = (1.1, 1.1),
                  τ = (10.0, 18.0),
                  nonlinearity = PopulationActionsParameters(
                                     ZeroedSigmoidNonlinearity(a=aE, θ=θE), 
                                     SigmoidNonlinearity(a=aI, θ=θI)
                  ),
                  stimulus = pops(CircleStimulusParameter;
                      strength = [stim_strength, stim_strength],
                      radius = [stim_radius, stim_radius],
                      time_windows = [[(0.0, stim_duration)], [(0.0, stim_duration)]],
                      baseline=[common_baseline, common_baseline]
                  ),
                  connectivity = FFTParameter(pops(ExpAbsSumDecayingConnectivityParameter;
                      amplitude = [Aee -Aei;
                                   Aie -Aii],
                      spread = [(See,) (Sei,);
                                (Sie,) (Sii,)]
                     )
                  ),
                  space = PeriodicLattice{Float64,N_ARR}(; n_points=(n_lattice,), 
                                                           extent=(x_lattice,)),
                  tspan = (0.0,stop_time),
                  algorithm=Euler(),
                  dt=1.0,
                  save_idxs=nothing,
                  step_reduction=nothing,
                  global_reduction=nothing,
                  callback=nothing,
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
        algorithm = algorithm,
        dt = dt,
        save_idxs = save_idxs,
        step_reduction = step_reduction,
        global_reduction = global_reduction,
        callback = callback,
        other_opts...
    )
end
