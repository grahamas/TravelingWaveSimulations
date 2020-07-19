prototypes_dict = Dict()
function get_prototype(prototype_name)
    return prototypes_dict[prototype_name]
end

const ABS_STOP=300.0
const X_PROP = 0.9
prototypes_dict["dos_ring"] = (
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
                  dt=0.1,
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
        dt=dt,
        algorithm = algorithm,
        save_idxs = save_idxs,
        step_reduction = step_reduction,
        global_reduction = global_reduction,
        callback = callback,
        other_opts...
    )
end

prototypes_dict["dos_ring_absconn"] = (args...;
                  Aee=70.0, See=25.0,
                  Aii=2.0, Sii=27.0,
                  Aie=35.0, Sie=25.0,
                  Aei=70.0, Sei=27.0,
                  connectivity = FFTParameter(pops(ExpAbsSumDecayingConnectivityParameter;
                      amplitude = [Aee -Aei;
                                   Aie -Aii],
                      spread = [(See,) (Sei,);
                                (Sie,) (Sii,)]
                     )
                  ),
				 kwargs...) -> prototypes_dict["dos_ring"](args...;
				 	Aee=Aee, See=See,
					Aii=Aii, Sii=Sii,
					Aie=Aie, Sie=Sie,
					Aei=Aei, Sei=Sei,
					connectivity=connectivity,
					kwargs...)
					

prototypes_dict["oscillating_pulse"] = (
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
                  algorithm=Tsit5(),
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
        save_idxs = save_idxs,
        step_reduction = step_reduction,
        global_reduction = global_reduction,
        callback = callback,
        other_opts...
    )
end

prototypes_dict["propagating_torus"] = (
                  N_ARR=2,N_CDT=2,P=2; 
                  stop_time=ABS_STOP,
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
                  α = (1.0, 1.0),
                  β = (1.0, 1.0),
                  τ = (3.0, 3.0),
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
                      spread = [(See,See) (Sei,Sei);
                                (Sie,Sie) (Sii,Sii)]
                     )
                  ),
                  space = PeriodicLattice{Float64,N_ARR}(; n_points=(n_lattice,n_lattice), 
                                                           extent=(x_lattice,x_lattice)),
                  tspan = (0.0,stop_time),
                  dt = 0.1,
                  algorithm=Tsit5(),
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
        dt = dt,
        algorithm = algorithm,
        save_idxs = save_idxs,
        step_reduction = step_reduction,
        global_reduction = global_reduction,
        callback = callback,
        other_opts...
    )
end

prototypes_dict["orientation_torus"] = (
                  N_ARR=2,N_CDT=2,P=2; 
                  stop_time=ABS_STOP,
                  Aee=16.0, See=2.5,
                  Aii=4.0, Sii=2.7,
                  Aie=27.0, Sie=2.5,
                  Aei=18.2, Sei=2.7,
                  Atheta=16.0,
                  Stheta_x=25.0, Stheta_theta=pi/12,
                  n_lattice=100, x_lattice=100.0, 
                  n_theta=8, x_theta=2pi,
                  θE=2.6,
                  θI=8.0,
                  aE=1.2,
                  aI=1.0,
                  stim_strength=1.2,
                  stim_x=10.0, stim_theta=0.1,
                  stim_duration=5.0,
                  common_baseline=0.0,
                  pop_names = ("E", "I"),
                  α = (1.5, 1.0),
                  β = (1.1, 1.1),
                  τ = (10.0, 18.0),
                  nonlinearity = PopulationActionsParameters(
                                     ZeroedSigmoidNonlinearity(a=aE, θ=θE), 
                                     SigmoidNonlinearity(a=aI, θ=θI)
                  ),
                  stimulus = pops(RectangleStimulusParameter;
                      strength = [stim_strength, stim_strength],
                      radius = [(stim_x,stim_theta), (stim_x,stim_theta)],
                      time_windows = [[(0.0, stim_duration)], [(0.0, stim_duration)]],
                      baseline=[common_baseline, common_baseline]
                  ),
                  connectivity = FFTParameter(pops(ExpAbsSumDecayingConnectivityParameter;
                      amplitude = [Aee -Aei;
                                   Aie -Aii],
                      spread = [(See,See) (Sei,Sei);
                                (Sie,Sie) (Sii,Sii)]
                     )
                  ),
                  space = PeriodicLattice{Float64,N_ARR}(; n_points=(n_lattice,n_lattice), 
                                                           extent=(x_lattice,x_lattice)),
                  tspan = (0.0,stop_time),
                  dt = 0.1,
                  algorithm=Tsit5(),
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
        dt = dt,
        algorithm = algorithm,
        save_idxs = save_idxs,
        step_reduction = step_reduction,
        global_reduction = global_reduction,
        callback = callback,
        other_opts...
    )
end

prototypes_dict["harris_ermentrout_pure_sigmoid"] = (
                  N_ARR=1,N_CDT=1,P=2; 
                  SNR_scale=80.0, stop_time=ABS_STOP,
                  tau=0.68,
                  sigma=0.8,
                  Aee=1.0, See=1.0,
                  Aii=0.25, Sii=sigma,
                  Aie=1.0, Sie=1.0,
                  Aei=1.5, Sei=sigma,
                  n_lattice=512, x_lattice=1400.0, 
                  firing_θE=6.0,
                  firing_θI=11.4,
                  blocking_θE=30.0,
                  blocking_θI=30.0,
                  stim_strength=6.0,
                  stim_radius=14.0,
                  stim_duration=1.0,
                  pop_names = ("E", "I"),
                  α = (1.0, 1.0),
                  β = (1.0, 1.0),
                  τ = (1.0, tau),
                  nonlinearity = pops(UnrectifiedSigmoidNonlinearity;
                      θ = [0.125, 0.4],
                      a = [50.0, 50.0],
                  ),
                  stimulus = pops(CircleStimulusParameter;
                      strength = [stim_strength, stim_strength],
                      radius = [stim_radius, stim_radius],
                      time_windows = [[(0.0, stim_duration)], [(0.0, stim_duration)]],
                      baseline=[0.0, 0.0]
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
                  dt=0.1,
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
        dt=dt,
        algorithm = algorithm,
        save_idxs = save_idxs,
        step_reduction = step_reduction,
        global_reduction = global_reduction,
        callback = callback,
        other_opts...
    )
end
