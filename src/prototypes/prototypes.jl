prototypes_dict = Dict()
function get_prototype(prototype_name)
    return prototypes_dict[prototype_name]
end

const ABS_STOP=300.0
const X_PROP = 0.9
prototypes_dict["ring_monotonic"] = (
                  N_ARR=1,N_CDT=1,P=2; 
                  SNR_scale=80.0, stop_time=ABS_STOP,
                  Aee=70.0, See=25.0,
                  Aii=2.0, Sii=27.0,
                  Aie=35.0, Sie=25.0,
                  Aei=70.0, Sei=27.0,
                  n_lattice=512, x_lattice=1400.0, 
                  aE=1.2, θE=6.0,
                  aI=1.0, θI=11.4,
                  stim_strength=6.0,
                  stim_radius=14.0,
                  stim_duration=7.0,
                  pop_names = ("E", "I"),
                  min_dist_for_propagation=x_lattice * 0.4,
                  const_jitter = (x_lattice / n_lattice) * 3,
                  vel_jitter = 1.5,
                  slope_min=1e-4,
                  α = (1.0, 1.0),
                  β = (1.0, 1.0),
                  τ = (3.0, 3.0),
                  nonlinearity = pops(RectifiedSigmoidNonlinearity;
                      θ = [θE, θI],
                      a = [aE, aI]
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
                  algorithm=Tsit5(),
                  save_idxs=nothing,
                  dt=0.1,
                  callback = (
                      is_propagated,
                    (min_dist_for_propagation=min_dist_for_propagation,
                     slope_min=slope_min, const_jitter=const_jitter,
                     vel_jitter=vel_jitter)
                  ),
                  global_reduction = already_reduced_to_min_propagation_cls,
                  save_on=false,
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
        global_reduction = global_reduction,
        callback = callback,
        save_on=save_on,
        other_opts...
    )
end

prototypes_dict["ring_blocking"] = (args...;
                  blocking_θI=10.0,
                  blocking_aI=1.0,
                  firing_θI=7.0,
                  firing_aI=1.0,
                  θE=6.0,
                  aE=1.2,
                  nonlinearity = PopulationActionsParameters(
                                     RectifiedSigmoidNonlinearity(a=aE, θ=θE),
                                     DifferenceOfSigmoidsParameter(
                                        firing_θ = firing_θI,
                                        firing_a = firing_aI,
                                        blocking_θ = blocking_θI,
                                        blocking_a = blocking_aI
                                     )
                  ),
                  kwargs...) -> begin
prototypes_dict["ring_monotonic"](args...;
                                    θE = θE,
                                    aE = aE,
                                    nonlinearity=nonlinearity,
                                    kwargs...)
end   

prototypes_dict["ring_normed_blocking"] = (args...;
                  blocking_θI=10.0,
                  blocking_aI=1.0,
                  firing_θI=7.0,
                  firing_aI=1.0,
                  θE=6.0,
                  aE=1.2,
                  nonlinearity = PopulationActionsParameters(
                                     RectifiedSigmoidNonlinearity(a=aE, θ=θE),
                                     NormedDifferenceOfSigmoidsParameter(
                                        firing_θ = firing_θI,
                                        firing_a = firing_aI,
                                        blocking_θ = blocking_θI,
                                        blocking_a = blocking_aI
                                     )
                  ),
                  kwargs...) -> begin
prototypes_dict["ring_monotonic"](args...;
                                    θE = θE,
                                    aE = aE,
                                    nonlinearity=nonlinearity,
                                    kwargs...)
end    

prototypes_dict["harris_ermentrout"] = (N_ARR=1,N_CDT=1,P=2; 
                  stop_time=ABS_STOP,
                  Aee=1., See=25.0,
                  Aii=0.25, Sii=27.0,
                  Aie=1., Sie=25.0,
                  Aei=1.5, Sei=27.0,
                  n_lattice=512, x_lattice=1400.0, 
                  aE=50.0, θE=0.125,
                  aI=50.0, θI=0.4,
                  stim_strength=0.5,
                  stim_radius=14.0,
                  stim_duration=7.0,
                  pop_names = ("E", "I"),
                  min_dist_for_propagation=x_lattice * 0.4,
                  const_jitter = (x_lattice / n_lattice) * 3,
                  vel_jitter = 1.5,
                  slope_min=1e-4,
                  α = (1.0, 1.0),
                  β = (1.0, 1.0),
                  τ = (1.0, 0.4),
                  nonlinearity = pops(SimpleSigmoidNonlinearity;
                      θ = [θE, θI],
                      a = [aE, aI]
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
                  algorithm=Tsit5(),
                  save_idxs=nothing,
                  callback = (
                      is_propagated,
                    (min_dist_for_propagation=min_dist_for_propagation,
                     slope_min=slope_min, const_jitter=const_jitter,
                     vel_jitter=vel_jitter)
                  ),
                  global_reduction = already_reduced_to_min_propagation_cls,
                  save_on=false,
                  other_opts...
            ) -> begin
    @show α
    Simulation(
        HarrisErmentrout2018{N_CDT,P}(;
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
        global_reduction = global_reduction,
        callback = callback,
        save_on=save_on,
        other_opts...
    )
end

prototypes_dict["full_dynamics_monotonic"] = (N_ARR=1,N_CDT=1,P=2; 
                  stop_time=ABS_STOP,
                  Aee=1., See=25.0,
                  Aii=0.25, Sii=27.0,
                  Aie=1., Sie=25.0,
                  Aei=1.5, Sei=27.0,
                  n_lattice=512, x_lattice=1400.0, 
                  aE=50.0, θE=0.125,
                  firing_θI = 0.2,
                  firing_aI = 0.2,
                  aI=firing_aI, θI=firing_θI,
                  stim_strength=0.5,
                  stim_radius=14.0,
                  stim_duration=7.0,
                  pop_names = ("E", "I"),
                  min_dist_for_propagation=x_lattice * 0.4,
                  const_jitter = (x_lattice / n_lattice) * 3,
                  vel_jitter = 1.5,
                  slope_min=1e-4,
                  αE = 0.4, αI = 0.7,
                  α=(αE, αI),
                  β = (1.0, 1.0),
                  τ = (1.0, 0.4),
                  nonlinearity = pops(SimpleSigmoidNonlinearity;
                      θ = [θE, θI],
                      a = [aE, aI]
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
                  algorithm=Tsit5(),
                  save_idxs=nothing,
                  callback = [(
                      is_propagated,
                    (min_dist_for_propagation=min_dist_for_propagation,
                     slope_min=slope_min, const_jitter=const_jitter,
                     vel_jitter=vel_jitter)
                  ),
                  (
                      is_spread,
                      (max_spread_proportion=0.75,
                      min_spread_time=8,
                      max_spread_value=1e-3)
                  )],
                  global_reduction=identity,
                  #global_reduction = already_reduced_to_min_propagation_cls,
                  save_on=false,
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
        global_reduction = global_reduction,
        callback = callback,
        save_on=save_on,
        other_opts...
    )
end

prototypes_dict["full_dynamics_blocking"] = (args...; 
            α=(0.4, 0.7),
            blocking_θI=0.5,
            blocking_aI=50.0,
            firing_θI=0.2,
            firing_aI=50.0,
            aE=50.0, θE=0.125,
            nonlinearity = PopulationActionsParameters(
                                RectifiedSigmoidNonlinearity(a=aE, θ=θE),
                                NormedDifferenceOfSigmoidsParameter(
                                    firing_θ = firing_θI,
                                    firing_a = firing_aI,
                                    blocking_θ = blocking_θI,
                                    blocking_a = blocking_aI
                                )
                ),
            kwargs...) -> begin
    prototypes_dict["full_dynamics_monotonic"](args...; 
        α, blocking_aI, blocking_θI, firing_aI, firing_θI,
        aE, θE,
        nonlinearity=nonlinearity, kwargs...)
end

prototypes_dict["full_dynamics_blocking_erf"] = (args...; 
            α=(0.4, 0.7),
            blocking_θI=0.5,
            blocking_aI=50.0,
            firing_θI=0.2,
            firing_aI=50.0,
            aE=50.0, θE=0.125,
            nonlinearity = PopulationActionsParameters(
                                ErfNonlinearity(a=aE, θ=θE),
                                DifferenceOfErfsParameter(
                                    firing_θ = firing_θI,
                                    firing_a = firing_aI,
                                    blocking_θ = blocking_θI,
                                    blocking_a = blocking_aI
                                )
                ),
            kwargs...) -> begin
    prototypes_dict["full_dynamics_monotonic"](args...; 
        α, blocking_aI, blocking_θI, firing_aI, firing_θI,
        aE, θE,
        nonlinearity=nonlinearity, kwargs...)
end

prototypes_dict["full_dynamics_monotonic_erf"] = (args...; 
            α=(0.4, 0.7),
            θI=0.2,
            aI=50.0,
            aE=50.0, θE=0.125,
            nonlinearity = pops(ErfNonlinearity;
                      θ = [θE, θI],
                      a = [aE, aI]
                  ),
            kwargs...) -> begin
    prototypes_dict["full_dynamics_monotonic"](args...; 
        α, blocking_aI, blocking_θI, firing_aI, firing_θI,
        aE, θE,
        nonlinearity=nonlinearity, kwargs...)
end


					

prototypes_dict["oscillating_pulse_monotonic"] = (
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
                                     RectifiedZeroedSigmoidNonlinearity(a=aE, θ=θE), 
                                     RectifiedSigmoidNonlinearity(a=aI, θ=θI)
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

prototypes_dict["propagating_torus_monotonic"] = (
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
                                     RectifiedZeroedSigmoidNonlinearity(a=aE, θ=θE), 
                                     RectifiedSigmoidNonlinearity(a=aI, θ=θI)
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

prototypes_dict["orientation_torus_monotonic"] = (
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
                                     RectifiedZeroedSigmoidNonlinearity(a=aE, θ=θE), 
                                     RectifiedSigmoidNonlinearity(a=aI, θ=θI)
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

prototypes_dict["harris_ermentrout_rectified"] = (
                  args...;
                  aE=50.0, θE=0.125,
                  aI=50.0, θI=0.4,
                  nonlinearity = pops(RectifiedSigmoidNonlinearity;
                      θ = [θE, θI],
                      a = [aE, aI],
                  ),
                  kwargs...
            ) -> begin
    prototypes_dict["harris_ermentrout"](args...; nonlinearity=nonlinearity, kwargs...)
end
