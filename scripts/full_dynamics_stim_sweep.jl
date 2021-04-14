# FIXME WIP unconfirmed correspondance to HE2018 -- currently testing
using DrWatson
using TravelingWaveSimulations

let stim_strength_range = 0.1:0.1:2.0,
    αE = 0.4, αI = 0.7,
    firing_θI=0.2, blocking_θI=0.5,
    mod_strs = ["stim_strength=$(stim_strength_range)",
    "firing_θI=0.2", "blocking_θI=0.5",
    "αE=0.4", "αI=0.7"
    ],
    mods = (stim_strength=stim_strength_range,
        α=(αE,αI), firing_θI=firing_θI, blocking_θI=blocking_θI,
        global_reduction = TravelingWaveSimulations.reduce_to_min_propagation_cls,
        save_on=false, maxiters=1e3
    ),
    experiment_name = "stim_sweep"
;

if !@isdefined(data_root)
    data_root = datadir()
end

@time iterate_prototype("full_dynamics_monotonic",
                        mods, mod_strs;
                        data_root=data_root,
                        experiment_name=experiment_name,
                        progress=false
                    )
@time iterate_prototype("full_dynamics_blocking",
                        mods, mod_strs;
                        data_root=data_root,
                        experiment_name=experiment_name,
                        progress=false
                    )

end
