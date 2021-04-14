# FIXME WIP unconfirmed correspondance to HE2018 -- currently testing
using DrWatson
using TravelingWaveSimulations

let mod_strs = ["tau=0.5:0.3:1.0", "sigma=0.5:0.3:1.6",
    "stim_strength=[0.1]", "stim_width=700.0"],
    mods = (tau=0.5:0.3:1.0, sigma=0.5:0.3:1.6,
    stim_strength=[0.1], stim_width=700.0,
    global_reduction = TravelingWaveSimulations.reduce_to_min_propagation_cls)
;

if !@isdefined(data_root)
    data_root = datadir()
end

iterate_prototype("full_dynamics_monotonic",
                        mods, mod_strs,
                        data_root=data_root)
iterate_prototype("full_dynamics_blocking",
                        mods, mod_strs,
                        data_root=data_root)

end
