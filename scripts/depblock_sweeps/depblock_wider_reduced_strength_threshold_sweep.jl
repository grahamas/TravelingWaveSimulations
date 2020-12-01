iterate_prototype("ring_normed_blocking",
                        ["Aee=40.0:48.0:350.0", "Aei=20.0:48.0:350.0",
                         "Aie=15.0:48.0:250.0", "Aii=1.0:60.0:250.0",
                         "θE=6.0", "firing_θI=7.0",
                         "blocking_θI=9.0:3.0:30.0",
                         #"stim_strength=[2.0, 2.2, 2.5, 2.7, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0, 13.0, 16.0, 19.0]",
                         "stim_strength=2.0:3.0:20.0",
                         "step_reduction=nothing",
                         "velocity_threshold=1e-4",
                         "saveat=false",
                         "n_traveling_frames_threshold=30",
                         "x_lattice=2800.0", "n_lattice=512"],
                        experiment_name="wider_reduced_strength_depthreshold_A",
                        data_root=data_root,
                        max_sims_in_mem=5000)