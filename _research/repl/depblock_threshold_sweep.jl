iterate_prototype("ring_blocking",
                        ["Aee=40.0:24.0:250.0", "Aei=20.0:32.0:250.0",
                         "Aie=15.0:32.0:250.0", "Aii=1.0:40.0:250.0",
                         "θE=6.0", "firing_θI=7.0",
                         "blocking_θI=9.0:1.0:20.0",
                         "step_reduction=nothing",
						 "velocity_threshold=1e-6",
						 "n_traveling_frames_threshold=60"],
                        experiment_name="depblock_A",
                        data_root=data_root)
