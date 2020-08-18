iterate_prototype("dos_ring_absconn",
                        ["Aee=40.0:12.0:250.0", "Aei=20.0:16.0:250.0",
                         "Aie=15.0:16.0:250.0", "Aii=1.0:20.0:250.0",
                         "firing_θE=6.0", "firing_θI=7.0",
                         "blocking_θE=25.0", "blocking_θI=10.0",
                         "step_reduction=nothing",
						 "velocity_threshold=1e-6",
						 "n_traveling_frames_threshold=10"],
                        data_root=data_root)
