iterate_prototype("ring_blocking", 
                        ["See=14.0:24.0:120.0", "Sei=14.0:24.0:120.0",
                         "Sie=14.0:24.0:120.0", "Sii=14.0:24.0:120.0",
                         "Aii=50.0",
                         "θE=6.0", 
                         "firing_θI=7.0",
                         "blocking_θI=10.0",
                         "save_everystep=false",
						 "velocity_threshold=1e-6",
						 "n_traveling_frames_threshold=60"],
                        experiment_name="report2_S_sweep",
                        data_root=data_root)
# FIXME step_reduction is broken so reverting to "nothing"
# This may be fine since I'm only saving the global reduction anyway
# However it will limit the maximum size simulation
