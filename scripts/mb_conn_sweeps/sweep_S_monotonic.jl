iterate_prototype("ring_monotonic", 
                        ["See=14.0:12.0:120.0", "Sei=14.0:12.0:120.0",
                         "Sie=14.0:12.0:120.0", "Sii=14.0:12.0:120.0",
                         "θE=6.0", 
                         "θI=7.0",
                         "Aii=50.0",
                         "Aei=70.0",
                         "Aie=70.0",
                         "Aee=70.0",
                         "save_everystep=false",
						 "velocity_threshold=1e-6",
						 "n_traveling_frames_threshold=60"],
                        experiment_name="report2_S_sweep",
                        data_root=data_root)
# FIXME step_reduction is broken so reverting to "nothing"
# This may be fine since I'm only saving the global reduction anyway
# However it will limit the maximum size simulation
