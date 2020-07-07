based_on_example(; example_name="reduced_line_dos_effectively_sigmoid", 
                        modifications=["See=14.0:16.0:320.0", "Sei=14.0:16.0:320.0",
                                       "Sie=14.0:16.0:320.0", "Sii=14.0:24.0:600.0",
                                       "firing_θE=6.0", "firing_θI=7.0",
                                       "blocking_θE=25.0", "blocking_θI=10.0",
                                       "n_lattice=512", "x_lattice=1400",
                                       "step_reduction=nothing"],
                        data_root=data_root)
# FIXME step_reduction is broken so reverting to "nothing"
# This may be fine since I'm only saving the global reduction anyway
# However it will limit the maximum size simulation
