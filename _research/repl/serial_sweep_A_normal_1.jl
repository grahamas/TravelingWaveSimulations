based_on_example_serial(; example_name="reduced_line_dos_effectively_sigmoid", 
                        modifications=["Aee=40.0:12.0:250.0", "Aei=20.0:16.0:250.0",
                                       "Aie=15.0:16.0:250.0", "Aii=1.0:20.0:250.0",
                                       "firing_θE=6.0", "firing_θI=7.0",
                                       "blocking_θE=25.0", "blocking_θI=25.0",
                                       #"save_everystep=false"])
                                       "step_reduction=nothing"])
# FIXME step_reduction is broken so reverting to "nothing"
# This may be fine since I'm only saving the global reduction anyway
# However it will limit the maximum size simulation
