based_on_example(; example_name="reduced_line_dos_effectively_sigmoid", 
                        modifications=["Aee=40.0:60.0:250.0", "Aei=20.0:80.0:200.0",
                                       "Aie=15.0:80.0:140.0", "Aii=1.0:100.0:200.0",
                                       "firing_θE=6.0", "firing_θI=7.0",
                                       "blocking_θE=25.0", "blocking_θI=10.0",
                                       "n_lattice=512", "x_lattice=1400",
                                       "step_reduction=nothing"],
                        data_root=data_root)
