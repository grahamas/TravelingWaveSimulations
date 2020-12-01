iterate_prototype("ring_normed_blocking",
                        ["Aee=40.0:48.0:250.0", "Aei=20.0:64.0:250.0",
                         "Aie=15.0:64.0:250.0", "Aii=1.0:80.0:250.0",
                         "θE=6.0", "firing_θI=7.0",
                         "blocking_θI=9.0:1.0:20.0",
                         "x_lattice=2800.0", "n_lattice=512"],
                        experiment_name="reduced_depthreshold_A",
                        data_root=data_root)
