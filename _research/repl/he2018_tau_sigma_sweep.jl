iterate_prototype("harris_ermentrout_monotonic",
                        ["tau=0.5:0.1:1.0", "sigma=0.5:0.1:1.6",
                         "stim_strength=[0.1]", "stim_width=700.0"],
                        data_root=data_root)
iterate_prototype("harris_ermentrout_blocking",
                        ["tau=0.5:0.1:1.0", "sigma=0.5:0.1:1.6",
                         "stim_strength=[0.1]", "stim_width=700.0"],
                        data_root=data_root)
