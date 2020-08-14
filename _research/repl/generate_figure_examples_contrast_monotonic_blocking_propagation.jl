include("figure_contrast_monotonic_blocking_propagation.jl")

# Definitely propagating
save_figure_example_contrast_monotonic_blocking_all((:Aei, :Aee), (:Aie, :Aii), [
                          (Aee=200.0, Aei=50.0, Aii=200, Aie=50.0),
                          (Aee=150.0, Aei=150.0, Aii=50, Aie=115.0),
                          (Aee=50.0, Aei=200.0, Aii=200, Aie=50.0)
                         ],  
                            fcmb_monotonic_A_fpath, fcmb_blocking_A_fpath, 
                            :has_propagation, "contrast_monotonic_blocking_A";
                            primary_squish_origin = SA[0.0, 244.0],
                            primary_squish_line = SA[-0.8125, 1.],
                            secondary_squish_origin = SA[0.0, 244.0],
                            secondary_squish_line = SA[-0.8125, 1.]
                             )

save_figure_example_contrast_monotonic_blocking_all((:Sei, :See), (:Sie, :Sii), [
                          (See=14.0, Sei=110.0, Sii=14., Sie=110.0),
                          (See=14.0, Sei=14.0, Sii=110., Sie=14.0),
                          #(See=30.0, Sei=50.0, Sii=50., Sie=90.0),
                          (See=50.0, Sei=40.0, Sii=40., Sie=50.0)
                         ],  
                            fcmb_monotonic_S_fpath, fcmb_blocking_S_fpath, 
                            :has_propagation, "contrast_monotonic_blocking_S",
                            primary_squish_origin = SA[0.0, 110.0],
                            primary_squish_line = SA[1., -2.0],
                            secondary_squish_origin = SA[0.0, 110.0],
                            secondary_squish_line = SA[1., -2.0]
                           )
#save_figure_contrast_monotonic_blocking((:Aie, :Aee), fcmb_monotonic_A_fpath, fcmb_blocking_A_fpath, :has_propagation, "contrast_monotonic_blocking_A")
#save_figure_contrast_monotonic_blocking((:Aei, :Aie), fcmb_monotonic_A_fpath, fcmb_blocking_A_fpath, :has_propagation, "contrast_monotonic_blocking_A")
#save_figure_contrast_monotonic_blocking((:Aii, :Aee), fcmb_monotonic_A_fpath, fcmb_blocking_A_fpath, :has_propagation, "contrast_monotonic_blocking_A")
#save_figure_contrast_monotonic_blocking((:Aii, :Aie), fcmb_monotonic_A_fpath, fcmb_blocking_A_fpath, :has_propagation, "contrast_monotonic_blocking_A")
#save_figure_contrast_monotonic_blocking((:Aii, :Aei), fcmb_monotonic_A_fpath, fcmb_blocking_A_fpath, :has_propagation, "contrast_monotonic_blocking_A")
#
#save_figure_contrast_monotonic_blocking((:Sei, :See), fcmb_monotonic_S_fpath, fcmb_blocking_S_fpath, :has_propagation, "contrast_monotonic_blocking_S")
#save_figure_contrast_monotonic_blocking((:Sie, :See), fcmb_monotonic_S_fpath, fcmb_blocking_S_fpath, :has_propagation, "contrast_monotonic_blocking_S")
#save_figure_contrast_monotonic_blocking((:Sei, :Sie), fcmb_monotonic_S_fpath, fcmb_blocking_S_fpath, :has_propagation, "contrast_monotonic_blocking_S")
#save_figure_contrast_monotonic_blocking((:Sii, :See), fcmb_monotonic_S_fpath, fcmb_blocking_S_fpath, :has_propagation, "contrast_monotonic_blocking_S")
#save_figure_contrast_monotonic_blocking((:Sii, :Sie), fcmb_monotonic_S_fpath, fcmb_blocking_S_fpath, :has_propagation, "contrast_monotonic_blocking_S")
#save_figure_contrast_monotonic_blocking((:Sii, :Sei), fcmb_monotonic_S_fpath, fcmb_blocking_S_fpath, :has_propagation, "contrast_monotonic_blocking_S")
