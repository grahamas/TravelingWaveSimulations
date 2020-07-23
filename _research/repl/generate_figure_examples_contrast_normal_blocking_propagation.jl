include("figure_contrast_normal_blocking_propagation.jl")

# Definitely propagating
save_figure_example_contrast_monotonic_blocking((:Aei, :Aee), [
                          (Aee=200.0, Aei=50.0, Aii=200, Aie=50.0),
                          (Aee=150.0, Aei=150.0, Aii=50, Aie=115.0),
                          (Aee=50.0, Aei=200.0, Aii=200, Aie=50.0)
                         ],  fcmb_monotonic_A_fpath, fcmb_blocking_A_fpath, :has_propagation, "contrast_monotonic_blocking_A")

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
