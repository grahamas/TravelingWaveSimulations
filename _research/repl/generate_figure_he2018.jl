include("figure_contrast_monotonic_blocking_propagation.jl")

he2018_monotonic_fpath = get_recent_simulation_data_path(joinpath(homedir(), "data", "harris_ermentrout_monotonic"))
he2018_blocking_fpath = get_recent_simulation_data_path(joinpath(homedir(), "data", "harris_ermentrout_blocking"))

save_figure_contrast_monotonic_blocking((:tau, :sigma), he2018_monotonic_fpath, he2018_blocking_fpath, :has_propagation, "he2018", slice_y_intersect=0.0, slice_slope=1.0)

