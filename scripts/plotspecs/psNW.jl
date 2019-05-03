plotspecs = [
  NonlinearityPlot(;
    fn_bounds = (-1,15)
    ),
  # SpaceTimePlot(),
  SubsampledPlot(
    plot_type=WaveStatsPlot,
    time_subsampler=Subsampler(
      Δ = 0.01,
      window = (1.2, 1.8)
    ),
    space_subsampler=Subsampler(
        window = (5.0,Inf)
      )
    ),
    peak_interpolation_n = 2
]
