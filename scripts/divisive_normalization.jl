using DrWatson
quickactivate(@__DIR__, "TravelingWaveSimulations")

using Simulation73, TravelingWaveSimulations, NeuralModels
using Lazy
ENV["GKSwstype"] = "100" # For headless plotting (on server)
ENV["MPLBACKEND"]="Agg"
using Plots
pyplot()

strength = [1.2, 1.2]
width = [28.1, 28.1]
SNR = [80.0, 80.0]
time_windows = [[(0.0, 55.0)], [(0.0, 55.0)]]
make_simulation = get_example("neuman_line")
stimulus = MultipleDifferentStimuli([
    pops(GaussianNoiseStimulus; SNR = SNR),
    pops(MultipleSameStimuli{SharpBumpStimulus};
        strength=strength,
        width=width,
        time_windows=time_windows,
        centers=centers
    )
])

neuman_line = get_example("neuman_line_two_stimuli")
