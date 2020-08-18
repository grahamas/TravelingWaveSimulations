# First activate the TravelingWaveSimulations project
using Pkg; Pkg.activate(@__DIR__)

# Second, load the TravelingWaveSimulations package
using TravelingWaveSimulations

# Finally, run a single simulation based on an example (located in src/examples)
based_on_example(; example_name="sigmoid_normal_fft",
                   modifications=["iiS=1.0", "x=1000.0"],
                   analyses=["animate"])

# You can specify any analysis in scripts/analyses, and any modification either 1) in scripts/modifications (ignore this option) OR any (any!) keyword argument you see in the example specification.  You can run multiple simulations at once by specifying e.g. "iiS=1.0:1.0:3.0"
