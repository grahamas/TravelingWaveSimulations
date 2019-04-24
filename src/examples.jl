module Examples

using WilsonCowanModel
using Simulation73
import DifferentialEquations: Euler

function neuman_line(;
        α = [1.1, 1.0],
        β = [1.1, 1.1],
        τ = [0.1, 0.18],
        a = [1.2, 1.0],
        θ = [2.6, 8.0],
        n_points=301,
        extent=100.0,
        strength = [1.2, 1.2],
        width = [2.81, 2.81],
        SNR = [80.0, 80.0],
        time_window = [(0.0, 0.55), (0.0, 0.55)],
        amplitude = [16.0 -18.2;
                     27.0 -4.0],
        spread = [2.5 2.7;
                  2.7 2.5],
        stop_time = 1.8,
        dt = 0.01,
        space_save_every=1,
        time_save_every=1,
        algorithm=Euler()
    )
    N=1; P=2;
    simulation = Simulation(;
      model = WCMSpatial{Float64,N,P}(;
        pop_names = ["E", "I"],
        α = α,
        β = β,
        τ = τ,
        space = Pops{P}(Segment{Float64}(; n_points=n_points, extent=extent)),
        nonlinearity = pops(SigmoidNonlinearity{Float64};
          a = a,
          θ = θ),
          # a=Float64[BV(1.2, (0.1, 2.0)), BV(1.0, (0.1, 2.0))],
          # θ=Float64[BV(8.0, (2.0, 9.0)), BV(2.6, (2.0,9.0))]),
        stimulus = pops(NoisyStimulus{Float64,N};
          strength=strength,
          time_window=time_window,
          width=width,
          SNR=SNR,
          stim_type=[SharpBumpStimulus{Float64,N}, SharpBumpStimulus{Float64,N}]),
        connectivity = pops(ShollConnectivity{Float64};
          amplitude = amplitude,
          # spread = Float64[BV(2.5, (1.0, 4.0)) BV(2.7, (1.0, 4.0));
          #            BV(2.7, (1.0, 4.0)) BV(2.5, (1.0, 4.0))])
          spread = spread)
        ),
      solver = Solver{Float64}(;
        stop_time = stop_time,
        dt = dt,
        space_save_every=space_save_every,
        time_save_every=time_save_every,
        #stiffness=:stiff
        algorithm=algorithm
      )
    )
    return simulation
end

function neuman_square(;
    α = [1.1, 1.0],
    β = [1.1, 1.1],
    τ = [0.1, 0.18],
    a = [1.2, 1.0],
    θ = [2.6, 8.0],
    n_points=(51,51),
    extent=(50.0,50.0),
    strength = [1.2, 1.2],
    width = [2.81, 2.81],
    SNR = [80.0, 80.0],
    time_window = [(0.0, 0.55), (0.0, 0.55)],
    amplitude = [16.0 -18.2;
                 27.0 -4.0],
    spread = [(2.5,2.5) (2.7,2.7);
              (2.7,2.7) (2.5,2.5)],
    stop_time = 1.8,
    dt = 0.01,
    space_save_every=1,
    time_save_every=1,
  algorithm=Euler()
  )
  N=2
  P=2
  simulation = Simulation(;
    model = WCMSpatial{Float64,N,P}(;
      pop_names = ["E", "I"],
      α = α,
      β = β,
      τ = τ,
      space = Pops{P}(Grid{Float64}(; n_points=n_points, extent=extent)),
      nonlinearity = pops(SigmoidNonlinearity{Float64};
        a = a,
        θ = θ),
        # a=Float64[BV(1.2, (0.1, 2.0)), BV(1.0, (0.1, 2.0))],
        # θ=Float64[BV(8.0, (2.0, 9.0)), BV(2.6, (2.0,9.0))]),
      stimulus = pops(NoisyStimulus{Float64,N};
        strength=strength,
        time_window=time_window,
        width=width,
        SNR=SNR,
        stim_type=[SharpBumpStimulus{Float64,N}, SharpBumpStimulus{Float64,N}]),
      connectivity = pops(GaussianConnectivity{Float64};
        amplitude = amplitude,
        # spread = Float64[BV(2.5, (1.0, 4.0)) BV(2.7, (1.0, 4.0));
        #            BV(2.7, (1.0, 4.0)) BV(2.5, (1.0, 4.0))])
        spread = spread)
      ),
    solver = Solver{Float64}(;
      stop_time = stop_time,
      dt = dt,
      space_save_every=space_save_every,
      time_save_every=time_save_every,
      #stiffness=:stiff
      algorithm=algorithm
      )
  )
end


end
