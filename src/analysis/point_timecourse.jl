# Analyse as if from a fixed electrode

function point_timecourse(exec::Execution, location_proportion)
    us = exec.solution.u
    t = timepoints(exec)
    xs = space(exec)
    n_xs = size(xs)
    x_dx = floor.(Int, n_xs .* location_proportion)
    @show x_dx
    x = xs[x_dx...]
    u = [[u[x_dx...]] for u in us]
    return BareSolution(u=u,x=[x],t=t)    
end

function point_average(u::AbstractArray, pt::Tuple, xs, σ)
    pt_dist = map(xs) do x
        -sum((x .- pt) ./ σ) .^ 2) / 2
    end
    unscaled = exp.(pt_dist)
    return unscaled ./ sum(unscaled)
end
function point_average_timecourse(exec::Execution, location_proportion, σ)
    us = exec.solution.u
    t = timepoints(exec)
    xs = space(exec)
    n_xs = size(xs)
    x_dx = floor.(Int, n_xs .* location_proportion)
    @show x_dx
    x = xs[x_dx...]
    u = [[point_average(u, x, xs, σ)] for u in us]
    return BareSolution(u=u,x=[x],t=t)  
end