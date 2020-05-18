function cb_E_is_fully_propagated_with_save_idxs(u,t,integrator)
    sub_u = population(u[integrator.opts.save_idxs],1);
    t > 5 && ((all(isapprox.(sub_u, 0.0, atol=1e-4)) || (sub_u[end-10] > 0.005)))
end
function cb_E_is_fully_propagated(u,t,integrator)
    pop = population(u,1)
    t > 5 && ((all(isapprox.(u, 0.0, atol=1e-4)) || (pop[end-10] > 0.005)))
end