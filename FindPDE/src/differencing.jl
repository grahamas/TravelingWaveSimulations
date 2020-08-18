
function finite_diff(u::AbstractArray{T}, dx, derivative_order, approximation_order=2) where T
    dim = 1
    # TODO: Test against FiniteDiff in PDE-FIND
    derivative_operator = CenteredDifference{dim}(derivative_order, approximation_order, dx, size(u, dim))
    boundary_operator = Neumann0BC(dx, 1) #FIXME: sends boundary to zero
    op = derivative_operator * boundary_operator
    return op * u
end

function manual_finite_diff(u::AbstractArray{T,1}, dx, d, approx_order=2) where T
    """
    Takes dth derivative data using 2nd order finite difference method (up to d=3)
    Works but with poor accuracy for d > 3
    
    Input:
    u = data to be differentiated
    dx = Grid spacing.  Assumes uniform spacing
    """
    @assert approx_order == 2
    
    n = size(u,1)
    ux = Array{T}(undef, n)
    
    if d == 1
        for i in 2:n-1
            ux[i] = (u[i+1]-u[i-1]) / (2*dx)
        end
        
        ux[1] = (-(3.0/2)*u[1] + 2*u[2] - u[3]/2) / dx
        ux[n] = ((3.0/2)*u[n] - 2*u[n-1] + u[n-2]/2) / dx
        return ux
    end
    
    if d == 2
        for i in 2:n-1
            ux[i] = (u[i+1]-2*u[i]+u[i-1]) / dx^2
        end
        
        ux[1] = (2*u[1] - 5*u[2] + 4*u[3] - u[4]) / dx^2
        ux[n] = (2*u[n] - 5*u[n-1] + 4*u[n-2] - u[n-3]) / dx^2
        return ux
    end
    
    if d == 3
        for i in 3:n-2
            ux[i] = (u[i+2]/2-u[i+1]+u[i-1]-u[i-2]/2) / dx^3
        end
        
        ux[1] = (-2.5*u[1]+9*u[2]-12*u[3]+7*u[4]-1.5*u[5]) / dx^3
        ux[2] = (-2.5*u[2]+9*u[3]-12*u[4]+7*u[5]-1.5*u[6]) / dx^3
        ux[n] = (2.5*u[n]-9*u[n-1]+12*u[n-2]-7*u[n-3]+1.5*u[n-4]) / dx^3
        ux[n-1] = (2.5*u[n-1]-9*u[n-2]+12*u[n-3]-7*u[n-4]+1.5*u[n-5]) / dx^3
        return ux
    end
    
    if d > 3
        error("manual_finite_diff unimplemented for d > 3")
    end
end

function manual_finite_diff(u::AbstractArray{T,2}, dx, d, approx_order=2) where T
    ux = Array{T}(undef, size(u)...)
    for i_col in 1:size(u,2)
        ux[:,i_col] = manual_finite_diff(u[:,i_col], dx, d, approx_order)
    end
    return ux
end
