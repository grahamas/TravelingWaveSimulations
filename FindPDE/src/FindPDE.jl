module FindPDE

using DiffEqOperators
using DataFrames
#using MLJ
using LinearAlgebra
using Combinatorics
using LaTeXStrings
using PrettyPrinting

export build_linear_system, train_STRidge

include("differencing.jl")
include("multinomials.jl")
include("descriptions.jl")
include("regression.jl")

function build_linear_system(matrix_u::Matrix{T}, dt, dx, derivative_order, polynomial_order, deriv_approx_order=2) where
    {T}
    # TODO: multidimensional space
    n_x, n_t = size(matrix_u)

    # Matrix where each column is a flattened derivative matrix
    # Col 1: u
    # Cols 2-(derivative_order+1): u_(x...)
    n_space_derivatives = derivative_order + 2
    space_derivatives = Matrix{T}(undef, prod(size(matrix_u)), n_space_derivatives) 
    space_derivatives_desc = Vector{AbstractDescription}(undef, n_space_derivatives)
    
    # Add constant column
    space_derivatives[:,1] .= 1.0
    constant_desc = AtomicDescription(:C)
    space_derivatives_desc[1] = constant_desc

    # Add u (derivative order 0) column
    space_derivatives[:,2] .= matrix_u[:]
    base_desc = AtomicDescription(:u)
    space_derivatives_desc[2] = base_desc

    time_derivative_order = 1
    time_derivative_approximation_order = deriv_approx_order
    u_t = (finite_diff(matrix_u', dt, time_derivative_order, 
                      time_derivative_approximation_order)')[:]
    u_t_desc = differentiate(base_desc, :t)

    for order in 1:derivative_order
        space_derivatives[:,order+2] .= finite_diff(matrix_u, dx, order, deriv_approx_order)[:]
        space_derivatives_desc[order+2] = differentiate(base_desc, :x, order)
    end
    @show sum(isnan.(space_derivatives))
   
    # FIXME don't send constant
    space_terms, space_terms_desc =  impoverished_multinomial_recombination(space_derivatives, space_derivatives_desc, polynomial_order)
#     all_terms = hcat(u_t, space_terms)
#     all_terms_desc = [u_t_desc, space_terms_desc...]
#     all_terms_df = convert(DataFrame, all_terms)
#     rename!(all_terms_df, Symbol.(all_terms_desc))
#     all_terms_df[!,:C] .= 1.0
    return u_t, space_terms, space_terms_desc
end

function build_linear_system(u::Array{T,3}, u_names::Vector{Symbol}, dt, dx, derivative_order, polynomial_order, space_deriv_approx_order=2, time_deriv_approx_order=2) where
    {T}
    # TODO: multidimensional space; Currently must be 1D
    # n_x: num points of 1d space
    # n_t: num points of time
    # n_p: number of "populations"
    n_x, n_t, n_p = size(u)
    @assert length(u_names) == size(u,3)
    u_descs = AtomicDescription.(u_names)

    # Matrix where each column is a flattened derivative matrix
    # Col 1: Constant
    # Col 2-P+1: u[:,:,i]
    # Cols P+1-(derivative_order+1): u1_(x...)
    # cols (derivative_order+1)-...: u2_(y...)
    n_space_derivatives = 1 + (derivative_order + 1) * n_p
    space_derivatives = Matrix{T}(undef, prod(size(u)[[1,2]]), n_space_derivatives) 
    space_derivatives_desc = Vector{AbstractDescription}(undef, n_space_derivatives)
    time_derivatives = Matrix{T}(undef, prod(size(u)[[1,2]]), n_p)
    time_derivatives_desc = Vector{AbstractDescription}(undef, n_p)
    
    # Add constant column
    space_derivatives[:,1] .= 1.0
    constant_desc = AtomicDescription(:C)
    space_derivatives_desc[1] = constant_desc

    # Add u (derivative order 0) column
    time_derivative_order = 1
    n_nonderivatives = u_p + 1
    for i_p in 1:n_p
        this_matrix = u[:,:,i_p]
        space_derivatives[:,1+i_p] .= this_matrix[:]
        space_derivatives_desc[2] = u_descs[i_p]
        n_derivatives_so_far = derivative_order * (i_p - 1)
        this_p_zero_idx = n_nonderivatives + n_derivatives_so_far
        for order in 1:derivative_order
            space_derivatives[:,order+this_p_zero_idx] .= finite_diff(this_matrix, dx, order, space_deriv_approx_order)[:]
            space_derivatives_desc[order+this_p_zero_idx] = differentiate(u_descs[i_p], :x, order)
        end
        
        time_derivatives[:, i_p] = (finite_diff(this_matrix', dt, time_derivative_order, 
                      time_deriv_approx_order)')[:]
        time_derivatives_desc[i_p] = differentiate(u_descs[i_p], :t)
    end

    # FIXME don't send constant
    space_terms, space_terms_desc =  multinomial_recombination(
        space_derivatives[:,2:end], space_derivatives_desc[2:end], 
        polynomial_order)
    space_terms = hcat(space_derivatives[:,1], space_terms)
    space_terms_desc = hcat(space_derivatives_desc[1], space_terms_desc)
    return time_derivatives, time_derivatives_desc, space_terms, space_terms_desc
end



end # module
