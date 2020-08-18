

function multinomial_powers(n_variables, max_order)
    # NOTE: n_variables should NOT include constant
    Iterators.flatten(map(1:(max_order+1)) do order
        overflow_size = n_variables + order-1
        map(combinations(0:overflow_size-1, n_variables-1)) do indices
            starts = [0, (indices .+ 1)...]
            stops = [indices..., overflow_size]
            stops .- starts
        end
    end)
end
function multinomial(bases::AbstractVector, powers::AbstractVector)
    prod(bases .^ powers)
end


function multinomial_recombination(data::AbstractMatrix{T}, data_description, max_exponent::Int) where T
    """
    You might be tempted to apply a list of transforms arbitrarily, but 
    then you might end up with equivalent transforms, so that the columns
    are not independent.

    Build a matrix with columns representing multinomials up to degree max_exponent
    of all columns in data.

    input:
        data: column 1 is U; remaining columns are transforms of U including
              derivatives
        data_description: labels for each column
        max_exponent: maximum polynomial order of resulting columns
    output: 
        multinomials: multinomial combinations of columns in input data
        multinomials_description: descriptions of columns in multinomials
    """
    len_vars, n_vars = size(data)
    l_powers = multinomial_powers(n_vars, max_exponent) |> collect # don't power constant column
    n_terms = length(l_powers)
    multinomial_data = Matrix{T}(undef, len_vars, n_terms)
    multinomial_data_description = Vector{AbstractDescription}(undef, n_terms)
    for (i_power, powers) in enumerate(l_powers)
        for i_row in 1:len_vars
            multinomial_data[i_row,i_power] = multinomial(data[i_row,:], powers')
        end
        multinomial_data_description[i_power] = multinomial(data_description, powers)
    end

    return multinomial_data, multinomial_data_description
end

# only raises u^p, not u_x^p
# function impoverished_multinomial_recombination(us::AbstractMatrix{T}, us_desc, dus::AbstractMatrix{T}, dus_desc,
#         max_exponent::Int) where T
#     """
#     You might be tempted to apply a list of transforms arbitrarily, but 
#     then you might end up with equivalent transforms, so that the columns
#     are not independent.

#     Build a matrix with columns representing multinomials up to degree max_exponent
#     of all columns in data.

#     input:
#         data: column 1 is C; column 2 is U; remaining columns are transforms of U including
#               derivatives
#         data_description: labels for each column
#         max_exponent: maximum polynomial order of resulting columns
#     output: 
#         multinomials: multinomial combinations of columns in input data
#         multinomials_description: descriptions of columns in multinomials
#     """
    
#     len_vars, n_derivatives = size(dus)
#     _, n_us = size(us)
#     n_terms = (n_derivatives + n_us) + n_derivatives * max_exponent + (max_exponent - 1)
#     error("unimplemented")
#     ### UNIMPLEMENTED; BROKEN
#     multinomial_data = Matrix{T}(undef, len_vars, n_terms)
#     multinomial_data[:,1:n_vars] .= data
#     multinomial_data_description = Vector{AbstractDescription}(undef, n_terms)
#     multinomial_data_description[1:n_vars] .= data_description
#     count = n_vars + 1
#     for exponent in 1:max_exponent
#         u_p = data[:,base_idx] .^ exponent
#         u_p_desc = data_description[base_idx] ^ exponent
#         if exponent != 1
#             multinomial_data[:,count] .= u_p
#             multinomial_data_description[count] = u_p_desc
#             count += 1
#         end
#         for i_derivative in first_derivative_idx:n_vars
#             multinomial_data[:,count] = u_p .* data[:,i_derivative]
#             multinomial_data_description[count] = u_p_desc * data_description[i_derivative]
#             count += 1
#         end
#     end

#     return multinomial_data, multinomial_data_description
# end


