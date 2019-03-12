# Rename to remove N redundancy
struct WCMSpatial{T,D,P,C<:AbstractConnectivity{T,D},
                            L<:AbstractNonlinearity{T},
                            S<:AbstractStimulus{T,D},
                            SP<:Pops{P,T,D}} <: Model{T,D,P}
    α::SVector{P,T}
    β::SVector{P,T}
    τ::SVector{P,T}
    space::SP
    connectivity::SMatrix{P,P,C}
    nonlinearity::SVector{P,L}
    stimulus::SVector{P,S}
    pop_names::SVector{P,String}
end

struct DualModel{T,D,P,M<:Model{T,D,P}} <: Model{T,D,P}
    model::M
end
function Base.getproperty(dm::DualModel, sym::Symbol)
    if sym == :model
        return getfield(dm, sym)
    else
        return getproperty(getfield(dm, :model), sym)
    end
end
DualModel{T,D,P,M}(; kwargs...) where {T,D,P,M<:Model{T,D,P}} = DualModel{T,D,P,M}(M(;kwargs...))
#space_array(dm::DualModel) = space_array(model(dm))
struct CalculatedDualModel{CM<:CalculatedType{<:Model}}
    calculated_model::CM
end
function Calculated(dm::DualModel)
    CalculatedDualModel(Calculated(dm.model))
end
get_values(cdm::CalculatedDualModel) = get_values(cdm.calculated_model)

function WCMSpatial{T,N,P}(; pop_names::Array{Str,1}, α::Array{T,1}, β::Array{T,1}, τ::Array{T,1},
        space::SP, connectivity::Array{C,2}, nonlinearity::Array{L,1}, stimulus::Array{S,1}) where {T,P,N,Str<:AbstractString,C<:AbstractConnectivity{T},L<:AbstractNonlinearity{T},S<:AbstractStimulus{T},SP<:AbstractSpace{T,N}}
    WCMSpatial{T,N,P,C,L,S,SP}(SVector{P,T}(α),SVector{P,T}(β),SVector{P,T}(τ),space,SMatrix{P,P,C}(connectivity),SVector{P,L}(nonlinearity),SVector{P,S}(stimulus),SVector{P,Str}(pop_names))
end

space_array(model::WCMSpatial) = Calculated(model.space).value


# * Calculated WC73 Simulation Type
mutable struct CalculatedWCMSpatial{T,D,P,C,L,S,CC<:CalculatedType{<:C},CL <: CalculatedType{<:L},CS <: CalculatedType} <: CalculatedType{WCMSpatial{T,D,P,C,L,S}}
    α::SVector{P,T}
    β::SVector{P,T}
    τ::SVector{P,T}
    connectivity::SMatrix{P,P,CC}
    nonlinearity::SVector{P,CL}
    stimulus::SVector{P,CS}
end

function CalculatedWCMSpatial(wc::WCMSpatial{T,N,P,C,L,S}) where {T<:Real,N,P,
                                                  C<:AbstractConnectivity{T},
                                                  L<:AbstractNonlinearity{T},
                                                  S<:AbstractStimulus{T,N}}
    calc_space = Calculated(wc.space)
    connectivity = Calculated.(wc.connectivity, Ref(calc_space))
    nonlinearity = Calculated.(wc.nonlinearity)
    stimulus = Calculated.(wc.stimulus,Ref(calc_space))
    CC = eltype(connectivity)
    CL = eltype(nonlinearity)
    CS = eltype(stimulus)
    CalculatedWCMSpatial{T,N,P,C,L,S,CC,CL,CS}(
        wc.α, wc.β, wc.τ,
        connectivity, nonlinearity, stimulus)
end

function Calculated(wc::WCMSpatial)
    CalculatedWCMSpatial(wc)
end

function get_values(cwc::CalculatedWCMSpatial{T,N}) where {T,N}
    (cwc.α, cwc.β, cwc.τ, get_value.(cwc.connectivity), get_value.(cwc.nonlinearity), get_value.(cwc.stimulus))
end

function update_from_p!(cwc::CalculatedWCMSpatial{T}, new_p::Array{T}, model, variable_map) where T
    # Use the variable model stored by p_search to create static model
    new_model = model_from_p(model, variable_map, new_p)

    # Update the calculated values from the new static model
    cwc.α = new_model.α
    cwc.β = new_model.β
    cwc.τ = new_model.τ
    cwc.connectivity = update(cwc.connectivity, new_model.connectivity, new_model.space)
    cwc.nonlinearity = update(cwc.nonlinearity, new_model.nonlinearity)
    cwc.stimulus = update(cwc.stimulus, new_model.stimulus, new_model.space)
end


macro make_make_function(num_dims)
    D = eval(num_dims)
    to_syms = [Symbol(:to,:_,i) for i in 1:D]
    from_syms = [Symbol(:from,:_,i) for i in 1:D]
    space_colons = [:(:) for i in 1:D]
    D_P = D + 1
    D_CONN_P = D + D + 2
    D_CONN = D + D
    tensor_prod_expr = @eval @macroexpand @tensor dA[$(to_syms...),i] = dA[$(to_syms...),i] + connectivity_tensor[$(to_syms...),$(from_syms...),i,j] * A[$(from_syms...),j]
    #tensor_prod_expr = @eval @macroexpand @tensor dA[$(to_syms...),i] = dA[$(to_syms...),i] + connectivity_sarray[i,j][$(to_syms...),$(from_syms...)] * A[$(from_syms...),j]
    #tensor_prod_expr = @eval @macroexpand @tensor dA[to_1, i] = dA[to_1, i] + connectivity_tensor[to_1, from_1, i, j] * A[from_1, j]
    esc(quote
        function make_calculated_function(cwc::CalculatedWCMSpatial{T,$D,P,C,L,S,CC,CL,CS}) where {T,P,C<:AbstractConnectivity{T},L<:AbstractNonlinearity{T},S<:AbstractStimulus{T},CC<:CalculatedType{<:C},CL <: CalculatedType{<:L},CS <: CalculatedType}
            (α, β, τ, connectivity_sarray, nonlinearity_objs, stimulus_objs) = get_values(cwc)
            connectivity_tensor = Array{T,$D_CONN_P}(undef, size(connectivity_sarray[1])..., P, P)
            indices = CartesianIndex.(Iterators.product(1:P, 1:P))
            for index in indices
                connectivity_tensor[$(space_colons...), $(space_colons...), index] .= connectivity_sarray[index]#, size(connectivity_sarray[1])..., 1, 1)
            end
            let α::SVector{P,T}=α, β::SVector{P,T}=β, τ::SVector{P,T}=τ, connectivity_tensor::Array{T,$D_CONN_P}=connectivity_tensor, nonlinearity_objs::SVector{P,CL}=nonlinearity_objs, stimulus_objs::SVector{P,CS}=stimulus_objs
                (dA::Array{T,$D_P}, A::Array{T,$D_P}, p::Union{Array{T,1},Nothing}, t::T) -> begin
                    @views for i in 1:P
                        stimulate!(dA[$(space_colons...),i], stimulus_objs[i], t) # I'll bet it goes faster if we pull this out of the loop
                    end
                    @views $(tensor_prod_expr)
                    # for j in 1:P
                    #     dA[:,i] .+= connectivity_tensor[:,:,i,j] * A[:,j]
                    # end
                    @views for i in 1:P
                        #@show sum(dA[$(space_colons...),i])
                        nonlinearity!(dA[$(space_colons...),i], nonlinearity_objs[i])
                        #@show sum(dA[$(space_colons...),i])
                        dA[$(space_colons...),i] .*= β[i] .* (1.0 .- A[$(space_colons...),i])
                        dA[$(space_colons...),i] .+= -α[i] .* A[$(space_colons...),i]
                        dA[$(space_colons...),i] ./= τ[i]
                    end
                    dA
                end
            end
        end
    end)
end

#@show @macroexpand @make_make_function 1
@make_make_function 1
@make_make_function 2
