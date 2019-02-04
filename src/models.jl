# Rename to remove N redundancy
struct WCMSpatial1D{T,N,P,C<:AbstractConnectivity{T},
                            L<:AbstractNonlinearity{T},
                            S<:AbstractStimulus{T},
                            SP<:PopSpace{T,N,P}} <: Model{T,N,P}
    α::SVector{P,T}
    β::SVector{P,T}
    τ::SVector{P,T}
    space::SP
    connectivity::SMatrix{P,P,C}
    nonlinearity::SVector{P,L}
    stimulus::SVector{P,S}
    pop_names::SVector{P,String}
end

function WCMSpatial1D{T,N,P}(; pop_names::Array{Str,1}, α::Array{T,1}, β::Array{T,1}, τ::Array{T,1},
        space::SP, connectivity::Array{C,2}, nonlinearity::Array{L,1}, stimulus::Array{S,1}) where {T,P,N,Str<:AbstractString,C<:AbstractConnectivity{T},L<:AbstractNonlinearity{T},S<:AbstractStimulus{T},SP<:Space{T,N}}
    WCMSpatial1D{T,N,P,C,L,S,SP}(SVector{P,T}(α),SVector{P,T}(β),SVector{P,T}(τ),space,SMatrix{P,P,C}(connectivity),SVector{P,L}(nonlinearity),SVector{P,S}(stimulus),SVector{P,Str}(pop_names))
end

space_array(model::WCMSpatial1D) = Calculated(model.space).value


# * Calculated WC73 Simulation Type
mutable struct CalculatedWCMSpatial1D{T,N,P,C,L,S,CC<:CalculatedType{<:C},CL <: CalculatedType{<:L},CS <: CalculatedType} <: CalculatedType{WCMSpatial1D{T,N,P,C,L,S}}
    α::SVector{P,T}
    β::SVector{P,T}
    τ::SVector{P,T}
    connectivity::SMatrix{P,P,CC}
    nonlinearity::SVector{P,CL}
    stimulus::SVector{P,CS}
end

function CalculatedWCMSpatial1D(wc::WCMSpatial1D{T,N,P,C,L,S}) where {T<:Real,N,P,
                                                  C<:AbstractConnectivity{T},
                                                  L<:AbstractNonlinearity{T},
                                                  S<:AbstractStimulus{T}}
    connectivity = Calculated.(wc.connectivity, Ref(wc.space))
    nonlinearity = Calculated.(wc.nonlinearity)
    stimulus = Calculated.(wc.stimulus,Ref(wc.space))
    CC = eltype(connectivity)
    CL = eltype(nonlinearity)
    CS = eltype(stimulus)
    CalculatedWCMSpatial1D{T,N,P,C,L,S,CC,CL,CS}(
        wc.α, wc.β, wc.τ,
        connectivity, nonlinearity, stimulus)
end

function Calculated(wc::WCMSpatial1D)
    CalculatedWCMSpatial1D(wc)
end

function get_values(cwc::CalculatedWCMSpatial1D{T,N}) where {T,N}
    (cwc.α, cwc.β, cwc.τ, get_value.(cwc.connectivity), get_value.(cwc.nonlinearity), get_value.(cwc.stimulus))
end

function update_from_p!(cwc::CalculatedWCMSpatial1D{T}, new_p::Array{T}, model, variable_map) where T
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

function make_calculated_function(cwc::CalculatedWCMSpatial1D{T,1,P,C,L,S,CC,CL,CS}) where {T,P,C<:AbstractConnectivity{T},L<:AbstractNonlinearity{T},S<:AbstractStimulus{T},CC<:CalculatedType{<:C},CL <: CalculatedType{<:L},CS <: CalculatedType}
    (α, β, τ, connectivity_mx, nonlinearity_objs, stimulus_objs) = get_values(cwc)

    let α::SVector{P,T}=α, β::SVector{P,T}=β, τ::SVector{P,T}=τ, connectivity_mx::SMatrix{P,P,Matrix{T}}=connectivity_mx, nonlinearity_objs::SVector{P,CL}=nonlinearity_objs, stimulus_objs::SVector{P,CS}=stimulus_objs
        (dA::Array{T,2}, A::Array{T,2}, p::Union{Array{T,1},Nothing}, t::T) -> (
            @views for i in 1:P
                @show i
                stimulate!(dA[:,i], stimulus_objs[i], t) # I'll bet it goes faster if we pull this out of the loop
                for j in 1:P
                    dA[:,i] .+= connectivity_mx[i,j] * A[:,j]
                end
                # dA[:,i] .+= sum(connectivity_mx[i,j] * A[:,j] for j in 1:2)
                nonlinearity!(dA[:,i], nonlinearity_objs[i])
                dA[:,i] .*= β[i] .* (1.0 .- A[:,i])
                dA[:,i] .+= -α[i] .* A[:,i]
                dA[:,i] ./= τ[i]
                @show sum(dA[:,i])
                @show maximum(dA[:,i])                
            end
        )
    end
end
