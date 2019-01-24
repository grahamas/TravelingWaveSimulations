module WCM

using Parameters
using CalculatedParameters
import CalculatedParameters: Calculated, update
using Simulating
using Modeling
using WCMConnectivity
using WCMNonlinearity
using WCMStimulus
using Meshes
#using Exploration
using DifferentialEquations
using Targets
using Records
import Records: required_modules
using StaticArrays
using Exploration
using Subsampling

# Rename to remove N redundancy
struct WCMSpatial1D{T,N,P,C<:Connectivity{T},
                            L<:Nonlinearity{T},S<:Stimulus{T},SP<:PopSpace{T,N,P}} <: Model{T,N,P}
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
        space::SP, connectivity::Array{C,2}, nonlinearity::Array{L,1}, stimulus::Array{S,1}) where {T,P,N,Str<:AbstractString,C<:Connectivity{T},L<:Nonlinearity{T},S<:Stimulus{T},SP<:Space{T,N}}
    WCMSpatial1D{T,N,P,C,L,S,SP}(SVector{P,T}(α),SVector{P,T}(β),SVector{P,T}(τ),space,SMatrix{P,P,C}(connectivity),SVector{P,L}(nonlinearity),SVector{P,S}(stimulus),SVector{P,Str}(pop_names))
end

space_array(model::WCMSpatial1D) = Calculated(model.space).value

# import Exploration: base_type
# function base_type(::Type{WCMSpatial1D{T1,T2,T3,T4,T5,T6,T7}}) where {T1,T2,T3,T4,T5,T6,T7}
#     BT1 = base_type(T1); BT2 = base_type(T2); BT3 = base_type(T3)
#     BT4 = base_type(T4); BT5 = base_type(T5); BT6 = base_type(T6)
#     BT7 = base_type(T7)#; BT8 = base_type(T8)
#     return WCMSpatial1D{BT1,BT2,BT3,BT4,BT5,BT6,BT7}
# end

function subsampling_time_idxs(t_target, solver)
    t_solver = time_span(solver)[1]:save_dt(solver):time_span(solver)[end]
    subsampling_idxs(t_target, t_solver)
end
function subsampling_space_idxs(x_target, model, solver)
    x_model = space_arr(model)[1:solver.space_save_every:end,1]
    subsampling_idxs(x_target, x_model)
end

export subsampling_time_idxs, subsampling_space_idxs



# * Calculated WC73 Simulation Type
mutable struct CalculatedWCMSpatial1D{T,N,P,C,L,S,CC<:CalculatedParam{C},CL <: CalculatedParam{L},CS <: CalculatedParam{S}} <: CalculatedParam{WCMSpatial1D{T,N,P,C,L,S}}
    α::SVector{P,T}
    β::SVector{P,T}
    τ::SVector{P,T}
    connectivity::SMatrix{P,P,CC}
    nonlinearity::SVector{P,CL}
    stimulus::SVector{P,CS}
end

function CalculatedWCMSpatial1D(wc::WCMSpatial1D{T,N,P,C,L,S}) where {T<:Real,N,P,
                                                  C<:Connectivity{T},
                                                  L<:Nonlinearity{T},
                                                  S<:Stimulus{T}}
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

function get_values(cwc::CalculatedWCMSpatial1D{T,N}) where {T,N}
    (cwc.α, cwc.β, cwc.τ, get_value.(cwc.connectivity), get_value.(cwc.nonlinearity), get_value.(cwc.stimulus))
end



import Exploration: make_problem_generator

function make_calculated_function(cwc::CalculatedWCMSpatial1D{T,1,P,C,L,S,CC,CL,CS}) where {T,P,C<:Connectivity{T},L<:Nonlinearity{T},S<:Stimulus{T},CC<:CalculatedParam{C},CL <: CalculatedParam{L},CS <: CalculatedParam{S}}
    (α, β, τ, connectivity_mx, nonlinearity_objs, stimulus_objs) = get_values(cwc)

    let α::SVector{P,T}=α, β::SVector{P,T}=β, τ::SVector{P,T}=τ, connectivity_mx::SMatrix{P,P,Matrix{T}}=connectivity_mx, nonlinearity_objs::SVector{P,CL}=nonlinearity_objs, stimulus_objs::SVector{P,CS}=stimulus_objs
        (dA::Array{T,2}, A::Array{T,2}, p::Union{Array{T,1},Nothing}, t::T) -> (
            @views for i in 1:P
                stimulate!(dA[:,i], stimulus_objs[i], t) # I'll bet it goes faster if we pull this out of the loop
                for j in 1:P
                    dA[:,i] .+= connectivity_mx[i,j] * A[:,j]
                end
                # dA[:,i] .+= sum(connectivity_mx[i,j] * A[:,j] for j in 1:2)
                nonlinearity!(dA[:,i], nonlinearity_objs[i])
                dA[:,i] .*= β[i] .* (1.0 .- A[:,i])
                dA[:,i] .+= -α[i] .* A[:,i]
                dA[:,i] ./= τ[i]
            end
        )
    end
end

function make_problem_generator(model::M, solver::SV, variable_map) where {T,M<:WCMSpatial1D{T},SV<:Solver{T}}
    tspan = time_span(solver)
    u0 = initial_value(model)

    cwc = Calculated(model)

    function problem_generator(prob, new_p::Array{T})
        update_from_p!(cwc, new_p, model, variable_map)
        WilsonCowan73! = make_calculated_function(cwc)
        ode_fn = convert(ODEFunction{true}, WilsonCowan73!)
        return ODEProblem(ode_fn, u0, tspan, new_p)
    end

    return problem_generator
end

function Simulating.generate_problem(model::M, solver::SV) where {T,M<:WCMSpatial1D{T},SV<:Solver{T}}
    tspan = time_span(solver)
    u0 = initial_value(model)

    cwc = Calculated(model)

    WilsonCowan73! = make_calculated_function(cwc)

    ode_fn = convert(ODEFunction{true}, WilsonCowan73!)
    return ODEProblem(ode_fn, u0, tspan, nothing)
end

export WCMSpatial1D, space_array
export base_type
export generate_problem

end