# NOTE: P is currently the TRAILING dimension. I haven't revised comments yet.

# Rename to remove N redundancy
struct WCMSpatial{T,N_CDT,P,
        SCALARS<:NTuple{P,T},
        CONN<:PopInteractParam{P},#<:MaybeNoOpComposite{T,<:AbstractConnectivityParameter{T,N_CDT}}},
        NONL<:PopAct{P,<:AbstractNonlinearity{T}},
        STIM<:PopActParam{P,<:AbstractStimulusParameter{T}}
    } <: AbstractODEModel{T,N_CDT,P}
    α::SCALARS
    β::SCALARS
    τ::SCALARS
    connectivity::CONN
    nonlinearity::NONL
    stimulus::STIM
    pop_names::NTuple{P,String}
    function WCMSpatial{N_CDT,P}(α::SCALARS,β::SCALARS,τ::SCALARS,
                        conn::CONN,nonl::NONL,stim::STIM,
                        pop_names::NTuple{P,String}) where {
            T,N_CDT,P,
            SCALARS<:NTuple{P,T},
            CONN<:PopInteractParam{P},#,<:MaybeNoOpComposite{T,AbstractConnectivityParameter{T,N_CDT}}},
            NONL<:PopAct{P,<:AbstractNonlinearity{T}},
            STIM<:PopActParam{P,<:AbstractStimulusParameter{T}}
        }
        new{T,N_CDT,P,SCALARS,CONN,NONL,STIM}(α,β,τ,conn,nonl,stim,pop_names)
    end
end
struct WCMSpatialAction{T,N_CDT,P,
        SCALARS<:NTuple{P,T},
        CONN<:PopInteract{P},#,<:AbstractConnectivityAction{T,N_CDT}},
        NONL<:PopAct{P},#,<:AbstractNonlinearity{T}},
        STIM<:PopAct{P}#,<:AbstractStimulusAction{T,N_CDT}}
        } <: AbstractSpaceAction{T,N_CDT}
    α::SCALARS
    β::SCALARS
    τ::SCALARS
    connectivity::CONN
    nonlinearity::NONL
    stimulus::STIM
    pop_names::NTuple{P,String}
end

function WCMSpatial{N_CDT,P}(; pop_names, α, 
        β, τ, 
        connectivity, nonlinearity, stimulus
       ) where {N_CDT,P}
    WCMSpatial{N_CDT,P}(
        α, β, τ,
        connectivity, nonlinearity,
        stimulus, pop_names
    )
end
function (wcm::WCMSpatial{T,N_CDT,P,SCALARS})(space::AbstractSpace{T}) where {T,N_CDT,P,SCALARS}
    conn = wcm.connectivity(space)
    nonl = wcm.nonlinearity
    stim = wcm.stimulus(space)
    CONN = typeof(conn)
    NONL = typeof(nonl)
    STIM = typeof(stim)
    return WCMSpatialAction{T,N_CDT,P,SCALARS,CONN,NONL,STIM}(wcm.α, wcm.β, wcm.τ,
    conn, nonl, stim, wcm.pop_names)
end
function (wcm::WCMSpatialAction{T,N_CDT,P})(dA,A,p,t) where {T,N_CDT,P}
    wcm.stimulus(dA, A, t)
    wcm.connectivity(dA, A, t)
    wcm.nonlinearity(dA, A, t)
    for i in 1:P
        dAi = population(dA,i); Ai = population(A,i)
        dAi .*= wcm.β[i] .* (1.0 .- Ai)
        dAi .+= -wcm.α[i] .* Ai
        dAi ./= wcm.τ[i]
    end
end
        
