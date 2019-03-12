
abstract type AbstractStimulus{T,D} <: AbstractParameter{T} end

function update(calc_arr::AA, new_arr::AbstractArray{S,1}, space::AbstractSpace{T}) where {T, S <: AbstractStimulus{T}, AA<:AbstractArray{<:CalculatedType,1}}
    [calc_arr[i].stimulus != new_arr[i] ? Calculated(new_arr[i], space) : calc_arr[i] for i in CartesianIndices(calc_arr)]
end

distance(x,y) = sum((x .- y).^2)

#region NoStimulus

struct NoStimulus{T,N} <: AbstractStimulus{T,N} end
struct CalculatedNoStimulus{T,N} <: CalculatedType{NoStimulus{T,N}}
    stimulus::NoStimulus{T,N}
end
function Calculated(stimulus::NoStimulus{T,N}, args...) where {T,N}
    CalculatedNoStimulus{T,N}(stimulus)
end
function stimulate!(val::AT, calcnostim::CalculatedNoStimulus{T,N}, t::T) where {T,N,AT<: AbstractArray{T,N}}
    val[:] .= 0
end

#endregion

#region GaussianNoiseStimulus

function gaussian_noise!(val::AT, mean::T, sd::T) where {T, AT<:AbstractArray{T}} # assumes signal power is 0db
    randn!(val)
    val .*= sd
    val .+= mean
end

struct GaussianNoiseStimulus{T,N} <: AbstractStimulus{T,N}
    mean::T
    SNR::T
end

function GaussianNoiseStimulus{T,N}(; SNR::T=0.0, mean::T=0.0) where {T,N}
    GaussianNoiseStimulus{T,N}(mean, SNR)
end

struct CalculatedGaussianNoiseStimulus{T,N} <: CalculatedType{GaussianNoiseStimulus{T,N}}
    stimulus::GaussianNoiseStimulus{T,N}
    calc_space::CalculatedType{<:AbstractSpace{T,N}}
    mean::T
    sd::T
end


function Calculated(wns::GaussianNoiseStimulus{T,N}, calc_space::CalculatedType{<:AbstractSpace{T,N}}) where {T,N}
    sd = sqrt(1/10 ^ (wns.SNR / 10))
    CalculatedGaussianNoiseStimulus{T,N}(wns, calc_space, wns.mean, sd)
end

function stimulate!(val::AT, wns::CalculatedGaussianNoiseStimulus{T}, t::T) where {T,AT<: AbstractArray{T}}
    gaussian_noise!(val, wns.mean, wns.sd) # Not actually time dependent
end

#endregion


abstract type TransientBumpStimulus{T,N} <: AbstractStimulus{T,N} end

function stimulate!(val::AT, bump::CTBS, t::T) where {T, N, AT <: AbstractArray{T,N}, CTBS<:CalculatedType{<:TransientBumpStimulus{T,N}}}
    if bump.onset <= t < bump.offset
        val .= bump.on_frame
    else
        val .= bump.off_frame
    end
end
function stimulate_add!(val::AT, bump::CTBS, t::T) where {T, N, AT<: AbstractArray{T,N}, CTBS<:CalculatedType{<:TransientBumpStimulus{T,N}}}
    if bump.onset <= t < bump.offset
        val .+= bump.on_frame
    end
end



#region SharpBumpStimulus

struct SharpBumpStimulus{T,N} <: TransientBumpStimulus{T,N}
    width::T
    strength::T
    time_window::Tuple{T,T}
    center::NTuple{N,T}
end

function SharpBumpStimulus{T,N}(; strength::T=nothing, width::T=nothing,
        duration=nothing, time_window=nothing, center=NTuple{N,T}(zero(T) for i in 1:N)) where {T,N}
    if time_window == nothing
        return SharpBumpStimulus{T,N}(width, strength, (zero(T), duration), center)
    else
        @assert duration == nothing
        return SharpBumpStimulus{T,N}(width, strength, time_window, center)
    end
end

# function make_bump_frame(sbs::SharpBumpStimulus, mesh_coords::AbstractArray{DistT,1}) where {DistT}
#     mid_dx = floor(Int, size(mesh_coords, 1) / 2)
#     mid_point = mesh_coords[mid_dx,1]
#     frame = zero(mesh_coords)
#     half_width = sbs.width / 2      # using truncated division
#     start_dx = findfirst(mesh_coords .>= mid_point - half_width)
#     stop_dx = findlast(mesh_coords .<= mid_point + half_width)
#     frame[start_dx:stop_dx] .= sbs.strength
#     return frame
# end

function make_bump_frame(sbs::SharpBumpStimulus, mesh::CalculatedType{<:Pops}) where {DistT}
    pop = one_pop(mesh)
    frame = zero(pop)
    half_width = sbs.width / 2.0
    frame[distance.(pop, Ref(sbs.center)) .<= half_width] .= sbs.strength
    return frame
end

struct CalculatedSharpBumpStimulus{T,N} <: CalculatedType{SharpBumpStimulus{T,N}}
    stimulus::SharpBumpStimulus{T}
    calc_space::CalculatedType{<:Pops}
    onset::T
    offset::T
    on_frame::Array{T,1}
    off_frame::Array{T,1}
end

function Calculated(tbs::SharpBumpStimulus{T,N}, calc_space::CalculatedType{<:Pops{P,T,N,<:AbstractSpace{T,N}}}) where {P,T,N}
    on_frame = make_bump_frame(tbs, calc_space)
    off_frame = zero(on_frame)
    onset = tbs.time_window[1]
    offset = tbs.time_window[2]
    return CalculatedSharpBumpStimulus{T,N}(tbs, calc_space, onset, offset, on_frame, off_frame)
end

#endregion

#region Sech2BumpStimulus
# struct Sech2BumpStimulus{T} <: TransientBumpStimulus{T}
#     width::T
#     strength::T
#     window::Tuple{T,T}
# end
#
# function sech2(x, A, a)
#     A * (sech(a * x))^2
# end
#
# function Sech2BumpStimulus{T,N}(; strength::T=nothing, width::T=nothing,
#         duration=nothing, time_window=nothing, center::NTuple{N,T}=NTuple{N,T}(zero(T) for i in 1:N)) where {T,N}
#     if time_window == nothing
#         return Sech2BumpStimulus{T}(width, strength, (zero(T), duration), center)
#     else
#         @assert duration == nothing
#         return Sech2BumpStimulus{T}(width, strength, time_window, center)
#     end
# end
#
# function make_bump_frame(s2bs::Sech2BumpStimulus, mesh_coords::AbstractArray{DistT,T}) where {DistT, T}
#     frame = sech2.(distance.(mesh_coords, Ref(s2bs.center)), s2bs.strength, s2bs.width)
#     return frame
# end
#
#
# struct CalculatedSech2BumpStimulus{T} <: CalculatedType{Sech2BumpStimulus{T}}
#     stimulus::Sech2BumpStimulus{T}
#     space::Pops
#     onset::T
#     offset::T
#     on_frame::Array{T,1}
#     off_frame::Array{T,1}
# end
#
# function Calculated(tbs::Sech2BumpStimulus{T}, space::Pops{P,T,1,<:AbstractSpace{T,1}}) where {P,T}
#     calculated_space = Calculated(space)
#     on_frame = make_bump_frame(tbs, calculated_space.value[:,1])
#     off_frame = zero(on_frame)
#     onset = tbs.window[1]
#     offset = tbs.window[2]
#     return CalculatedSech2BumpStimulus{T}(tbs, space, onset, offset, on_frame, off_frame)
# end

#endregion

#region NoisySharpBumpStimulus

struct NoisyStimulus{T,N,STIM} <: AbstractStimulus{T,N}
    noise::GaussianNoiseStimulus{T,N}
    stim::STIM
end
function NoisyStimulus{T,N}(; stim_type::Type=SharpBumpStimulus{T,N}, SNR::T, mean::T=0.0, kwargs...) where {T,N}
    NoisyStimulus{T,N,stim_type}(
            GaussianNoiseStimulus{T,N}(SNR = SNR, mean=mean),
            stim_type(; kwargs...)
        )
end
struct CalculatedNoisyStimulus{T,N,STIM} <: CalculatedType{NoisyStimulus{T,N,STIM}}
    stimulus::NoisyStimulus{T,N,STIM}
    calc_noise::CalculatedGaussianNoiseStimulus{T,N}
    calc_stim::CalculatedType{STIM}
end
Calculated(ns::NoisyStimulus{T,N,STIM}, calc_space::CalculatedType{<:Pops{P,T,N,<:AbstractSpace{T,N}}}) where {P,T,N,STIM} = CalculatedNoisyStimulus{T,N,STIM}(ns, Calculated(ns.noise, calc_space), Calculated(ns.stim, calc_space))

function stimulate!(val::AT, stim_obj::CalculatedNoisyStimulus{T},t::T) where {T, AT<: AbstractArray{T}}
    stimulate!(val, stim_obj.calc_noise,t)
    stimulate_add!(val, stim_obj.calc_stim,t)
end

#endregion


CalculatedTypes.get_value(c::CalculatedType{<:AbstractStimulus}) = c
