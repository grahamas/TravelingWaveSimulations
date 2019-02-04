
abstract type AbstractStimulus{T} <: AbstractParameter{T} end

function update(calc_arr::AA, new_arr::AbstractArray{S,1}, space::Space{T}) where {T, S <: AbstractStimulus{T}, CS <: CalculatedType{S}, AA<:AbstractArray{CS,1}}
    [calc_arr[i].stimulus != new_arr[i] ? Calculated(new_arr[i], space) : calc_arr[i] for i in CartesianIndices(calc_arr)]
end

#region NoStimulus

struct NoStimulus{T} <: AbstractStimulus{T} end
struct CalculatedNoStimulus{T} <: CalculatedType{NoStimulus{T}}
    stimulus::NoStimulus{T}
end
function Calculated(stimulus::NoStimulus{T}, args...) where T
    CalculatedNoStimulus{T}(stimulus)
end
function stimulate!(val::AT, calcnostim::CalculatedNoStimulus{T}, t::T) where {T,AT<: AbstractArray{T}}
    val[:] .= 0
end

#endregion

#region GaussianNoiseStimulus

function gaussian_noise!(val::AT, mean::T, sd::T) where {T, AT<:AbstractArray{T}} # assumes signal power is 0db
    randn!(val)
    val .*= sd
    val .+= mean
end

struct GaussianNoiseStimulus{T} <: AbstractStimulus{T}
    mean::T
    SNR::T
end

function GaussianNoiseStimulus{T}(; SNR::T=0.0, mean::T=0.0) where {T}
    GaussianNoiseStimulus{T}(mean, SNR)
end

struct CalculatedGaussianNoiseStimulus{T} <: CalculatedType{GaussianNoiseStimulus{T}}
    stimulus::GaussianNoiseStimulus{T}
    space::Space{T}
    mean::T
    sd::T
end


function Calculated(wns::GaussianNoiseStimulus{T}, space::Space{T,N}) where {T,N}
    sd = sqrt(1/10 ^ (wns.SNR / 10))
    CalculatedGaussianNoiseStimulus{T}(wns, space, wns.mean, sd)
end

function stimulate!(val::AT, wns::CalculatedGaussianNoiseStimulus{T}, t::T) where {T,AT<: AbstractArray{T}}
    gaussian_noise!(val, wns.mean, wns.sd) # Not actually time dependent
end

#endregion


abstract type TransientBumpStimulus{T} <: AbstractStimulus{T} end

function stimulate!(val::AT, bump::CTBS, t::T) where {T, AT <: AbstractArray{T,1}, CTBS<:CalculatedType{<:TransientBumpStimulus}}
    if bump.onset <= t < bump.offset
        val .= bump.on_frame
    else
        val .= bump.off_frame
    end
end
function stimulate_add!(val::AT, bump::CTBS, t::T) where {T, AT<: AbstractArray{T,1}, CTBS<:CalculatedType{<:TransientBumpStimulus}}
    if bump.onset <= t < bump.offset
        val .+= bump.on_frame
    end
end







#region SharpBumpStimulus

struct SharpBumpStimulus{T} <: TransientBumpStimulus{T}
    width::T
    strength::T
    window::Tuple{T,T}
end

function SharpBumpStimulus{T}(; strength::T=nothing, width::T=nothing,
        duration=nothing, window=nothing) where {T}
    if window == nothing
        return SharpBumpStimulus{T}(width, strength, (zero(T), duration))
    else
        @assert duration == nothing
        return SharpBumpStimulus{T}(width, strength, window)
    end
end

function SharpBumpStimulus(p)
    SharpBumpStimulus(p[:(Stimulus.width)], p[:(Stimulus.strength)], p[:(Stimulus.window)])
end



function make_bump_frame(sbs::SharpBumpStimulus, mesh_coords::AbstractArray{DistT,1}) where {DistT}
    mid_dx = floor(Int, size(mesh_coords, 1) / 2)
    mid_point = mesh_coords[mid_dx,1]
    frame = zero(mesh_coords)
    half_width = sbs.width / 2      # using truncated division
    start_dx = findfirst(mesh_coords .>= mid_point - half_width)
    stop_dx = findlast(mesh_coords .<= mid_point + half_width)
    frame[start_dx:stop_dx] .= sbs.strength
    return frame
end


struct CalculatedSharpBumpStimulus{T} <: CalculatedType{SharpBumpStimulus{T}}
    stimulus::SharpBumpStimulus{T}
    space::PopSegment{T}
    onset::T
    offset::T
    on_frame::Array{T,1}
    off_frame::Array{T,1}
end

function Calculated(tbs::SharpBumpStimulus{T}, space::PopSegment{T}) where T
    calculated_space = Calculated(space)
    on_frame = make_bump_frame(tbs, calculated_space.value[:,1])
    off_frame = zero(on_frame)
    onset = tbs.window[1]
    offset = tbs.window[2]
    return CalculatedSharpBumpStimulus{T}(tbs, space, onset, offset, on_frame, off_frame)
end

#endregion

#region Sech2BumpStimulus
struct Sech2BumpStimulus{T} <: TransientBumpStimulus{T}
    width::T
    strength::T
    window::Tuple{T,T}
end

function sech2(x, A, a)
    A * (sech(a * x))^2
end

function Sech2BumpStimulus{T}(; strength::T=nothing, width::T=nothing,
        duration=nothing, window=nothing) where {T}
    if window == nothing
        return Sech2BumpStimulus{T}(width, strength, (zero(T), duration))
    else
        @assert duration == nothing
        return Sech2BumpStimulus{T}(width, strength, window)
    end
end

function make_bump_frame(s2bs::Sech2BumpStimulus, mesh_coords::AbstractArray{DistT,1}) where {DistT}
    mid_dx = floor(Int, size(mesh_coords, 1) / 2)
    mid_point = mesh_coords[mid_dx,1]
    frame = sech2.(mesh_coords, s2bs.strength, s2bs.width)
    return frame
end


struct CalculatedSech2BumpStimulus{T} <: CalculatedType{Sech2BumpStimulus{T}}
    stimulus::Sech2BumpStimulus{T}
    space::PopSegment{T}
    onset::T
    offset::T
    on_frame::Array{T,1}
    off_frame::Array{T,1}
end

function Calculated(tbs::Sech2BumpStimulus{T}, space::PopSegment{T}) where T
    calculated_space = Calculated(space)
    on_frame = make_bump_frame(tbs, calculated_space.value[:,1])
    off_frame = zero(on_frame)
    onset = tbs.window[1]
    offset = tbs.window[2]
    return CalculatedSech2BumpStimulus{T}(tbs, space, onset, offset, on_frame, off_frame)
end

#endregion

#region NoisySharpBumpStimulus

struct NoisyStimulus{T,STIM} <: AbstractStimulus{T}
    noise::GaussianNoiseStimulus{T}
    stim::STIM
end
function NoisyStimulus{T}(; stim_type::Type=SharpBumpStimulus{T}, SNR::T, mean::T=0.0, kwargs...) where T
    NoisyStimulus{T,stim_type}(
            GaussianNoiseStimulus{T}(SNR = SNR, mean=mean),
            stim_type(; kwargs...)
        )
end
struct CalculatedNoisyStimulus{T,STIM} <: CalculatedType{NoisyStimulus{T,STIM}}
    stimulus::NoisyStimulus{T,STIM}
    calc_noise::CalculatedGaussianNoiseStimulus{T}
    calc_stim::CalculatedType{STIM}
end
Calculated(ns::NoisyStimulus{T,STIM}, space::PopSegment{T}) where {T,STIM} = CalculatedNoisyStimulus{T,STIM}(ns, Calculated(ns.noise, space), Calculated(ns.stim, space))

function stimulate!(val::AT, stim_obj::CalculatedNoisyStimulus{T},t::T) where {T, AT<: AbstractArray{T,1}}
    stimulate!(val, stim_obj.calc_noise,t)
    stimulate_add!(val, stim_obj.calc_stim,t)
end

#endregion


CalculatedTypes.get_value(c::CalculatedType{<:AbstractStimulus}) = c
