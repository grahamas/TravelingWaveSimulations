
abstract type AbstractStimulus{T} <: AbstractParameter{T} end

function update(calc_arr::AA, new_arr::AbstractArray{S,1}, space::Space{T}) where {T, S <: AbstractStimulus{T}, CS <: CalculatedType{S}, AA<:AbstractArray{CS,1}}
    [calc_arr[i].stimulus != new_arr[i] ? Calculated(new_arr[i], space) : calc_arr[i] for i in CartesianIndices(calc_arr)]
end

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

#region SharpBumpStimulus

struct SharpBumpStimulus{T} <: AbstractStimulus{T}
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

function Calculated(sbs::SharpBumpStimulus{T}, space::PopSegment{T}) where T
    calculated_space = Calculated(space)
    on_frame = make_sharp_bump_frame(calculated_space.value[:,1], sbs.width, sbs.strength)
    off_frame = zero(on_frame)
    onset = sbs.window[1]
    offset = sbs.window[2]
    return CalculatedSharpBumpStimulus{T}(sbs, space, onset, offset, on_frame, off_frame)
end

struct CalculatedSharpBumpStimulus{T} <: CalculatedType{SharpBumpStimulus{T}}
    stimulus::SharpBumpStimulus{T}
    space::PopSegment{T}
    onset::T
    offset::T
    on_frame::Array{T,1}
    off_frame::Array{T,1}
end

function make_sharp_bump_frame(mesh_coords::AbstractArray{DistT,1}, width::DistT, strength::T) where {DistT,T}
    mid_dx = floor(Int, size(mesh_coords, 1) / 2)
    mid_point = mesh_coords[mid_dx,1]
    frame = zero(mesh_coords)
    half_width = width / 2      # using truncated division
    start_dx = findfirst(mesh_coords .>= mid_point - half_width)
    stop_dx = findlast(mesh_coords .<= mid_point + half_width)
    frame[start_dx:stop_dx] .= strength
    return frame
end

function stimulate!(val::AT, sharp_bump::CalculatedSharpBumpStimulus, t::T) where {T, AT <: AbstractArray{T,1}}
    if sharp_bump.onset <= t < sharp_bump.offset
        val .= sharp_bump.on_frame
    else
        val .= sharp_bump.off_frame
    end
end
function stimulate_add!(val::AT, sharp_bump::CalculatedSharpBumpStimulus, t::T) where {T, AT<: AbstractArray{T,1}}
    if sharp_bump.onset <= t < sharp_bump.offset
        val .+= sharp_bump.on_frame
    end
end

#endregion

#region NoisySharpBumpStimulus

struct NoisySharpBumpStimulus{T} <: AbstractStimulus{T}
    noise::GaussianNoiseStimulus{T}
    bump::SharpBumpStimulus{T}
end
function NoisySharpBumpStimulus{T}(; strength::T, window::Tuple{T,T}, width::T, SNR::T) where T
    NoisySharpBumpStimulus{T}(
            GaussianNoiseStimulus{T}(SNR = SNR),
            SharpBumpStimulus{T}(strength=strength, window=window,
                width=width)
        )
end
struct CalculatedNoisySharpBumpStimulus{T} <: CalculatedType{NoisySharpBumpStimulus{T}}
    stimulus::NoisySharpBumpStimulus{T}
    calc_noise::CalculatedGaussianNoiseStimulus{T}
    calc_bump::CalculatedSharpBumpStimulus{T}
end
Calculated(nsbs::NoisySharpBumpStimulus{T}, space::PopSegment{T}) where T = CalculatedNoisySharpBumpStimulus{T}(nsbs, Calculated(nsbs.noise, space), Calculated(nsbs.bump, space))

function stimulate!(val::AT, stim_obj::CalculatedNoisySharpBumpStimulus{T},t::T) where {T, AT<: AbstractArray{T,1}}
    stimulate!(val, stim_obj.calc_noise,t)
    stimulate_add!(val, stim_obj.calc_bump,t)
end

#endregion


CalculatedTypes.get_value(c::CalculatedGaussianNoiseStimulus{T}) where T = c
CalculatedTypes.get_value(c::CalculatedNoisySharpBumpStimulus{T}) where T = c
CalculatedTypes.get_value(c::CalculatedSharpBumpStimulus{T}) where T = c
