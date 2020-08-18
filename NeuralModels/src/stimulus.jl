
### Transient bumps ###
# Subtypes of TransientBumpStimulusParameter generate TransientBumpStimulusActions (NOT subtypes)
abstract type AbstractTransientBumpStimulusParameter{T} <: AbstractStimulusParameter{T} end
struct TransientBumpStimulusAction{T,N,FRAME,WINDOWS} <: AbstractStimulusAction{T,N}
    baseline::T
    bump_frame::FRAME
    time_windows::WINDOWS
    function TransientBumpStimulusAction(baseline::T, bump_frame::FRAME,time_windows::WINDOWS) where {T,N,FRAME<:AbstractArray{T,N},WINDOWS}
        new{T,N,FRAME,WINDOWS}(baseline,bump_frame,time_windows)
    end
end
function (bump_param::AbstractTransientBumpStimulusParameter{T})(space::AbstractSpace{T,N}) where {T,N}
    bump_frame = on_frame(bump_param, space)
    TransientBumpStimulusAction(bump_param.baseline, bump_frame, bump_param.time_windows)
end
function (bump::TransientBumpStimulusAction{T,N})(val::AbstractArray, ignored, t) where {T,N}
    for window in bump.time_windows
        if window[1] <= t < window[2]
            val .+= bump.bump_frame
        else
            val .+= bump.baseline
        end
    end
end


struct CircleStimulusParameter{T} <: AbstractTransientBumpStimulusParameter{T}
    radius::T
    strength::T
    time_windows::Array{Tuple{T,T},1}
    center::Union{NTuple,Nothing}
    baseline::T
end
center(sbs::CircleStimulusParameter) = sbs.center
export center

function CircleStimulusParameter(; strength, radius,
        duration=nothing, time_windows=nothing, center=nothing, baseline=0.0)
    if time_windows == nothing
        return CircleStimulusParameter(radius, strength, [(zero(typeof(strength)), duration)], center, baseline)
    else
        @assert duration == nothing
        return CircleStimulusParameter(radius, strength, time_windows, center, baseline)
    end
end
distance(x1::NTuple{N},x2::NTuple{N}) where N = sqrt(sum((x1 .- x2) .^ 2))
function on_frame(sbs::CircleStimulusParameter{T}, space::AbstractSpace{T,N_ARR,N_CDT}) where {T,N_ARR,N_CDT}
    coords = coordinates(space)
    frame = zero(space) .+ sbs.baseline
    center_coordinates = sbs.center == nothing ? Tuple(zeros(T,N_CDT)) : sbs.center
    distances = distance.(coords, Ref(center_coordinates))
    on_center = (distances .< sbs.radius) .| (distances .≈ sbs.radius)
    frame[on_center] .= sbs.strength
    return frame
end

struct RectangleStimulusParameter{T,N} <: AbstractTransientBumpStimulusParameter{T}
    widths::NTuple{N,T}
    strength::T
    time_windows::Array{Tuple{T,T},1}
    center::NTuple{N,T}
    baseline::T
end
center(sbs::RectangleStimulusParameter) = sbs.center

function RectangleStimulusParameter(; strength, widths,
                                    duration=nothing, time_windows=nothing, center=zero.(widths), baseline=0.0)
    if time_windows == nothing
        return RectangleStimulusParameter(widths, strength, [(zero(typeof(strength)), duration)], center, baseline)
    else
        @assert duration == nothing
        return RectangleStimulusParameter(widths, strength, time_windows, center, baseline)
    end
end
function on_frame(sbs::RectangleStimulusParameter{T,N_CDT}, space::AbstractSpace{T,N_ARR,N_CDT}) where {T,N_ARR,N_CDT}
    coords = coordinates(space)
    frame = zero(space) .+ sbs.baseline
    center_coordinates = sbs.center == nothing ? Tuple(zeros(T,N_CDT)) : sbs.center
    on_center = map(coords) do coord
        diffs = abs.(coord .- sbs.center)
        all((diffs .< sbs.widths) .| (diffs .≈ sbs.widths))
    end
    frame[on_center] .= sbs.strength
    return frame
end
################################
