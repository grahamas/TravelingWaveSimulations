export tw_metrics,
    SolitaryWaveformMetrics, WavefrontMetrics,
    ValuedSpace

const MaybeData{T} = Union{T,Missing}
abstract type AbstractWaveform{T_LOC,T_VAL} end
abstract type AbstractWaveformMetrics{T} end

struct Value{T_LOC,T_VAL}
    loc::T_LOC
    val::T_VAL
end
Value(loc::Int, val::T) where T = (@assert loc > 0; Value{Int,T}(loc,val))
Base.:(==)(a::Value, b::Value) = a.val == b.val
Base.:(==)(a::Value, b::Number) = a.val == b
Base.isless(a::Value, b::Value) = a.val < b.val
Base.isless(a::Value, b::Number) = a.val < b
Base.isless(a::Number, b::Value) = a < b.val

# function Base.minimum(arr::Array{<:Value})
#     val, idx = findmin([elt.val for elt in arr]) 
#     return arr[idx]
# end

function linear_interpolate((x1,x2), (y1,y2), x)
    @assert x1 <= x <= x2 || x2 <= x <= x1
    if x1 == x2
        return y1 + (y2 - y1) / 2 # FIXME: if flat, just split the difference
    end
    y1 + (x - x1) * ((y2 - y1) / (x2 - x1)) 
end


function linear_interpolate_loc((left,right)::Tuple{Value,Value}, val)
    loc = linear_interpolate((left.val,right.val), (left.loc,right.loc), val)
    Value(loc, val)
end
function linear_interpolate_val((left,right)::Tuple{Value,Value}, loc)
    val = linear_interpolate((left.loc,right.loc), (left.val,right.val), loc)
    Value(loc, val)
end

all_but_last(itr, n=1) = Iterators.take(itr, length(itr)-n)

struct ValuedSpace{T,C,N,AT<:AbstractArray{T,N},AC<:AbstractArray{C,N}} <: AbstractArray{T,N}
    values::AT
    coordinates::AC
    function ValuedSpace(values::AT, coordinates::AC) where {T,C,N,AT<:AbstractArray{T,N},AC<:AbstractArray{C,N}}
        @assert size(values) == size(coordinates)
        return new{T,C,N,AT,AC}(values, coordinates)
    end
end


Base.BroadcastStyle(::Type{<:ValuedSpace}) = Broadcast.ArrayStyle{ValuedSpace}()

Base.size(vs::ValuedSpace) = size(vs.values)
Base.getindex(vs::ValuedSpace, inds::Vararg{Int,N}) where N = vs.values[inds...]#ValuedSpace(vs.values[I], vs.coordinates[I])
getvalue(vs::ValuedSpace, idx::Union{CartesianIndex,Int}) = Value(vs.coordinates[idx], vs.values[idx])
#Base.getindex(vs::ValuedSpace, idx::Union{CartesianIndex,Int}) = ValuedSpace([vs.values[idx]], [vs.coordinates[idx]])
function Base.getindex(vs::ValuedSpace, fidx::AbstractFloat)
    lidx = findlast(vs.coordinates .< fidx)
    if lidx == length(vs)
        return getvalue(vs,lidx)
    end
    val = linear_interpolate_val((getvalue(vs,lidx), getvalue(vs,lidx+1)), fidx)
    return val
end
function getslice(vs::ValuedSpace, (left,right)::Tuple{Value,Value})
    return getslice(vs, (left.loc, right.loc))
end
function getslice(vs::ValuedSpace, (left_loc,right_loc)::Tuple{AbstractFloat,AbstractFloat})
    xs = vs.coordinates
    maybe_left = findlast(xs .< left_loc)
    left_idx = if maybe_left === nothing
        findfirst(xs .>= left_loc)
    else
        maybe_left + 1
    end
    maybe_right = findfirst(xs .> right_loc)
    right_idx = if maybe_right === nothing
        findlast(xs .<= right_loc)
    else
        maybe_right - 1
    end
    return ValuedSpace(vs.values[left_idx:right_idx], vs.coordinates[left_idx:right_idx])
end
Base.setindex!(vs::ValuedSpace, val, inds::Vararg{Int,N}) where N = vs.values[inds...] = val
Base.showarg(io::IO, vs::ValuedSpace, toplevel) = print(io, typeof(vs), ":\n\tValues: $(vs.values)\n\tCoordinates: $(vs.coordinates)")

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{ValuedSpace}}, ::Type{ElType}) where {ElType}
    # Scan the inputs for the ValuedSpace:
    vs = find_valuedspace(bc)
    # Use the char field of A to create the output
    ValuedSpace(similar(Array{ElType}, axes(bc)), vs.coordinates)
end

"`A = find_valuedspace(bc)` returns the first ValuedSpace among the arguments."
find_valuedspace(bc::Base.Broadcast.Broadcasted) = find_valuedspace(bc.args)
find_valuedspace(args::Tuple) = find_valuedspace(find_valuedspace(args[1]), Base.tail(args))
find_valuedspace(x) = x
find_valuedspace(a::ValuedSpace, rest) = a
find_valuedspace(::Any, rest) = find_valuedspace(rest)


function zero_crossing(A::ValuedSpace, i, j, atol=1e-4)
    if (A[i] <= 0 && A[j] >= 0) || (A[i] >= 0 && A[j] <= 0)
        return linear_interpolate_loc((getvalue(A,i),getvalue(A,j)), 0.0)
    elseif (A[i] < -atol && A[j] >= -atol) || (A[i] > atol && A[j] <= atol)
        return getvalue(A,j)
    else
        return nothing
    end
end

Base.diff(vs::ValuedSpace) = ValuedSpace(diff(vs.values) ./ diff(vs.coordinates), collect(all_but_last(vs.coordinates)) .+ (0.5 .* diff(vs.coordinates)))

@noinline function translate(val::Value, vs::ValuedSpace)
    vs[val.loc]
end
translate(::Nothing, args...) = nothing

struct Scored{OBJ,SCR}
    obj::OBJ
    score::SCR
end
Scored(obj::OBJ) where OBJ = Scored(obj, score(obj))

sigmoid(x::Real) = one(x) / (one(x) + exp(-x))
drop_proportion(apex::Value, nadir::Value) = min(abs(apex.val - nadir.val) / abs(apex.val), 1.0)
drop_duration(apex::Value, nadir::Value) = abs(apex.loc - nadir.loc)
drop_score(apex::Value{T}, nadir::Value{T}, θ) where T = drop_proportion(apex, nadir) * sigmoid((drop_duration(apex, nadir) - θ))
valid_score(score) = if isnan(score)
        return 0.0
    elseif score < 0
        # Confirm score range
        @warn "negative score"
        return 0.0
    else
        return score
    end 
# RIGHT WAVEFRONT FIXME
struct Wavefront{T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}} <: AbstractWaveform{T_LOC,T_VAL}
    left::V
    slope::V
    right::V
end
function translate(wf::Wavefront, args...)
    Wavefront(
        translate(wf.left, args...),
        translate(wf.slope, args...),
        translate(wf.right, args...)
        )
end
function score(wf::Wavefront, width_θ::T=5) where T
    apex = max(wf.left, wf.right)
    left_score = drop_score(apex, wf.left, width_θ)
    right_score = drop_score(apex, wf.right, width_θ)
    score = max(right_score, left_score) # TODO account for slope
    return valid_score(score)
end
struct WavefrontMetrics{T} <: AbstractWaveformMetrics{T}
    left_height::T
    right_height::T
    slope::T
    slope_loc::T
    width::T
end
function metrics(wf::Wavefront)
    WavefrontMetrics(
        wf.left.val,
        wf.right.val,
        wf.slope.val,
        wf.slope.loc,
        wf.right.loc - wf.left.loc        
    ) 
end

struct SolitaryWaveform{T_LOC,T_VAL,V<:Value{T_LOC,T_VAL}} <: AbstractWaveform{T_LOC,T_VAL}
    left::V
    left_slope::V
    apex::V
    right_slope::V
    right::V
end
function SolitaryWaveform(left::WF, right::WF) where {T_LOC,T_VAL,V<:Value{T_LOC,T_VAL},WF<:Wavefront{T_LOC,T_VAL,V}}
    SolitaryWaveform{T_LOC,T_VAL,V}(
        left.left,
        left.slope,
        left.right, # == right.left
        right.slope,
        right.right
    )
end
function score(sw::SolitaryWaveform, width_θ::T=10) where T
    apex = sw.apex
    left_score = drop_score(apex, sw.left, width_θ)
    right_score = drop_score(apex, sw.right, width_θ)
    score = right_score * left_score # TODO account for slope
    return valid_score(score)
end

struct SolitaryWaveformMetrics{T} <: AbstractWaveformMetrics{T}
    left_baseline_height::T
    right_baseline_height::T
    right_slope_loc::T
    apex_height::T
    apex_loc::T
    between_slopes_width::T
end
function metrics(sw::SolitaryWaveform)
    SolitaryWaveformMetrics(
        sw.left.val,
        sw.right.val,
        sw.right_slope.loc,
        sw.apex.val,
        sw.apex.loc,
        sw.right.loc - sw.left.loc
    )
end
    
choose_best(::Nothing, s::Scored) = s
choose_best(s1::Sd, s2::Sd) where {OBJ,Sd <: Scored{OBJ}} = (s1.score >= s2.score) ? s1 : s2
function choose_best(arr1::AbstractArray{<:Sd}, arr2::AbstractArray{<:Sd}) where {OBJ,Sd <: Scored{OBJ}}
    sum((elt) -> elt.score, arr1) > sum((elt) -> elt.score, arr2) ? arr1 : arr2
end




# Divide between peaks of second derivative (inflection points)
# Center on peaks of first derivative

function detect_all_fronts(valued_space::ValuedSpace)
    d_values = diff(valued_space)
    dd_values = diff(d_values)
    fronts = Wavefront{Float64,Float64}[]
    left_boundary = getvalue(valued_space, 1)
    steepest_slope = nothing
    for idx=all_but_last(eachindex(d_values),2)
        right_boundary = translate(zero_crossing(d_values, idx, idx+1, 1e-4), valued_space)
        if right_boundary !== nothing
            this_front = getslice(d_values, (left_boundary, right_boundary))
            _, midx = findmax(abs.(this_front))
            steepest_slope = getvalue(this_front, midx)
            push!(fronts, Wavefront(left_boundary,
                                    steepest_slope,
                                    right_boundary)
            )
            steepest_slope = nothing
            left_boundary = right_boundary
            right_boundary = nothing
        end
    end
    right_boundary = getvalue(valued_space, length(valued_space))
    this_front = getslice(d_values, (left_boundary, right_boundary))
    _, midx = findmax(abs.(this_front))
    steepest_slope = getvalue(this_front, midx)
    push!(fronts, Wavefront(left_boundary,
                            steepest_slope,
                            right_boundary)
    )
    return fronts
end

eat_left(left, right) = Wavefront(left.left, right.slope, right.right)
eat_right(left, right) = Wavefront(left.left, left.slope, right.right)
function consolidate_fronts(fronts::AbstractVector{<:Wavefront{Float64}}, vs::ValuedSpace, slope_min=1e-4)
    if length(fronts) == 0
        return fronts
    end
    new_fronts = Wavefront[]
    first_suff_dx = 1
    while abs(fronts[first_suff_dx].slope.val) < slope_min
        if first_suff_dx < length(fronts)
            first_suff_dx += 1
        else
            break
        end
    end
    
    current_front = eat_left(fronts[1], fronts[first_suff_dx])
    if length(fronts) > first_suff_dx
        for i in first_suff_dx+1:length(fronts)
            next_front = fronts[i]
            if abs(next_front.slope.val) < slope_min
                current_front = eat_right(current_front, next_front)
            else
                push!(new_fronts, current_front)
                current_front = next_front
            end
        end
    end
    push!(new_fronts, current_front)
    return new_fronts
end

# TODO deal with multipop
substantial_fronts(multipop::AbstractMatrix, xs::AbstractVector, slope_min=1e-4) = substantial_fronts(multipop[:,1], xs)
function substantial_fronts(frame::AbstractVector, xs::AbstractVector, slope_min=1e-4)
    vs = ValuedSpace(frame, xs)
    all_fronts = detect_all_fronts(vs)
    consolidated = consolidate_fronts(all_fronts, vs, slope_min)
    return consolidated
end

function scored_rightmost(::Type{<:WavefrontMetrics}, fronts::AbstractArray{<:Wavefront})#::Scored{<:Wavefront{T,T},T} where T
    if length(fronts) == 0
        return fronts
    else
        return Scored(fronts[end])
    end
end
function scored_rightmost(::Type{<:SolitaryWaveformMetrics}, fronts::AbstractArray{<:Wavefront})
    for idx = length(fronts):-1:2
        if fronts[idx].slope.val < 0 && fronts[idx-1].slope.val > 0
            return Scored(SolitaryWaveform(fronts[idx-1], fronts[idx]))
        end
    end
    return nothing
end
scored_rightmost_wavefront(multipop::AbstractArray{T,2}, xs) where T = scored_rightmost_wavefront(multipop[:,1], xs)

function metrics_df(metrics_type::Type{<:AbstractWaveformMetrics}, scored_wave_arr::AbstractArray, t)
    waveform_metrics = [metrics(s.obj) for s in scored_wave_arr if s !== nothing]
    scores = [s.score for s in scored_wave_arr if s !== nothing]
    if length(scores) < 3 || mean(scores) < 1e-2 # FIXME magic number
        return nothing
    end
    wave_metric_syms = fieldnames(metrics_type)
    wave_metrics_df = DataFrame(Dict(zip(wave_metric_syms, [[getproperty(wave, sym) for wave in waveform_metrics] for sym in wave_metric_syms])))
    wave_metrics_df.score = scores
    wave_metrics_df.t = t[scored_wave_arr .!== nothing]
    return wave_metrics_df
end

# TODO: deal with non-traveling waves (return nothing)
function traveling_wave_metrics_linear_regression(wave_metrics_df::DataFrame)
    scores = wave_metrics_df.score
    wave_metrics_syms = filter(name -> name ∉ [:score, :t], names(wave_metrics_df))
    mostly_extant_cols = filter(wave_metrics_syms) do colname
        (count(ismissing, wave_metrics_df[:,colname])/size(wave_metrics_df,1) < 0.01 && # FIXME magic numbers
            sum(wave_metrics_df[:,colname]) > 5 &&
            length(wave_metrics_df[:,colname]) > 4) # FIXME magic numbers
    end
    results = map(mostly_extant_cols) do metric_sym
        fmla = Term(metric_sym) ~ Term(:t) + ConstantTerm(1)
        glm_fit = try
            glm(fmla, wave_metrics_df, Normal(), IdentityLink(); wts=scores)
        catch e
            @show metric_sym
            @show wave_metrics_df[:,metric_sym]
            nothing
        end
        return metric_sym => glm_fit
    end |> Dict
    return (results, scores, wave_metrics_df)
end
traveling_wave_metrics_linear_regression(::Nothing) = (nothing, nothing, nothing)

function tw_metrics(metrics_type::Type{<:AbstractWaveformMetrics}, exec::Execution)
    u = exec.solution.u
    t = timepoints(exec)
    x = [x[1] for x in space(exec).arr]
    tw_metrics(metrics_type, u, t, x)
end
function tw_metrics(metrics_type::Type{<:AbstractWaveformMetrics}, frames, t, x)
    fronts = substantial_fronts.(frames, Ref(x))
    traveling_wave_metrics_linear_regression(metrics_df(metrics_type, scored_rightmost.(Ref(metrics_type), fronts), t))
end


@recipe function f(vs::ValuedSpace)
    (vs.coordinates, vs.values)
end

@recipe function f(wf_arr::Array{<:Wavefront}, vs=nothing)
    seriestype := :scatter
    @series begin
        color := :green
        marker := :star
        markersize := 10
        [wf.left.loc for wf in wf_arr], [wf.left.val for wf in wf_arr]
    end
    @series begin
        color := :red
        [wf.right.loc for wf in wf_arr], [wf.right.val for wf in wf_arr]
    end
    @series begin
        if vs !== nothing
            color := :blue
            [wf.slope.loc for wf in wf_arr], [vs[wf.slope.loc].val for wf in wf_arr]
        else
            color := :blue
            [wf.slope.loc for wf in wf_arr], [wf.slope.val for wf in wf_arr]
        end
    end
end
