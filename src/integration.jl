"""
    trapz(x::AbstractVector{S}, y::AbstractVector{T}) where {S<:Real, T<:Real}

Numerically integrates the values in `y` with respect to `x` using the trapezoidal rule.
"""
function trapz(x::AbstractVector{S}, y::AbstractVector{T}) where {S<:Real, T<:Real}
    N = length(x)
    @assert N>=2
    @assert length(y) == N
    t = zero(promote_type(T, S))
    @inbounds begin
        t += (x[2] - x[1]) * y[1]
        for k in eachindex(x)[2:end-1]
            t += y[k] * (x[k+1] - x[k-1])
        end
        t += (x[end] - x[end-1]) * y[end]
    end
    return 0.5 * t
end

"""
    trapz(x::AbstractVector{S}, f::T) where {S<:Real, T<:Function}

Numerically integrates `f(k::Int, xx::S)` for `xx = x[k]` with respect to `x` using the trapezoidal rule.
"""
function trapz(x::AbstractVector{S}, f::T) where {S<:Real, T<:Function}
    N = length(x)
    @assert N>=2
    y = k -> f(k, x[k])

    # get type of solution, and also check bounds
    y1, yN = y(1), y(N)
    t = zero(promote_type(S, typeof(y1)))

    @inbounds begin
        t += (x[2] - x[1]) * y1
        for k in eachindex(x)[2:end-1]
            t += y(k) * (x[k+1] - x[k-1])
        end
        t += (x[end] - x[end-1]) * yN
    end
    return 0.5 * t
end

"""
    cumtrapz(x::AbstractVector{S}, y::AbstractVector{T}) where {S<:Real, T<:Real}

Computes the cumulative integral of the values in `y` with respect to `x` using the trapezoidal rule.
"""
function cumtrapz(x::AbstractVector{S}, y::AbstractVector{T}) where {S<:Real, T<:Real}
    N = length(x)
    @assert N>=2
    @assert length(y) == N
    retarr = Vector{promote_type(S, T)}(undef, N)
    @inbounds begin
        retarr[1] = 0.0
        for k in eachindex(x)[2:end]
            retarr[k] = retarr[k-1] + (y[k] + y[k-1]) * (x[k] - x[k-1])
        end
        retarr .*= 0.5
    end
    return retarr
end

"""
    cumtrapz(x::AbstractVector{S}, f::T) where {S<:Real, T<:Function}

Computes the cumulative integral of `f(k::Int, xx::S)` for `xx = x[k]` with respect to `x` using the trapezoidal rule.
"""
function cumtrapz(x::AbstractVector{S}, f::T) where {S<:Real, T<:Function}
    N = length(x)
    @assert N>=2
    y = k -> f(k, x[k])

    # get type of solution, and also check bounds
    y1, _ = y(1), y(N)
    retarr = Vector{promote_type(S, typeof(y1))}(undef, N)
    @inbounds begin
        retarr[1] = 0.0
        for k in eachindex(x)[2:end]
            retarr[k] = retarr[k-1] + (y(k) + y(k-1)) * (x[k] - x[k-1])
        end
        retarr .*= 0.5
    end
    return retarr
end