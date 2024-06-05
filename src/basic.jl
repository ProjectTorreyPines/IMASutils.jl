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