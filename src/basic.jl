"""
    argmin_abs(x::AbstractVector{<:Real}, x0::Real)

Non-allocating version of `argmin(abs.(x .- x0))`.

Returns the index of the element in `x` closest to `x0` in absolute difference.

Throws an error if `x` is empty.
"""
function argmin_abs(x::AbstractVector{<:Real}, x0::Real)
    if isempty(x)
        throw(ArgumentError("collection must be non-empty"))
    end

    minval = NaN
    mini = 0

    @inbounds for i in eachindex(x)
        xi = x[i]
        d = abs(xi - x0)

        # skip if distance is NaN
        if isnan(d)
            continue
        end

        if isnan(minval) || d < minval
            minval = d
            mini = i
        end
    end

    if mini == 0
        throw(ArgumentError("all elements result in NaN distance from x0"))
    end

    return mini
end