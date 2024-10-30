"""
    mirror_bound(x::T, l::T, u::T) where {T<:Real}

Return tuple with value of x bounded between l and u
The bounding is done by mirroring the value at the bound limits.
"""
function mirror_bound(x::T, l::T, u::T) where {T<:Real}
    d = (u - l) / 2.0
    c = (u + l) / 2.0
    x0 = (x .- c) / d
    while abs(x0) > 1.0
        if x0 < 1.0
            x0 = -2.0 - x0
        else
            x0 = 2.0 - x0
        end
    end
    return x0 * d + c
end