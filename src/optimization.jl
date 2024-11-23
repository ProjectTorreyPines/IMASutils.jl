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

"""
    findall_interior_argmin(A::AbstractMatrix)

Returns vector of tuples with all (i,j) points where A[i,j] is
the lowest value of its eight neighbors A[i-1:i+1, j-1:j+1]
"""
findall_interior_argmin(A::AbstractMatrix) = _findall_interior_extrema(minimum, A)

"""
    findall_interior_argmax(A::AbstractMatrix)

Returns vector of tuples with all (i,j) points where A[i,j] is
the largest value of its eight neighbors A[i-1:i+1, j-1:j+1]
"""
findall_interior_argmax(A::AbstractMatrix) = _findall_interior_extrema(maximum, A)

function _findall_interior_extrema(fext::F, A::AbstractMatrix) where {F<:Function}
    #@assert fext in (minimum, maximum)
    Ni, Nj = size(A)
    @assert Ni > 2
    @assert Nj > 2
    pts = Tuple{Int, Int}[]
    for j in 2:(Nj-1)
        for i in 2:(Ni-1)
            if @inbounds A[i, j] == fext(@view(A[i-1:i+1, j-1:j+1]))
                push!(pts, (i, j))
            end
        end
    end
    return pts
end

"""
    find_interior_argmin(A::AbstractMatrix)

Returns tuple of (i,j) where A[i,j] is the lowest value of its eight neighbors A[i-1:i+1, j-1:j+1]

Errors if more than one argmin found
"""
find_interior_argmin(A::AbstractMatrix) = _find_interior_extrema(minimum, A)

"""
    find_interior_argmax(A::AbstractMatrix)

Returns tuple of (i,j) where A[i,j] is the largest value of its eight neighbors A[i-1:i+1, j-1:j+1]

Errors if more than one argmax found
"""
find_interior_argmax(A::AbstractMatrix) = _find_interior_extrema(maximum, A)

function _find_interior_extrema(fext::F, A::AbstractMatrix) where {F<:Function}
    @assert fext in (minimum, maximum)
    Ni, Nj = size(A)
    local im, ij
    found = false
    for j in 2:(Nj-1)
        for i in 2:(Ni-1)
            if @inbounds A[i, j] == fext(@view(A[i-1:i+1, j-1:j+1]))
                if !found
                    found = true
                    im, ij = i, j
                else
                    error("Multiple interior local $fext found at (i, j) = $((im, ij)) and $((i, j))")
                end
            end
        end
    end
    !found && error("No interior local $fext found")
    return im, ij
end