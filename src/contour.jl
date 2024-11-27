
# This keeps track of which edges are connected for a given setup of vertices
# Binary on the right denotes which corners are greater than the level value
#    ordered CCW from lower-left
# The result shows which edges are connected, in CCW fashion assuming values
#    inside the contour are greater than the level value
# Note that 1010 and 0101 are ambiguous saddle points with two connections
#   default assumes the value at the cell center is less than the level value
#   alternate assumes the value at the cell center is greater than the level value
const EDGE_TABLE = @SVector[
    SVector{2, Int}[],                                 # 0000
    [SVector{2, Int}(1, 4)],                           # 1000
    [SVector{2, Int}(2, 1)],                           # 0100
    [SVector{2, Int}(2, 4)],                           # 1100
    [SVector{2, Int}(3, 2)],                           # 0010
    [SVector{2, Int}(1, 4), SVector{2, Int}(3, 2)],    # 1010 default
    [SVector{2, Int}(3, 1)],                           # 0110
    [SVector{2, Int}(3, 4)],                           # 1110
    [SVector{2, Int}(4, 3)],                           # 0001
    [SVector{2, Int}(1, 3)],                           # 1001
    [SVector{2, Int}(4, 3), SVector{2, Int}(2, 1)],    # 0101 default
    [SVector{2, Int}(2, 3)],                           # 1101
    [SVector{2, Int}(4, 2)],                           # 0011
    [SVector{2, Int}(1, 2)],                           # 1011
    [SVector{2, Int}(4, 1)],                           # 0111
    SVector{2, Int}[],                                 # 1111
    [SVector{2, Int}(3, 4), SVector{2, Int}(1, 2)],    # 1010 alternate
    [SVector{2, Int}(4, 1), SVector{2, Int}(2, 3)]     # 0101 alternate
]

# Helper function to interpolate between two points
@inline function interpolate(p1::SVector{2, T},
                             p2::SVector{2, T},
                             v1::T, v2::T, c::T) where {T<:Real}
    t = (c - v1) / (v2 - v1)
    return p1 + t * (p2 - p1)
end

# Get value at vertices for cell with (i, j) in lower left
@inline function get_vertices(values::AbstractMatrix{T}, i::Int, j::Int) where {T<:Real}
    return values[i, j], values[i+1, j], values[i+1, j+1], values[i, j+1]
end

# Get value at vertices for cell with (i, j) in lower left, assuming inbounds
@inline function get_vertices_inbounds(values::AbstractMatrix{T}, i::Int, j::Int) where {T<:Real}
    @inbounds v1, v2, v3, v4 = values[i, j], values[i+1, j], values[i+1, j+1], values[i, j+1]
    return v1, v2, v3, v4
end

# Determine which edge case these vertices and level value correspond to
@inline function get_case(v1::T, v2::T, v3::T, v4::T, level::T) where {T<:Real}
    case_index = (v1 > level ? 1 : 0) + (v2 > level ? 2 : 0) +
                 (v3 > level ? 4 : 0) + (v4 > level ? 8 : 0) + 1

    if case_index in (6, 11)
        vm = 0.25 * (v1 + v2 + v3 + v4)
        (vm > level) && (case_index = (case_index == 6) ? 17 : 18)
    end
    return case_index
end

# Get the edges for a given case with these vertices and level value
@inline function get_edges(v1::T, v2::T, v3::T, v4::T, level::T) where {T<:Real}
    case_index = get_case(v1, v2, v3, v4, level)
    return get_edges(case_index)
end
@inline get_edges(case_index::Int) = EDGE_TABLE[case_index]


# Compute the intersecting point of the contour with this edge
@inline function get_point(edge::Int, x_coords::AbstractVector{T}, y_coords::AbstractVector{T},
                           i::Int, j::Int, v1::T, v2::T, v3::T, v4::T, level::T) where {T<:Real}
    if edge == 1
        p1 = SVector(x_coords[i], y_coords[j])
        p2 = SVector(x_coords[i+1], y_coords[j])
        v_start, v_end = v1, v2
    elseif edge == 2
        p1 = SVector(x_coords[i+1], y_coords[j])
        p2 = SVector(x_coords[i+1], y_coords[j+1])
        v_start, v_end = v2, v3
    elseif edge == 3
        p1 = SVector(x_coords[i+1], y_coords[j+1])
        p2 = SVector(x_coords[i], y_coords[j+1])
        v_start, v_end = v3, v4
    elseif edge == 4
        p1 = SVector(x_coords[i], y_coords[j+1])
        p2 = SVector(x_coords[i], y_coords[j])
        v_start, v_end = v4, v1
    end

    return interpolate(p1, p2, v_start, v_end, level)
end

# return the two interecting points for the contour in this cell
@inline function get_segment(edges, x_coords, y_coords, i, j, v1, v2, v3, v4, c)
    point1 = get_point(edges[1], x_coords, y_coords, i, j, v1, v2, v3, v4, c)
    point2 = get_point(edges[2], x_coords, y_coords, i, j, v1, v2, v3, v4, c)
    return (point1, point2)
end

# Find the next cell in the forward direction
function next_forward(i, j, edges)
    fedge = edges[2]
    if fedge == 1
        return i, j-1
    elseif fedge == 2
        return i+1, j
    elseif fedge == 3
        return i, j+1
    elseif fedge == 4
        return i-1, j
    end
end

# Find the next cell in the backward direction
function next_backward(i, j, edges)
    bedge = edges[1]
    if bedge == 1
        return i, j-1
    elseif bedge == 2
        return i+1, j
    elseif bedge == 3
        return i, j+1
    elseif bedge == 4
        return i-1, j
    end
end


"""
    contour_cache(values::Matrix{T}; aggression_level::Int = 2) where {T<:Real}

For a given matrix `values`, create cache for x_contour and y_contour
`aggression_level` determines how large to make the cache:

1. assume contour no bigger than an inscribed ellipse

2. assume contour no bigger than the other boundary

3. assume contour could go through every single cell once
"""
function contour_cache(values::Matrix{T}; aggression_level::Int = 2) where {T<:Real}
    @assert aggression_level in (1, 2, 3)
    nx, ny = size(values)
    if aggression_level == 1
        # assume inscribed ellipse length
        Ncache = ceil(Int, π * sqrt(2.0 * (nx^2 + ny^2)))
    elseif aggression_level == 2
        # assume rectangle
        Ncache = 2 * (nx + ny)
    else #aggression_level == 3
        # assume area filling
        Ncache = nx * ny
    end
    return Vector{T}(undef, Ncache), Vector{T}(undef, Ncache)
end

"""
    contour_from_midplane(values::Matrix{T},
                          x_coords::AbstractVector{T},
                          y_coords::AbstractVector{T},
                          level::T, xaxis::T, yaxis::T, vaxis::T) where {T<:Real}

Find a contour of (x_coords, y_coords, values) at value=level that crosses y=yaxis at the smallest x > xaxis

This will correspond to a closed surface around the axis if it exists, otherwise it will give upto one of possibly multiple open surfaces

Returns vectors for x_contour and y_contour
"""
function contour_from_midplane(values::Matrix{T},
                          x_coords::AbstractVector{T},
                          y_coords::AbstractVector{T},
                          level::T, xaxis::T, yaxis::T, vaxis::T) where {T<:Real}
    x_cache, y_cache = contour_cache(values)
    x_contour, y_contour = contour_from_midplane!(x_cache, y_cache, values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
    return collect(x_contour), collect(y_contour)
end

# give an empty view of the cache vectors if no contour found
@inline function empty_contour(x_cache::Vector{T}, y_cache::Vector{T}) where {T<:Real}
    x_empty = @view x_cache[1:0]
    y_empty = @view y_cache[1:0]
    return x_empty, y_empty
end

# compare (x, y) of each point within atol and rtol
@inline function approx_point(a::SVector{2, T}, b::SVector{2, T}; atol=eps(T), rtol=sqrt(eps(T))) where {T<:Real}
    return isapprox(a[1], b[1]; atol, rtol) && isapprox(a[2], b[2]; atol, rtol)
end

"""
    contour_from_midplane!(x_cache::Vector{T}, y_cache::Vector{T},
                           values::Matrix{T},
                           x_coords::AbstractVector{T},
                           y_coords::AbstractVector{T},
                           level::T, xaxis::T, yaxis::T, vaxis::T;
                           atol::T=eps(T), rtol::T=sqrt(eps(T))) where {T<:Real}

Find a contour of (x_coords, y_coords, values) at value=level that crosses y=yaxis at the smallest x > xaxis

This will correspond to a closed surface around the axis if it exists, otherwise it will give upto one of possibly multiple open surfaces

The contour is computed in-place using x_cache and y_cache, and returned as `view`s of those caches to avoid allocations

atol and rtol are used for checking if the contour closes on itself in every cell
"""
function contour_from_midplane!(x_cache::Vector{T}, y_cache::Vector{T},
                                values::Matrix{T},
                                x_coords::AbstractVector{T},
                                y_coords::AbstractVector{T},
                                level::T, xaxis::T, yaxis::T, vaxis::T;
                                atol::T=eps(T), rtol::T=sqrt(eps(T))) where {T<:Real}

    x_cache .= NaN
    y_cache .= NaN
    nx, ny = size(values)

    direction = level < vaxis ? :decreasing : :increasing

    ia = searchsortedlast(x_coords, xaxis) #argmin(abs(x - xaxis) for x in x_coords)
    ja = searchsortedlast(y_coords, yaxis) #argmin(abs(y - yaxis) for y in y_coords)
    inv_dc = 1.0 / (level - vaxis)

    # check normalized value around axis
    # greater than 1 is outside the current contour
    v1, v2, v3, v4 = get_vertices(values, ia, ja)
    v1norm = (v1 - vaxis) * inv_dc
    v2norm = (v2 - vaxis) * inv_dc
    v3norm = (v3 - vaxis) * inv_dc
    v4norm = (v4 - vaxis) * inv_dc
    if v1norm > 1 && v2norm > 1 && v3norm > 1 && v4norm > 1
        # contour is inside one grid cell
        return empty_contour(x_cache, y_cache)
    elseif v1norm >=1 || v2norm >= 1 || v3norm >= 1 || v4norm >= 1
        # contour goes through first cell
        istart = ia
    else
        # traverse outward in x until we find a value outside the contour
        # then the contour must go through that cell
        istart = 0
        for i in (ia+1):nx
            vnorm = (values[i, ja] - vaxis) * inv_dc
            if vnorm >= 1.0
                # if =1, contour crossing happens in next cell
                # if >1, then it was in the last cell
                istart = (vnorm == 1.0) ? i : i-1
                break
            end
        end
        istart == 0 && return empty_contour(x_cache, y_cache)
    end

    jstart = ja
    status = 0

    # Store the first two points in the forward direction
    v1, v2, v3, v4 = get_vertices_inbounds(values, istart, jstart)
    first_edges = get_edges(v1, v2, v3, v4, level)[1] # either it's open or the first one is the closed one
    first_point, second_point = get_segment(first_edges, x_coords, y_coords, istart, jstart, v1, v2, v3, v4, level)
    nforward = 1
    x_cache[nforward], y_cache[nforward] = first_point
    nforward = 2
    x_cache[nforward], y_cache[nforward] = second_point
    last_point = second_point

    # Traverse contour in forward direction
    i, j = next_forward(istart, jstart, first_edges)
    while status == 0
        if i < 1 || j < 1 || i >= nx || j >= ny
            # We've reached the boundary, so almost certainly an open contour
            status = 1
            break
        end

        v1, v2, v3, v4 = get_vertices_inbounds(values, i, j)
        edges = get_edges(v1, v2, v3, v4, level)

        point1, point2 = get_segment(edges[1], x_coords, y_coords, i, j, v1, v2, v3, v4, level)
        if approx_point(point1, last_point; atol, rtol)
            # First (usually only) segment in cell connects to previous segment
            nforward += 1
            x_cache[nforward], y_cache[nforward] = point2
            last_point = point2
            i, j = next_forward(i, j, edges[1])
        elseif length(edges) == 2
            # There's a saddle point, so check the second segment in the cell
            point1, point2 = get_segment(edges[2], x_coords, y_coords, i, j, v1, v2, v3, v4, level)
            if approx_point(point1, last_point; atol, rtol)
                # this connects
                nforward += 1
                x_cache[nforward], y_cache[nforward] = point2
                last_point = point2
                i, j = next_forward(i, j, edges[2])
            else
                error("Failed to connect to last point")
            end
        else
            error("Failed to connect to last point")
        end

        if approx_point(point2, first_point; atol, rtol)
            # This is a closed contour
            x_cache[nforward], y_cache[nforward] = x_cache[1], y_cache[1]
            status = 2
        end
    end

    # If we're on an open contour, so we want to trace it out in the reverse direction
    # the first two points we're already recorded last time
    if status == 1
        i, j = next_backward(istart, jstart, first_edges)
        last_forward = SVector(x_cache[nforward], y_cache[nforward])
    end
    nbackward = 0
    last_point = first_point
    while status == 1

        if i < 1 || j < 1 || i >= nx || j >= ny
            # We've reached other end of an open contour
            status = 3
            break
        end

        # Traverse backward
        v1, v2, v3, v4 = get_vertices_inbounds(values, i, j)
        edges = get_edges(v1, v2, v3, v4, level)

        # points in reverse
        point2, point1 = get_segment(edges[1], x_coords, y_coords, i, j, v1, v2, v3, v4, level)
        if approx_point(point1, last_point; atol, rtol)
            # First (usually only) segment in cell connects to previous segment
            nbackward += 1
            x_cache[nforward + nbackward], y_cache[nforward + nbackward] = point2
            last_point = point2
            i, j = next_backward(i, j, edges[1])
        elseif length(edges) == 2
            # There's a saddle point, so check the second segment in the cell

            # points in reverse
            point2, point1 = get_segment(edges[2], x_coords, y_coords, i, j, v1, v2, v3, v4, level)
            if approx_point(point1, last_point; atol, rtol)
                # this connects
                nbackward += 1
                x_cache[nforward + nbackward], y_cache[nforward + nbackward] = point2
                last_point = point2
                i, j = next_backward(i, j, edges[2])
            else
                error("Failed to connect to last point")
            end
        else
            error("Failed to connect to last point")
        end

        if approx_point(point2, last_forward; atol, rtol)
            # We have a closed contour, so it must be tangent to the boundary
            x_cache[nforward + nbackward], y_cache[nforward + nbackward] = x_cache[nforward], y_cache[nforward]
            status = 4
        end
    end

    # Reverse forward points nforward and nforward+1 are adjacent
    x_forward = @view x_cache[1:nforward]
    y_forward = @view y_cache[1:nforward]
    reverse!(x_forward)
    reverse!(y_forward)

    x_contour = @view x_cache[1:(nforward + nbackward)]
    y_contour = @view y_cache[1:(nforward + nbackward)]

    if direction === :decreasing
        # all the points will be clockwise, so reverse them
        reverse!(x_contour)
        reverse!(y_contour)
    end # if :increasing, they're already CCW

    return x_contour, y_contour

end