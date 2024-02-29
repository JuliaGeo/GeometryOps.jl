# # Orientation

export isclockwise, isconcave


# #### `isclockwise`
# The orientation of a geometry is whether it runs clockwise or counter-clockwise.

# This is defined for linestrings, linear rings, or vectors of points.  

# #### `isconcave`
# A polygon is concave if it has at least one interior angle greater than 180 degrees, 
# meaning that the interior of the polygon is not a convex set.

# These are all adapted from Turf.jl.

# The may not necessarily be what want in the end but work for now!

"""
    isclockwise(line::Union{LineString, Vector{Position}})::Bool

Take a ring and return true or false whether or not the ring is clockwise or
counter-clockwise.

## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

ring = GI.LinearRing([(0, 0), (1, 1), (1, 0), (0, 0)])
GO.isclockwise(ring)

# output
true
```
"""
isclockwise(geom)::Bool = isclockwise(GI.trait(geom), geom)

function isclockwise(::AbstractCurveTrait, line)::Bool
    sum = 0.0
    prev = GI.getpoint(line, 1)
    for p in GI.getpoint(line)
        # sum will be zero for the first point as x is subtracted from itself
        sum += (GI.x(p) - GI.x(prev)) * (GI.y(p) + GI.y(prev))
        prev = p
    end

    return sum > 0.0
end

"""
    isconcave(poly::Polygon)::Bool

Take a polygon and return true or false as to whether it is concave or not.

## Examples
```jldoctest
import GeoInterface as GI, GeometryOps as GO

poly = GI.Polygon([[(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]])
GO.isconcave(poly)

# output
false
```
"""
function isconcave(poly)::Bool
    sign = false

    exterior = GI.getexterior(poly)

    # FIXME handle not closed polygons
    GI.npoint(exterior) <= 4 && return false
    n = GI.npoint(exterior) - 1

    for i in 1:n
        j = ((i + 1) % n) === 0 ? 1 : (i + 1) % n
        m = ((i + 2) % n) === 0 ? 1 : (i + 2) % n

        pti = GI.getpoint(exterior, i)
        ptj = GI.getpoint(exterior, j)
        ptm = GI.getpoint(exterior, m)

        dx1 = GI.x(ptm) - GI.x(ptj)
        dy1 = GI.y(ptm) - GI.y(ptj)
        dx2 = GI.x(pti) - GI.x(ptj)
        dy2 = GI.y(pti) - GI.y(ptj)

        cross = (dx1 * dy2) - (dy1 * dx2)

        if i === 0
            sign = cross > 0
        elseif sign !== (cross > 0)
            return true
        end
    end

    return false
end

#=
This is commented out.
```julia
"""
    isparallel(line1::LineString, line2::LineString)::Bool

Return `true` if each segment of `line1` is parallel to the correspondent segment of `line2`

## Examples
```julia
import GeoInterface as GI, GeometryOps as GO
julia> line1 = GI.LineString([(9.170356, 45.477985), (9.164434, 45.482551), (9.166644, 45.484003)])
GeoInterface.Wrappers.LineString{false, false, Vector{Tuple{Float64, Float64}}, Nothing, Nothing}([(9.170356, 45.477985), (9.164434, 45.482551), (9.166644, 45.484003)], nothing, nothing)

julia> line2 = GI.LineString([(9.169356, 45.477985), (9.163434, 45.482551), (9.165644, 45.484003)])
GeoInterface.Wrappers.LineString{false, false, Vector{Tuple{Float64, Float64}}, Nothing, Nothing}([(9.169356, 45.477985), (9.163434, 45.482551), (9.165644, 45.484003)], nothing, nothing)

julia> 
GO.isparallel(line1, line2)
true
```
"""
function isparallel(line1, line2)::Bool
    seg1 = linesegment(line1)
    seg2 = linesegment(line2)

    for i in eachindex(seg1)
        coors2 = nothing
        coors1 = seg1[i]
        coors2 = seg2[i]
        _isparallel(coors1, coors2) == false && return false
    end
    return true
end

@inline function _isparallel(p1, p2)
    slope1 = bearing_to_azimuth(rhumb_bearing(GI.x(p1), GI.x(p2)))
    slope2 = bearing_to_azimuth(rhumb_bearing(GI.y(p1), GI.y(p2)))

    return slope1 === slope2
end
```
=#

# This is actual code:

_isparallel(((ax, ay), (bx, by)), ((cx, cy), (dx, dy))) = 
    _isparallel(bx - ax, by - ay, dx - cx, dy - cy)

_isparallel(Δx1, Δy1, Δx2, Δy2) = (Δx1 * Δy2 == Δy1 * Δx2)  


