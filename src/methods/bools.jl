# # Boolean conditions

export isclockwise, isconcave
export point_on_line, point_in_polygon, point_in_ring
export line_on_line, line_in_polygon, polygon_in_polygon

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
