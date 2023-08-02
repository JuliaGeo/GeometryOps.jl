# # Containment

export contains

"""
    contains(ft1::AbstractGeometry, ft2::AbstractGeometry)::Bool

Return true if the second geometry is completely contained by the first geometry.
The interiors of both geometries must intersect and, the interior and boundary of the secondary (geometry b)
must not intersect the exterior of the primary (geometry a).
`contains` returns the exact opposite result of `within`.

## Examples

```jldoctest
line = GI.LineString([[1, 1], [1, 2], [1, 3], [1, 4]])
point = Point([1, 2])
contains(line, point)
# output
true
```
"""
contains(g1, g2)::Bool = within(g2, g1)


# This currently works for point-in-linestring or point-in-polygon.

# More GeometryBasics code

# _cross(p1, p2, p3) = (GI.x(p1) - GI.x(p3)) * (GI.y(p2) - GI.y(p3)) - (GI.x(p2) - GI.x(p3)) * (GI.y(p1) - GI.y(p3))

# """
#     contains(pointlist, point)::Bool

# Returns `true` if `point` is contained in `pointlist` (geometrically, not as a set)
# ,and  `false` otherwise.
# """
# contains(pointlist, point) = contains(GI.trait(pointlist), GI.trait(point), pointlist, point)

# Implementation of a point-in-polygon algorithm
# from Luxor.jl.  This is the Hormann-Agathos (2001) algorithm.

# For the source, see [the code from Luxor.jl](https://github.com/JuliaGraphics/Luxor.jl/blob/66d60fb51f6b1bb38690fe8dcc6c0084eeb80710/src/polygons.jl#L190-L229).

#function contains(::Union{GI.LineStringTrait, GI.LinearRingTrait}, ::GI.PointTrait, pointlist, point)
#    n = GI.npoint(pointlist)
#    c = false
#    q1 = GI.getpoint(pointlist, 1)
#    q2 = GI.getpoint(pointlist, 1)
#    @inbounds for (counter, current_point) in enumerate(Iterators.drop(GI.getpoint(pointlist), 1))
#        q1 = q2
#        ## if reached last point, set "next point" to first point.
#        ##
#        if counter == (n-1)
#            q2 = GI.getpoint(pointlist, 1)
#        else
#            q2 = current_point
#        end
#        if GI.x(q1) == GI.x(point) && GI.x(q1) == GI.y(point)
#            ## allowonedge || error("isinside(): VertexException a")
#            continue
#        end
#        if GI.y(q2) == GI.y(point)
#            if GI.x(q2) == GI.x(point)
#                ## allowonedge || error("isinside(): VertexException b")
#                continue
#            elseif (GI.y(q1) == GI.y(point)) && ((GI.x(q2) > GI.x(point)) == (GI.x(q1) < GI.x(point)))
#                ## allowonedge || error("isinside(): EdgeException")
#                continue
#            end
#        end
#        if (GI.y(q1) < GI.y(point)) != (GI.y(q2) < GI.y(point)) # crossing
#            if GI.x(q1) >= GI.x(point)
#                if GI.x(q2) > GI.x(point)
#                    c = !c
#                elseif ((_cross(q1, q2, point) > 0) == (GI.y(q2) > GI.y(q1)))
#                    c = !c
#                end
#            elseif GI.x(q2) > GI.x(point)
#                if ((_cross(q1, q2, point) > 0) == (GI.y(q2) > GI.y(q1)))
#                    c = !c
#                end
#            end
#        end
#    end
#    return c

#end
#function contains(poly::Polygon{2, T1}, point::Point{2, T2}) where {T1, T2}
#    c = contains(poly.exterior, point)
#    for interior in poly.interiors
#        if contains(interior, point)
#            return false
#        end
#    end
#    return c
#end

## TODOs: implement contains for mesh, simplex, and 3d objects (eg rect, triangle, etc.)

#contains(mp::MultiPolygon{2, T1}, point::Point{2, T2}) where {T1, T2} = any((contains(poly, point) for poly in mp.polygons))
