# # Disjointness checks

"""
    disjoint(geom1, geom2)::Bool

Return `true` if the intersection of the two geometries is an empty set.

# Examples

```jldoctest
import GeometryOps as GO, GeoInterface as GI

poly = GI.Polygon([[(-1, 2), (3, 2), (3, 3), (-1, 3), (-1, 2)]])
point = (1, 1)
GO.disjoint(poly, point)

# output
true
```
"""
disjoint(g1, g2)::Bool = disjoint(trait(g1), g1, trait(g2), g2)
disjoint(::FeatureTrait, g1, ::Any, g2)::Bool = disjoint(GI.geometry(g1), g2)
disjoint(::Any, g1, t2::FeatureTrait, g2)::Bool = disjoint(g1, geometry(g2))
disjoint(::PointTrait, g1, ::PointTrait, g2)::Bool = !point_equals_point(g1, g2)
disjoint(::PointTrait, g1, ::LineStringTrait, g2)::Bool = !point_on_line(g1, g2)
disjoint(::PointTrait, g1, ::PolygonTrait, g2)::Bool = !point_in_polygon(g1, g2)
disjoint(::LineStringTrait, g1, ::PointTrait, g2)::Bool = !point_on_line(g2, g1)
disjoint(::LineStringTrait, g1, ::LineStringTrait, g2)::Bool = !line_on_line(g1, g2)
disjoint(::LineStringTrait, g1, ::PolygonTrait, g2)::Bool = !line_in_polygon(g2, g1)
disjoint(::PolygonTrait, g1, ::PointTrait, g2)::Bool = !point_in_polygon(g2, g1)
disjoint(::PolygonTrait, g1, ::LineStringTrait, g2)::Bool = !line_in_polygon(g2, g1)
disjoint(::PolygonTrait, g1, ::PolygonTrait, g2)::Bool = polygon_disjoint(g2, g1)

function polygon_disjoint(poly1, poly2)
    for point in GI.getpoint(poly1)
        point_in_polygon(point, poly2) && return false
    end
    for point in GI.getpoint(poly2)
        point_in_polygon(point, poly1) && return false
    end
    return !intersects(poly1, poly2)
end
