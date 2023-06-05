# # Disjointness checks

"""
    disjoint(geom1, geom2)::Bool

Return `true` if the intersection of the two geometries is an empty set.

# Examples
```jldoctest
julia> poly = Polygon([[[-1, 2], [3, 2], [3, 3], [-1, 3], [-1, 2]]])
Polygon(Array{Array{Float64,1},1}[[[-1.0, 2.0], [3.0, 2.0], [3.0, 3.0], [-1.0, 3.0], [-1.0, 2.0]]])

julia> point = Point([1, 1])
Point([1.0, 1.0])

julia> disjoint(poly, point)
true
```
"""
disjoint(t1::FeatureTrait, g1, t2, g2)::Bool = disjoint(GI.geometry(g1), g2)
disjoint(t1, g1, t2::FeatureTrait, g2)::Bool = disjoint(g1, geometry(g2))
disjoint(t1::PointTrait, g1, t2::PointTrait, g2)::Bool = !point_equals_point(g1, g2)
disjoint(t1::PointTrait, g1, t2::LineStringTrait, g2)::Bool = !point_on_line(g1, g2)
disjoint(t1::PointTrait, g1, t2::PolygonTrait, g2)::Bool = !point_in_polygon(g1, g2)
disjoint(t1::LineStringTrait, g1, t2::PointTrait, g2)::Bool = !point_on_line(g2, g1)
disjoint(t1::LineStringTrait, g1, t2::LineStringTrait, g2)::Bool = !line_on_line(g1, g2)
disjoint(t1::LineStringTrait, g1, t2::PolygonTrait, g2)::Bool = !line_in_polygon(g2, g1)
disjoint(t1::PolygonTrait, g1, t2::PointTrait, g2)::Bool = !point_in_polygon(g2, g1)
disjoint(t1::PolygonTrait, g1, t2::LineStringTrait, g2)::Bool = !line_in_polygon(g2, g1)
disjoint(t1::PolygonTrait, g1, t2::PolygonTrait, g2)::Bool = !poly_in_poly(g2, g1)
