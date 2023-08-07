# # Containment/withinness

export within


"""
    within(geom1, geom)::Bool

Return `true` if the first geometry is completely within the second geometry.
The interiors of both geometries must intersect and, the interior and boundary of the primary (geometry a)
must not intersect the exterior of the secondary (geometry b).
`within` returns the exact opposite result of `contains`.

## Examples
```jldoctest setup=:(using GeometryOps, GeometryBasics)
julia> line = LineString([[1, 1], [1, 2], [1, 3], [1, 4]])
LineString(Array{Float64,1}[[1.0, 1.0], [1.0, 2.0], [1.0, 3.0], [1.0, 4.0]])

julia> point = Point([1, 2])
Point([1.0, 2.0])

julia> within(point, line)
true
```
"""
within(g1, g2)::Bool = within(trait(g1), g1, trait(g2), g2)::Bool
within(::GI.FeatureTrait, g1, ::Any, g2)::Bool = within(GI.geometry(g1), g2)
within(::Any, g1, t2::GI.FeatureTrait, g2)::Bool = within(g1, geometry(g2))
within(::GI.PointTrait, g1, ::GI.LineStringTrait, g2)::Bool = point_on_line(g1, g2; ignore_end_vertices=true)
within(::GI.PointTrait, g1, ::GI.PolygonTrait, g2)::Bool = point_in_polygon(g1, g2; ignore_boundary=true)
within(::GI.LineStringTrait, g1, ::GI.PolygonTrait, g2)::Bool = line_in_polygon(g1, g2)
within(::GI.LineStringTrait, g1, ::GI.LineStringTrait, g2)::Bool = line_on_line(g1, g2)
within(::GI.PolygonTrait, g1, ::GI.PolygonTrait, g2)::Bool = polygon_in_polygon(g1, g2)
