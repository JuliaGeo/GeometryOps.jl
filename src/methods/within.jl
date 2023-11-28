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
import GeometryOps as GO, GeoInterface as GI

line = GI.LineString([(1, 1), (1, 2), (1, 3), (1, 4)])
point = (1, 2)
GO.within(point, line)

# output
true
```
"""
# Syntactic sugar
within(g1, g2)::Bool = within(trait(g1), g1, trait(g2), g2)::Bool
within(::GI.FeatureTrait, g1, ::Any, g2)::Bool = within(GI.geometry(g1), g2)
within(::Any, g1, t2::GI.FeatureTrait, g2)::Bool = within(g1, GI.geometry(g2))

# Points within geometries
# Point in point
within(
    ::GI.PointTrait, g1,
    ::GI.PointTrait, g2,
) = equals(g1, g2)
# Point in (on) line string
within(
    ::GI.PointTrait, g1,
    ::GI.LineStringTrait, g2,
) = _point_curve_process(
    g1, g2;
    process = within_process,
    exclude_boundaries = true,
    repeated_last_coord = false,
)
# Point in (on) curve
within(
    ::GI.PointTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _point_curve_process(
    g1, g2;
    process = within_process,
    exclude_boundaries = true,
    repeated_last_coord = true,
)

# Point in polygon
within(
    ::GI.PointTrait, g1,
    ::GI.PolygonTrait, g2,
) = _point_polygon_process(
    g1, g2;
    process = within_process, exclude_boundaries = true,
)

# Lines in geometries
within(
    ::GI.LineStringTrait, line,
    ::GI.LineStringTrait, curve,
) = _line_curve_process(
    line, curve;
    process = within_process,
    exclude_boundaries = false,
    closed_line = false,
    closed_curve = false,
)

within(
    ::GI.LineStringTrait, line,
    ::GI.LinearRingTrait, ring,
) = _line_curve_process(
    line, ring;
    process = within_process,
    exclude_boundaries = false,
    closed_line = false,
    closed_curve = true,
)

within(
    ::GI.LineStringTrait, line,
    ::GI.PolygonTrait, polygon,
) = _line_polygon_process(
    line, polygon;
    process = within_process,
    exclude_boundaries = false,
    close = false,
)

# Rings in geometries
within(
    ::GI.LinearRingTrait, line,
    ::GI.LineStringTrait, curve,
) = _line_curve_process(
    line, curve;
    process = within_process,
    exclude_boundaries = false,
    closed_line = true,
    closed_curve = false,
)

within(
    ::GI.LinearRingTrait, line,
    ::GI.LinearRingTrait, ring,
) = _line_curve_process(
    line, ring;
    process = within_process,
    exclude_boundaries = false,
    closed_line = true,
    closed_curve = true,
)

within(
    ::GI.LinearRingTrait, line,
    ::GI.PolygonTrait, polygon,
) = _line_polygon_process(
    line, polygon;
    process = within_process,
    exclude_boundaries = true,
    close = true,
)

# Polygons within geometries
within(::GI.PolygonTrait, g1, ::GI.PolygonTrait, g2)::Bool = polygon_in_polygon(g1, g2)

# Everything not specified
# TODO: Add multipolygons
within(::GI.AbstractTrait, g1, ::GI.AbstractCurveTrait, g2)::Bool = false


