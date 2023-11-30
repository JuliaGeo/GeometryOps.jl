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
# Syntactic sugar
disjoint(g1, g2)::Bool = disjoint(trait(g1), g1, trait(g2), g2)
disjoint(::FeatureTrait, g1, ::Any, g2)::Bool = disjoint(GI.geometry(g1), g2)
disjoint(::Any, g1, t2::FeatureTrait, g2)::Bool = disjoint(g1, geometry(g2))
# Point disjoint geometries

# Point disjoint from point
disjoint(
    ::GI.PointTrait, g1,
    ::GI.PointTrait, g2,
) = !equals(g1, g2)

# Point in from line string
disjoint(
    ::GI.PointTrait, g1,
    ::GI.LineStringTrait, g2,
) = _point_curve_process(
    g1, g2;
    process = disjoint_process,
    exclude_boundaries = false,
    repeated_last_coord = false,
)

# Point disjoint from curve
disjoint(
    ::GI.PointTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _point_curve_process(
    g1, g2;
    process = disjoint_process,
    exclude_boundaries = false,
    repeated_last_coord = true,
)

# Point disjoint from polygon
disjoint(
    ::GI.PointTrait, g1,
    ::GI.PolygonTrait, g2,
) = _point_polygon_process(
    g1, g2;
    process = disjoint_process, exclude_boundaries = false,
)

# Geometry with point
disjoint(
    trait1::GI.AbstractTrait, g1,
    trait2::GI.PointTrait, g2,
) = disjoint(trait2, g2, trait1, g1)

# Lines disjoint from geometries

# Lines disjoint from lines
disjoint(
    ::GI.LineStringTrait, line,
    ::GI.LineStringTrait, curve,
) = _line_curve_process(
    line, curve;
    process = disjoint_process,
    exclude_boundaries = false,
    closed_line = false,
    closed_curve = false,
)

# Lines disjoint from rings
disjoint(
    ::GI.LineStringTrait, line,
    ::GI.LinearRingTrait, ring,
) = _line_curve_process(
    line, ring;
    process = disjoint_process,
    exclude_boundaries = false,
    closed_line = false,
    closed_curve = true,
)

# Lines disjoint from polygons
disjoint(
    ::GI.LineStringTrait, line,
    ::GI.PolygonTrait, polygon,
) = _line_polygon_process(
    line, polygon;
    process = disjoint_process,
    exclude_boundaries = false,
    close = false,
)

# Rings disjoint from geometries

# Rings disjoint from lines
disjoint(
    ::GI.LinearRingTrait, line,
    ::GI.LineStringTrait, curve,
) = _line_curve_process(
    line, curve;
    process = disjoint_process,
    exclude_boundaries = false,
    closed_line = true,
    closed_curve = false,
)

# Rings disjoint from rings
disjoint(
    ::GI.LinearRingTrait, line,
    ::GI.LinearRingTrait, ring,
) = _line_curve_process(
    line, ring;
    process = disjoint_process,
    exclude_boundaries = false,
    closed_line = true,
    closed_curve = true,
)

# Rings disjoint from polygons
disjoint(
    ::GI.LinearRingTrait, line,
    ::GI.PolygonTrait, polygon,
) = _line_polygon_process(
    line, polygon;
    process = disjoint_process,
    exclude_boundaries = false,
    close = true,
)

# Polygons disjoint from polygons
function disjoint(
    ::GI.PolygonTrait, poly1,
    ::GI.PolygonTrait, poly2;
)
    if disjoint(GI.getexterior(poly1), poly2)
        return true
    else
        for hole in GI.gethole(poly1)
            if within(poly2, hole)
                return true
            end
        end
    end
    return false
end

