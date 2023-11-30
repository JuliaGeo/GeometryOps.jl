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

# Lines in/on geometries

# Lines in/on lines
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

# Lines in/on rings
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

# Lines in/on polygons
within(
    ::GI.LineStringTrait, line,
    ::GI.PolygonTrait, polygon,
) = _line_polygon_process(
    line, polygon;
    process = within_process,
    exclude_boundaries = false,
    close = false,
)

# Rings in/on geometries

# Rings in/on lines
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

# Rings in/on rings
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

# Rings in/on polygons
within(
    ::GI.LinearRingTrait, line,
    ::GI.PolygonTrait, polygon,
) = _line_polygon_process(
    line, polygon;
    process = within_process,
    exclude_boundaries = false,
    close = true,
)

# Polygons within polygons
function within(
    ::GI.PolygonTrait, poly1,
    ::GI.PolygonTrait, poly2;
)
    if _line_polygon_process(
        GI.getexterior(poly1), poly2;
        process = within_process,
        exclude_boundaries = false,
        close = true,
        line_is_poly_ring = true
    )
        for hole in GI.gethole(poly2)
            if _line_polygon_process(
                hole, poly1;
                process = within_process,
                exclude_boundaries = false,
                close = true,
                line_is_poly_ring = true
            )
                return false
            end
        end
        return true
    end
    return false
end


# Everything not specified
# TODO: Add multipolygons
within(::GI.AbstractTrait, g1, ::GI.AbstractCurveTrait, g2)::Bool = false


