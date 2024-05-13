# # Within

export within

#=
## What is within?

The within function checks if one geometry is inside another geometry. This
requires that the two interiors intersect and that the interior and
boundary of the first geometry is not in the exterior of the second geometry.

To provide an example, consider these two lines:
```@example within
import GeometryOps as GO 
import GeoInterface as GI
using Makie
using CairoMakie

l1 = GI.LineString([(0.0, 0.0), (1.0, 0.0), (0.0, 0.1)])
l2 = GI.LineString([(0.25, 0.0), (0.75, 0.0)])
f, a, p = lines(GI.getpoint(l1), color = :blue)
scatter!(GI.getpoint(l1), color = :blue)
lines!(GI.getpoint(l2), color = :orange)
scatter!(GI.getpoint(l2), color = :orange)
f
```
We can see that all of the points and edges of l2 are within l1, so l2 is
within l1, but l1 is not within l2
```@example within
GO.within(l1, l2)  # false
GO.within(l2, l1)  # true
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait.

Each of these calls a method in the geom_geom_processors file. The methods in
this file determine if the given geometries meet a set of criteria. For the
`within` function and arguments g1 and g2, this criteria is as follows:
    - points of g1 are allowed to be in the interior of g2 (either through
    overlap or crossing for lines)
    - points of g1 are allowed to be on the boundary of g2
    - points of g1 are not allowed to be in the exterior of g2
    - at least one point of g1 is required to be in the interior of g2
    - no points of g1 are required to be on the boundary of g2
    - no points of g1 are required to be in the exterior of g2

The code for the specific implementations is in the geom_geom_processors file.
=#

const WITHIN_POINT_ALLOWS = (in_allow = true, on_allow = false, out_allow = false)
const WITHIN_CURVE_ALLOWS = (over_allow = true, cross_allow = true, on_allow = true, out_allow = false)
const WITHIN_POLYGON_ALLOWS = (in_allow = true, on_allow = true, out_allow = false)
const WITHIN_REQUIRES = (in_require = true, on_require = false, out_require = false)
const WITHIN_EXACT = (exact = _False(),)

"""
    within(geom1, geom2)::Bool

Return `true` if the first geometry is completely within the second geometry.
The interiors of both geometries must intersect and the interior and boundary of
the primary geometry (geom1) must not intersect the exterior of the secondary
geometry (geom2).

Furthermore, `within` returns the exact opposite result of `contains`.

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
within(g1, g2) = _within(trait(g1), g1, trait(g2), g2)

# # Convert features to geometries
_within(::GI.FeatureTrait, g1, ::Any, g2) = within(GI.geometry(g1), g2)
_within(::Any, g1, t2::GI.FeatureTrait, g2) = within(g1, GI.geometry(g2))
_within(::FeatureTrait, g1, ::FeatureTrait, g2) = within(GI.geometry(g1), GI.geometry(g2))


# # Points within geometries

# Point is within another point if those points are equal.
_within(
    ::GI.PointTrait, g1,
    ::GI.PointTrait, g2,
) = equals(g1, g2)

#= Point is within a linestring if it is on a vertex or an edge of that line,
excluding the start and end vertex if the line is not closed. =#
_within(
    ::GI.PointTrait, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _point_curve_process(
    g1, g2;
    WITHIN_POINT_ALLOWS...,
    closed_curve = false,
)

# Point is within a linearring if it is on a vertex or an edge of that ring.
_within(
    ::GI.PointTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _point_curve_process(
    g1, g2;
    WITHIN_POINT_ALLOWS...,
    closed_curve = true,
)

#= Point is within a polygon if it is inside of that polygon, excluding edges,
vertices, and holes. =#
_within(
    ::GI.PointTrait, g1,
    ::GI.PolygonTrait, g2,
) = _point_polygon_process(
    g1, g2;
    WITHIN_POINT_ALLOWS...,
    WITHIN_EXACT...,
)

# No geometries other than points can be within points
_within(
    ::Union{GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    ::GI.PointTrait, g2,
) = false


# # Lines within geometries

#= Linestring is within another linestring if their interiors intersect and no
points of the first line are in the exterior of the second line. =#
_within(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _line_curve_process(
    g1, g2;
    WITHIN_CURVE_ALLOWS...,
    WITHIN_REQUIRES...,
    WITHIN_EXACT...,
    closed_line = false,
    closed_curve = false,
)

#= Linestring is within a linear ring if their interiors intersect and no points
of the line are in the exterior of the ring. =#
_within(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    WITHIN_CURVE_ALLOWS...,
    WITHIN_REQUIRES...,
    WITHIN_EXACT...,
    closed_line = false,
    closed_curve = true,
)

#= Linestring is within a polygon if their interiors intersect and no points of
the line are in the exterior of the polygon, although they can be on an edge. =#
_within(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    WITHIN_POLYGON_ALLOWS...,
    WITHIN_REQUIRES...,
    WITHIN_EXACT...,
    closed_line = false,
)


# # Rings covered by geometries

#= Linearring is within a linestring if their interiors intersect and no points
of the ring are in the exterior of the line. =#
_within(
    ::GI.LinearRingTrait, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _line_curve_process(
    g1, g2;
    WITHIN_CURVE_ALLOWS...,
    WITHIN_REQUIRES...,
    WITHIN_EXACT...,
    closed_line = true,
    closed_curve = false,
)

#= Linearring is within another linearring if their interiors intersect and no
points of the first ring are in the exterior of the second ring. =#
_within(
    ::GI.LinearRingTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    WITHIN_CURVE_ALLOWS...,
    WITHIN_REQUIRES...,
    WITHIN_EXACT...,
    closed_line = true,
    closed_curve = true,
)

#= Linearring is within a polygon if their interiors intersect and no points of
the ring are in the exterior of the polygon, although they can be on an edge. =#
_within(
    ::GI.LinearRingTrait, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    WITHIN_POLYGON_ALLOWS...,
    WITHIN_REQUIRES...,
    WITHIN_EXACT...,
    closed_line = true,
)


# # Polygons within geometries

#= Polygon is within another polygon if the interior of the first polygon 
intersects with the interior of the second and no points of the first polygon
are outside of the second polygon. =#
_within(
    ::GI.PolygonTrait, g1,
    ::GI.PolygonTrait, g2,
) = _polygon_polygon_process(
    g1, g2;
    WITHIN_POLYGON_ALLOWS...,
    WITHIN_REQUIRES...,
    WITHIN_EXACT...,
)

# Polygons cannot be within any curves
_within(
    ::GI.PolygonTrait, g1,
    ::GI.AbstractCurveTrait, g2,
) = false


# # Geometries within multi-geometry/geometry collections

#= Geometry is within a multi-geometry or a collection if the geometry is within
at least one of the collection elements. =#
function _within(
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g2,
)
    for sub_g2 in GI.getgeom(g2)
        within(g1, sub_g2) && return true
    end
    return false
end

# # Multi-geometry/geometry collections within geometries

#= Multi-geometry or a geometry collection is within a geometry if all
elements of the collection are within the geometry. =#
function _within(
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g1,
    ::GI.AbstractGeometryTrait, g2,
)
    for sub_g1 in GI.getgeom(g1)
        !within(sub_g1, g2) && return false
    end
    return true
end
