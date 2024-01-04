# # Crosses

export crosses

#=
## What is crosses?

The crosses function checks if one geometry is crosses another geometry.
A geometry can only cross another geometry if they are either two lines, or if
the two geometries have different dimensionalities. If checking two lines, they
must meet in one point. If checking two geometries of different dimensions, the
interiors must meet in at least one point and at least one of the geometries
must have a point outside of the other geometry.

Note that points can't cross any geometries, despite different dimension, due to
their inability to be both interior and exterior to any other shape.

To provide an example, consider these two lines:
```@example crosses
using GeometryOps
using GeometryOps.GeometryBasics
using Makie
using CairoMakie

l1 = Line([Point(0.0, 0.0), Point(1.0, 0.0)])
l2 = Line([Point(0.5, 1.0), Point(0.5, -1.0)])

f, a, p = lines(l1)
lines!(l2)
```
We can see that these two lines cross at their midpoints.
```@example crosses
crosses(l1, l2)  # true
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait.

Each of these calls a method in the geom_geom_processors file. The methods in
this file determine if the given geometries meet a set of criteria. For the
`crosses` function and arguments g1 and g2, this criteria is as follows:
    - points of g1 are allowed to be in the interior of g2 (only through
    crossing and NOT overlap for lines)
    - points of g1 are allowed to be on the boundary of g2
    - points of g1 are allowed to be in the exterior of g2
    - at least one point of g1 are required to be in the interior of g2
    - no points of g1 are required to be on the boundary of g2
    - at least one point of g1 are required to be in the exterior of g2

The code for the specific implementations is in the geom_geom_processors file.
=#

const CROSSES_CURVE_ALLOWS = (over_allow = false, cross_allow = true, on_allow = true, out_allow = true)
const CROSSES_POLYGON_ALLOWS = (in_allow = true, on_allow = true, out_allow = true)
const CROSSES_REQUIRES = (in_require = true, on_require = false, out_require = true)

"""
    crosses(geom1, geom2)::Bool

Return `true` if the first geometry crosses the second geometry. If they are
both lines, they must meet in one point. Otherwise, they must be of different
dimensions, the interiors must intersect, and the interior of the first geometry
must intersect the exterior of the secondary geometry.

## Examples
```jldoctest setup=:(using GeometryOps, GeometryBasics)
import GeometryOps as GO, GeoInterface as GI

l1 = GI.Line([(0.0, 0.0), (1.0, 0.0)])
l2 = GI.Line([(0.5, 1.0), (0.5, -1.0)])

GO.crosses(l1, l2)
# output
true
```
"""
crosses(g1, g2) = _crosses(trait(g1), g1, trait(g2), g2)

# # Convert features to geometries
_crosses(::GI.FeatureTrait, g1, ::Any, g2) = crosses(GI.geometry(g1), g2)
_crosses(::Any, g1, t2::GI.FeatureTrait, g2) = crosses(g1, GI.geometry(g2))


# # Non-specified geometries 

# Points and geometries with the same dimensions D where D ≂̸ 1 default to false
_crosses(
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g2,
) = false


# # Lines cross geometries

#= Linestring crosses another linestring if the intersection of the two lines
is exlusivly points (only cross intersections) =#
_crosses(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _line_curve_process(
    g1, g2;
    CROSSES_CURVE_ALLOWS...,
    CROSSES_REQUIRES...,
    closed_line = false,
    closed_curve = false,
)

#= Linestring crosses a linearring if the intersection of the line and ring is
exlusivly points (only cross intersections) =#
_crosses(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    CROSSES_CURVE_ALLOWS...,
    CROSSES_REQUIRES...,
    closed_line = false,
    closed_curve = true,
)

#= Linestring crosses a polygon if at least some of the line interior is in the
polygon interior and some of the line interior is exterior to the polygon. =#
_crosses(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    CROSSES_POLYGON_ALLOWS...,
    CROSSES_REQUIRES...,
    closed_line = false,
)


# # Rings cross geometries

#= Linearring crosses a linestring if the intersection of the line and ring is
exlusivly points (only cross intersections) =#
_crosses(
    trait1::GI.LinearRingTrait, g1,
    trait2::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _crosses(trait2, g2, trait1, g1)

#= Linearring crosses another ring if the intersection of the two rings is
exlusivly points (only cross intersections) =#
_crosses(
    ::GI.LinearRingTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    CROSSES_CURVE_ALLOWS...,
    CROSSES_REQUIRES...,
    closed_line = true,
    closed_curve = true,
)

#= Linearring crosses a polygon if at least some of the ring interior is in the
polygon interior and some of the ring interior is exterior to the polygon. =#
_crosses(
    ::GI.LinearRingTrait, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    CROSSES_POLYGON_ALLOWS...,
    CROSSES_REQUIRES...,
    closed_line = true,
)


# # Polygons cross geometries

#= Polygon crosses a curve if at least some of the curve interior is in the
polygon interior and some of the curve interior is exterior to the polygon.=#
_crosses(
    trait1::GI.PolygonTrait, g1,
    trait2::GI.AbstractCurveTrait, g2
) = _crosses(trait2, g2, trait1, g1)


# # Geometries cross multi-geometry/geometry collections

#= Geometry crosses a multi-geometry or a collection if the geometry crosses
one of the elements of the collection. =#
function _crosses(
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g2,
)
    for sub_g2 in GI.getgeom(g2)
        crosses(g1, sub_g2) && return true
    end
    return false
end

# # Multi-geometry/geometry collections cross geometries

#= Multi-geometry or a geometry collection crosses a geometry one elements of
the collection crosses the geometry. =#
function _crosses(
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g1,
    ::GI.AbstractGeometryTrait, g2,
)
    for sub_g1 in GI.getgeom(g1)
        crosses(sub_g1, g2) && return true
    end
    return false
end