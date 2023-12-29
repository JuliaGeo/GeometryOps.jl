# # Touches

export touches

#=
## What is touches?

The touches function checks if one geometry touches another geometry. In other
words, the interiors of the two geometries don't interact, but one of the
geometries must have a boundary point that interacts with either the other
geometies interior or boundary.


To provide an example, consider these two lines:
```@example touches
using GeometryOps
using GeometryOps.GeometryBasics
using Makie
using CairoMakie

l1 = Line([Point(0.0, 0.0), Point(1.0, 0.0)])
l2 = Line([Point(1.0, 0.0), Point(1.0, -1.0)])

f, a, p = lines(l1)
lines!(l2)
```
We can see that these two lines touch only at their endpoints.
```@example touches
touches(l1, l2)  # true
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait.

Each of these calls a method in the geom_geom_processors file. The methods in
this file determine if the given geometries meet a set of criteria. For the
`touches` function and arguments g1 and g2, this criteria is as follows:
    - points of g1 are not allowed to be in the interior of g2
    - points of g1 are allowed to be on the boundary of g2
    - points of g1 are allowed to be in the exterior of g2
    - no points of g1 are required to be in the interior of g2
    - at least one point of g1 is required to be on the boundary of g2
    - no points of g1 are required to be in the exterior of g2

The code for the specific implementations is in the geom_geom_processors file.
=#

const TOUCHES_POINT_ALLOWED = (in_allow = false, on_allow = true, out_allow = false)
const TOUCHES_CURVE_ALLOWED = (over_allow = false, cross_allow = false, on_allow = true, out_allow = true)
const TOUCHES_POLYGON_ALLOWS = (in_allow = false, on_allow = true, out_allow = true)
const TOUCHES_REQUIRED = (in_require = false, on_require = true, out_require = false)

"""
    touches(geom1, geom2)::Bool

Return `true` if the first geometry touches the second geometry. In other words,
the two interiors cannot interact, but one of the geometries must have a
boundary point that interacts with either the other geometies interior or
boundary.

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
touches(g1, g2)::Bool = _touches(trait(g1), g1, trait(g2), g2)

# # Convert features to geometries
_touches(::GI.FeatureTrait, g1, ::Any, g2) = touches(GI.geometry(g1), g2)
_touches(::Any, g1, t2::GI.FeatureTrait, g2) = touches(g1, GI.geometry(g2))


# # Point touches geometries

# Point cannot touch another point as if they are equal, interiors interact
_touches(
    ::GI.PointTrait, g1,
    ::GI.PointTrait, g2,
) = false

# Point touches a linestring if it equal to the first of last point of the line
function _touches(
    ::GI.PointTrait, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
)
    n = GI.npoint(g2)
    p1 = GI.getpoint(g2, 1)
    pn = GI.getpoint(g2, n)
    equals(p1, pn) && return false
    return equals(g1, p1) || equals(g1, pn)
end

# Point cannot 'touch' a linearring given that the ring has no boundary points
_touches(
    ::GI.PointTrait, g1,
    ::GI.LinearRingTrait, g2,
) = false

# Point touches a polygon if it is on the boundary of that polygon
_touches(
    ::GI.PointTrait, g1,
    ::GI.PolygonTrait, g2,
) = _point_polygon_process(
    g1, g2;
    TOUCHES_POINT_ALLOWED...,
)

#= Geometry touches a point if the point is on the geometry boundary. =#
_touches(
    trait1::Union{GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    trait2::GI.PointTrait, g2,
) = _touches(trait2, g2, trait1, g1)


# # Lines touching geometries

#= Linestring touches another line if at least one bounday point interacts with
the bounday of interior of the other line, but the interiors don't interact. =#
_touches(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _line_curve_process(
    g1, g2;
    TOUCHES_CURVE_ALLOWED...,
    TOUCHES_REQUIRED...,
    closed_line = false,
    closed_curve = false,
)


#= Linestring touches a linearring if at least one of the boundary points of the
line interacts with the linear ring, but their interiors can't interact. =#
_touches(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    TOUCHES_CURVE_ALLOWED...,
    TOUCHES_REQUIRED...,
    closed_line = false,
    closed_curve = true,
)

#= Linestring touches a polygon if at least one of the boundary points of the
line interacts with the boundary of the polygon. =#
_touches(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    TOUCHES_POLYGON_ALLOWS...,
    TOUCHES_REQUIRES...,
    closed_line = false,
)


# # Rings touch geometries

#= Linearring touches a linestring if at least one of the boundary points of the
line interacts with the linear ring, but their interiors can't interact. =#
_touches(
    trait1::GI.LinearRingTrait, g1,
    trait2::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _touches(trait2, g2, trait1, g1)

#= Linearring cannot touch another linear ring since they are both exclusively
made up of interior points and no bounday points =#
_touches(
    ::GI.LinearRingTrait, g1,
    ::GI.LinearRingTrait, g2,
) = false

#= Linearring touches a polygon if at least one of the points of the ring
interact with the polygon bounday and non are in the polygon interior. =#
_touches(
    ::GI.LinearRingTrait, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    TOUCHES_POLYGON_ALLOWS...,
    TOUCHES_REQUIRES...,
    closed_line = true,
)


# # Polygons touch geometries

#= Polygon touches a curve if at least one of the curve bounday points interacts
with the polygon's bounday and no curve points interact with the interior.=#
_touches(
    trait1::GI.PolygonTrait, g1,
    trait2::GI.AbstractCurveTrait, g2
) = _touches(trait2, g2, trait1, g1)


# # Geometries touch multi-geometry/geometry collections

#= Geometry touch a multi-geometry or a collection if the geometry touches at
least one of the elements of the collection. =#
function _touches(
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g2,
)
    for sub_g2 in GI.getgeom(g2)
        touches(g1, sub_g2) && return true
    end
    return false
end

# # Multi-geometry/geometry collections cross geometries

#= Multi-geometry or a geometry collection touches a geometry if at least one
elements of the collection touches the geometry. =#
function _touches(
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g1,
    ::GI.AbstractGeometryTrait, g2,
)
    for sub_g1 in GI.getgeom(g1)
        touches(sub_g1, g2) && return true
    end
    return false
end