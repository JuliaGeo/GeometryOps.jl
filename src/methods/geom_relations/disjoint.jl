# # Disjoint

export disjoint
#=
## What is disjoint?

The disjoint function checks if one geometry is outside of another geometry,
without sharing any boundaries or interiors.

To provide an example, consider these two lines:
```@example disjoint
import GeometryOps as GO
import GeoInterface as GI
using Makie
using CairoMakie

l1 = GI.LineString([(0.0, 0.0), (1.0, 0.0), (0.0, 0.1)])
l2 = GI.LineString([(2.0, 0.0), (2.75, 0.0)])
f, a, p = lines(GI.getpoint(l1), color = :blue)
scatter!(GI.getpoint(l1), color = :blue)
lines!(GI.getpoint(l2), color = :orange)
scatter!(GI.getpoint(l2), color = :orange)
f
```
We can see that none of the edges or vertices of l1 interact with l2 so they are
disjoint.
```@example disjoint
GO.disjoint(l1, l2)  # returns true
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait.

Each of these calls a method in the geom_geom_processors file. The methods in
this file determine if the given geometries meet a set of criteria. For the
`disjoint` function and arguments g1 and g2, this criteria is as follows:
    - points of g1 are not allowed to be in the interior of g2
    - points of g1 are not allowed to be on the boundary of g2
    - points of g1 are allowed to be in the exterior of g2
    - no points required to be in the interior of g2
    - no points of g1 are required to be on the boundary of g2
    - no points of g1 are required to be in the exterior of g2

The code for the specific implementations is in the geom_geom_processors file.
=#

const DISJOINT_ALLOWS = (in_allow = false, on_allow = false, out_allow = true)
const DISJOINT_CURVE_ALLOWS = (over_allow = false, cross_allow = false, on_allow = false, out_allow = true)
const DISJOINT_REQUIRES = (in_require = false, on_require = false, out_require = false)
const DISJOINT_EXACT = (exact = _False(),)

"""
    disjoint(geom1, geom2)::Bool

Return `true` if the first geometry is disjoint from the second geometry.

Return `true` if the first geometry is disjoint from the second geometry. The
interiors and boundaries of both geometries must not intersect.

## Examples
```jldoctest setup=:(using GeometryOps, GeoInterface)
import GeometryOps as GO, GeoInterface as GI

line = GI.LineString([(1, 1), (1, 2), (1, 3), (1, 4)])
point = (2, 2)
GO.disjoint(point, line)

# output
true
```
"""
disjoint(g1, g2) = _disjoint(trait(g1), g1, trait(g2), g2)

# # Convert features to geometries
_disjoint(::FeatureTrait, g1, ::Any, g2) = disjoint(GI.geometry(g1), g2)
_disjoint(::Any, g1, ::FeatureTrait, g2) = disjoint(g1, geometry(g2))
_disjoint(::FeatureTrait, g1, ::FeatureTrait, g2) = disjoint(GI.geometry(g1), GI.geometry(g2))

# # Point disjoint geometries

# Point is disjoint from another point if the points are not equal.
_disjoint(
    ::GI.PointTrait, g1,
    ::GI.PointTrait, g2,
) = !equals(g1, g2)

# Point is disjoint from a linestring if it is not on the line's edges/vertices.
_disjoint(
    ::GI.PointTrait, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _point_curve_process(
    g1, g2;
    DISJOINT_ALLOWS...,
    closed_curve = false,
)

# Point is disjoint from a linearring if it is not on the ring's edges/vertices. 
_disjoint(
    ::GI.PointTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _point_curve_process(
    g1, g2;
    DISJOINT_ALLOWS...,
    closed_curve = true,
)

#= Point is disjoint from a polygon if it is not on any edges, vertices, or
within the polygon's interior. =#
_disjoint(
    ::GI.PointTrait, g1,
    ::GI.PolygonTrait, g2,
) = _point_polygon_process(
    g1, g2;
    DISJOINT_ALLOWS...,
    DISJOINT_EXACT...,
)

#= Geometry is disjoint from a point if the point is not in the interior or on
the boundary of the geometry. =#
_disjoint(
    trait1::Union{GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    trait2::GI.PointTrait, g2,
) = _disjoint(trait2, g2, trait1, g1)


# # Lines disjoint geometries

#= Linestring is disjoint from another line if they do not share any interior
edge/vertex points or boundary points. =#
_disjoint(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _line_curve_process(
    g1, g2;
    DISJOINT_CURVE_ALLOWS...,
    DISJOINT_REQUIRES...,
    DISJOINT_EXACT...,
    closed_line = false,
    closed_curve = false,
)

#= Linestring is disjoint from a linearring if they do not share any interior
edge/vertex points or boundary points. =#
_disjoint(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    DISJOINT_CURVE_ALLOWS...,
    DISJOINT_REQUIRES...,
    DISJOINT_EXACT...,
    closed_line = false,
    closed_curve = true,
)

#= Linestring is disjoint from a polygon if the interior and boundary points of
the line are not in the polygon's interior or on the polygon's boundary. =# 
_disjoint(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    DISJOINT_ALLOWS...,
    DISJOINT_REQUIRES...,
    DISJOINT_EXACT...,
    closed_line = false,
)

#= Geometry is disjoint from a linestring if the line's interior and boundary
points don't intersect with the geometrie's interior and boundary points. =#
_disjoint(
    trait1::Union{GI.LinearRingTrait, GI.PolygonTrait}, g1,
    trait2::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _disjoint(trait2, g2, trait1, g1)


# # Rings disjoint geometries

#= Linearrings is disjoint from another linearring if they do not share any
interior edge/vertex points or boundary points.=#
_disjoint(
    ::GI.LinearRingTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    DISJOINT_CURVE_ALLOWS...,
    DISJOINT_REQUIRES...,
    DISJOINT_EXACT...,
    closed_line = true,
    closed_curve = true,
)

#= Linearring is disjoint from a polygon if the interior and boundary points of
the ring are not in the polygon's interior or on the polygon's boundary. =# 
_disjoint(
    ::GI.LinearRingTrait, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    DISJOINT_ALLOWS...,
    DISJOINT_REQUIRES...,
    DISJOINT_EXACT...,
    closed_line = true,
)

# # Polygon disjoint geometries

#= Polygon is disjoint from another polygon if they do not share any edges or
vertices and if their interiors do not intersect, excluding any holes. =#
_disjoint(
    ::GI.PolygonTrait, g1,
    ::GI.PolygonTrait, g2,
) = _polygon_polygon_process(
    g1, g2;
    DISJOINT_ALLOWS...,
    DISJOINT_REQUIRES...,
    DISJOINT_EXACT...,
)


# # Geometries disjoint multi-geometry/geometry collections

#= Geometry is disjoint from a multi-geometry or a collection if all of the
elements of the collection are disjoint from the geometry. =#
function _disjoint(
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g2,
)
    disjoint_partial(x) = disjoint(g1, x)
    for sub_g2 in GI.getgeom(g2)
        !disjoint_partial(sub_g2) && return false
    end
    return true
end

# # Multi-geometry/geometry collections coveredby geometries

#= Multi-geometry or a geometry collection is covered by a geometry if all
elements of the collection are covered by the geometry. =#
function _disjoint(
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g1,
    ::GI.AbstractGeometryTrait, g2,
)
    disjoint_partial(x) = disjoint(x, g2)
    for sub_g1 in GI.getgeom(g1)
        !disjoint_partial(sub_g1) && return false
    end
    return true
end
