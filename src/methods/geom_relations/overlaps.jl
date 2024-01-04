# # Overlaps

export overlaps

#=
## What is overlaps?

The overlaps function checks if two geometries overlap. Two geometries overlap
if they have the same dimension, and if they overlap then their interiors
interact, but they both also need interior points exterior to the other
geometry. 

Note that this means it is impossible for a single point to overlap with a
single point and a line only overlaps with another line if only a section of
each line is colinear (crosses don't count for interior points interacting). 

To provide an example, consider these two lines:
```@example overlaps
using GeometryOps
using GeometryOps.GeometryBasics
using Makie
using CairoMakie

l1 = GI.LineString([(0.0, 0.0), (0.0, 10.0)])
l2 = GI.LineString([(0.0, -10.0), (0.0, 3.0)])
f, a, p = lines(GI.getpoint(l1), color = :blue)
scatter!(GI.getpoint(l1), color = :blue)
lines!(GI.getpoint(l2), color = :orange)
scatter!(GI.getpoint(l2), color = :orange)
```
We can see that the two lines overlap in the plot:
```@example overlaps
overlap(l1, l2)
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait.

Each of these calls a method in the geom_geom_processors file. The methods in
this file determine if the given geometries meet a set of criteria. For the
`overlaps` function and arguments g1 and g2, this criteria is as follows:
    - points of g1 are allowed to be in the interior of g2
    - points of g1 are allowed to be on the boundary of g2
    - points of g1 are allowed to be in the exterior of g2
    - at least one point of g1 is required to be in the interior of g2
    - at least one point of g2 is required to be in the interior of g1
    - no points of g1 is required to be on the boundary of g2
    - at least one point of g1 is required to be in the exterior of g2
    - at least one point of g2 is required to be in the exterior of g1

The code for the specific implementations is in the geom_geom_processors file.
=#

const OVERLAPS_CURVE_ALLOWS = (over_allow = true, cross_allow = false, on_allow = true, out_allow = true)
const OVERLAPS_POLYGON_ALLOWS = (in_allow = true, on_allow = true, out_allow = true)
const OVERLAPS_REQUIRES = (in_require = true, on_require = false, out_require = true)

"""
    overlaps(geom1, geom2)::Bool

Compare two Geometries of the same dimension and return true if their interiors
interact, but they both also have interior points exterior to the other
geometry. Lines crossing doesn't count for interiors interacting as overlaps
of curves must be of dimension one. 

## Examples
```jldoctest
import GeometryOps as GO, GeoInterface as GI
poly1 = GI.Polygon([[(0,0), (0,5), (5,5), (5,0), (0,0)]])
poly2 = GI.Polygon([[(1,1), (1,6), (6,6), (6,1), (1,1)]])

GO.overlaps(poly1, poly2)
# output
true
```
"""
overlaps(g1, g2)::Bool = _overlaps(GI.trait(g1), g1, GI.trait(g2), g2)


# # Convert features to geometries
_overlaps(::GI.FeatureTrait, g1, ::Any, g2) = overlaps(GI.geometry(g1), g2)
_overlaps(::Any, g1, t2::GI.FeatureTrait, g2) = overlaps(g1, GI.geometry(g2))


# # Non-specified geometries 

# Geometries of different dimensions and points cannot overlap and return false
_overlaps(
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g2,
) = false


# # Lines cross curves

#= Linestring overlaps with another linestring when they share co-linear
segments (interiors interacting), but both have interior points exterior to the
other line. =#
_overlaps(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _line_curve_process(
    g1, g2;
    OVERLAPS_CURVE_ALLOWS...,
    OVERLAPS_REQUIRES...,
    closed_line = false,
    closed_curve = false,
) && _line_curve_process(
        g2, g1;
        OVERLAPS_CURVE_ALLOWS...,
        OVERLAPS_REQUIRES...,
        closed_line = false,
        closed_curve = false,
    )

#= Linestring overlaps with a linearring when they share co-linear segments
(interiors interacting), but both have interior points exterior to the other. =#
_overlaps(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    OVERLAPS_CURVE_ALLOWS...,
    OVERLAPS_REQUIRES...,
    closed_line = false,
    closed_curve = true,
) && _line_curve_process(
        g2, g1;
        OVERLAPS_CURVE_ALLOWS...,
        OVERLAPS_REQUIRES...,
        closed_line = true,
        closed_curve = false,
    )


# # Rings cross curves

#= Linearring overlaps with a linestring when they share co-linear segments
(interiors interacting), but both have interior points exterior to the other. =#
_overlaps(
    trait1::GI.LinearRingTrait, g1,
    trait2::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _overlaps(trait2, g2, trait1, g1)

#= Linearring overlaps with another linearring when they share co-linear
segments (interiors interacting), but both have interior points exterior to the
other line. =#
_overlaps(
    ::GI.LinearRingTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    OVERLAPS_CURVE_ALLOWS...,
    OVERLAPS_REQUIRES...,
    closed_line = true,
    closed_curve = true,
) && _line_curve_process(
        g2, g1;
        OVERLAPS_CURVE_ALLOWS...,
        OVERLAPS_REQUIRES...,
        closed_line = true,
        closed_curve = true,
    )


# # Polygons cross polygons

#= Polygon overlaps with another polygon when their interiors intersect, but
both have interior points exterior to the other polygon. =#
_overlaps(
    ::GI.PolygonTrait, g1,
    ::GI.PolygonTrait, g2,
) = _polygon_polygon_process(
    g1, g2;
    OVERLAPS_POLYGON_ALLOWS...,
    OVERLAPS_REQUIRES...,
) && _polygon_polygon_process(
    g2, g1;
    OVERLAPS_POLYGON_ALLOWS...,
    OVERLAPS_REQUIRES...,
)

# # Geometries disjoint multi-geometry/geometry collections

# Multipoints overlap with other multipoints if only some sub-points are shared
function _overlaps(
    ::GI.MultiPointTrait, g1,
    ::GI.MultiPointTrait, g2,
)
    one_diff = false  # assume that all the points are the same
    one_same = false  # assume that all points are different
    for p1 in GI.getpoint(g1)
        match_point = false
        for p2 in GI.getpoint(g2)
            if equals(p1, p2)  # Point is shared
                one_same = true
                match_point = true
                break
            end
        end
        one_diff |= !match_point  # Point isn't shared
        one_same && one_diff && return true
    end
    return false
end

#= Geometry overlaps a multi-geometry or a collection if the geometry overlaps
at least one of the elements of the collection. =#
function _overlaps(
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g2,
)
    for sub_g2 in GI.getgeom(g2)
        overlaps(g1, sub_g2) && return true
    end
    return false
end

# # Multi-geometry/geometry collections cross geometries

#= Multi-geometry or a geometry collection overlaps a geometry if at least one
elements of the collection overlaps the geometry. =#
function _overlaps(
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g1,
    ::GI.AbstractGeometryTrait, g2,
)
    for sub_g1 in GI.getgeom(g1)
        overlaps(sub_g1, g2) && return true
    end
    return false
end