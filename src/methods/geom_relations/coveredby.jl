# # CoveredBy

export coveredby

#=
## What is coveredby?

The coveredby function checks if one geometry is covered by another geometry.
This is an extension of within that does not require the interiors of the two
geometries to intersect, but still does require that the interior and boundary
of the first geometry isn't outside of the second geometry. 

To provide an example, consider this point and line:
```@example coveredby
import GeometryOps as GO
import GeoInterface as GI
using Makie
using CairoMakie

p1 = (0.0, 0.0)
l1 = GI.Line([p1, (1.0, 1.0)])
f, a, p = lines(GI.getpoint(l1))
scatter!(p1, color = :red)
f
```
As we can see, `p1` is on the endpoint of l1. This means it is not `within`, but
it does meet the definition of `coveredby`.
```@example coveredby
GO.coveredby(p1, l1)  # true
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait.

Each of these calls a method in the geom_geom_processors file. The methods in
this file determine if the given geometries meet a set of criteria. For the
`coveredby` function and arguments g1 and g2, this criteria is as follows:
    - points of g1 are allowed to be in the interior of g2 (either through
    overlap or crossing for lines)
    - points of g1 are allowed to be on the boundary of g2
    - points of g1 are not allowed to be in the exterior of g2
    - no points of g1 are required to be in the interior of g2
    - no points of g1 are required to be on the boundary of g2
    - no points of g1 are required to be in the exterior of g2

The code for the specific implementations is in the geom_geom_processors file.
=#
const COVEREDBY_ALLOWS = (in_allow = true, on_allow = true, out_allow = false)
const COVEREDBY_CURVE_ALLOWS = (over_allow = true, cross_allow = true, on_allow = true, out_allow = false)
const COVEREDBY_CURVE_REQUIRES = (in_require = false, on_require = false, out_require = false)
const COVEREDBY_POLYGON_REQUIRES = (in_require = true, on_require = false, out_require = false,)
const COVEREDBY_EXACT = (exact = _False(),)

"""
    coveredby(g1, g2)::Bool

Return `true` if the first geometry is completely covered by the second
geometry. The interior and boundary of the primary geometry (g1) must not
intersect the exterior of the secondary geometry (g2).

Furthermore, `coveredby` returns the exact opposite result of `covers`. They are
equivalent with the order of the arguments swapped.

## Examples
```jldoctest setup=:(using GeometryOps, GeometryBasics)
import GeometryOps as GO, GeoInterface as GI
p1 = GI.Point(0.0, 0.0)
p2 = GI.Point(1.0, 1.0)
l1 = GI.Line([p1, p2])

GO.coveredby(p1, l1)
# output
true
```
"""
coveredby(g1, g2) = _coveredby(trait(g1), g1, trait(g2), g2)

# # Convert features to geometries
_coveredby(::GI.FeatureTrait, g1, ::Any, g2) = coveredby(GI.geometry(g1), g2)
_coveredby(::Any, g1, t2::GI.FeatureTrait, g2) = coveredby(g1, GI.geometry(g2))
_coveredby(::FeatureTrait, g1, ::FeatureTrait, g2) = coveredby(GI.geometry(g1), GI.geometry(g2))

# # Points coveredby geometries

# Point is coveredby another point if those points are equal
_coveredby(
    ::GI.PointTrait, g1,
    ::GI.PointTrait, g2,
) = equals(g1, g2)

# Point is coveredby a line/linestring if it is on a line vertex or an edge
_coveredby(
    ::GI.PointTrait, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _point_curve_process(
    g1, g2;
    COVEREDBY_ALLOWS...,
    closed_curve = false,
)

# Point is coveredby a linearring if it is on a vertex or an edge of ring
_coveredby(
    ::GI.PointTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _point_curve_process(
    g1, g2;
    COVEREDBY_ALLOWS...,
    closed_curve = true,
)

# Point is coveredby a polygon if it is inside polygon, including edges/vertices
_coveredby(
    ::GI.PointTrait, g1,
    ::GI.PolygonTrait, g2,
) = _point_polygon_process(
    g1, g2;
    COVEREDBY_ALLOWS...,
    COVEREDBY_EXACT...,
)

# Points cannot cover any geometry other than points
_coveredby(
    ::Union{GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    ::GI.PointTrait, g2,
) = false


# # Lines coveredby geometries

#= Linestring is coveredby a line if all interior and boundary points of the
first line are on the interior/boundary points of the second line. =#
_coveredby(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _line_curve_process(
    g1, g2;
    COVEREDBY_CURVE_ALLOWS...,
    COVEREDBY_CURVE_REQUIRES...,
    COVEREDBY_EXACT...,
    closed_line = false,
    closed_curve = false,
)

#= Linestring is coveredby a ring if all interior and boundary points of the
line are on the edges of the ring. =#
_coveredby(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    COVEREDBY_CURVE_ALLOWS...,
    COVEREDBY_CURVE_REQUIRES...,
    COVEREDBY_EXACT...,
    closed_line = false,
    closed_curve = true,
)

#= Linestring is coveredby a polygon if all interior and boundary points of the
line are in the polygon interior or on its edges, including hole edges. =#
_coveredby(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    COVEREDBY_ALLOWS...,
    COVEREDBY_CURVE_REQUIRES...,
    COVEREDBY_EXACT...,
    closed_line = false,
)

# # Rings covered by geometries

#= Linearring is covered by a line if all vertices and edges of the ring are on
the edges and vertices of the line. =#
_coveredby(
    ::GI.LinearRingTrait, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _line_curve_process(
    g1, g2;
    COVEREDBY_CURVE_ALLOWS...,
    COVEREDBY_CURVE_REQUIRES...,
    COVEREDBY_EXACT...,
    closed_line = true,
    closed_curve = false,
)

#= Linearring is covered by another linear ring if all vertices and edges of the
first ring are on the edges/vertices of the second ring. =#
_coveredby(
    ::GI.LinearRingTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    COVEREDBY_CURVE_ALLOWS...,
    COVEREDBY_CURVE_REQUIRES...,
    COVEREDBY_EXACT...,
    closed_line = true,
    closed_curve = true,
)

#= Linearring is coveredby a polygon if all vertices and edges of the ring are
in the polygon interior or on the polygon edges, including hole edges. =# 
_coveredby(
    ::GI.LinearRingTrait, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    COVEREDBY_ALLOWS...,
    COVEREDBY_CURVE_REQUIRES...,
    COVEREDBY_EXACT...,
    closed_line = true,
)


# # Polygons covered by geometries

#= Polygon is covered by another polygon if if the interior and edges of the
first polygon are in the second polygon interior or on polygon edges, including
hole edges.=#
_coveredby(
    ::GI.PolygonTrait, g1,
    ::GI.PolygonTrait, g2,
) = _polygon_polygon_process(
    g1, g2;
    COVEREDBY_ALLOWS...,
    COVEREDBY_POLYGON_REQUIRES...,
    COVEREDBY_EXACT...,
)

# Polygons cannot covered by any curves
_coveredby(
    ::GI.PolygonTrait, g1,
    ::GI.AbstractCurveTrait, g2,
) = false


# # Geometries coveredby multi-geometry/geometry collections

#= Geometry is covered by a multi-geometry or a collection if one of the elements
of the collection cover the geometry. =#
function _coveredby(
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g2,
)
    for sub_g2 in GI.getgeom(g2)
        coveredby(g1, sub_g2) && return true
    end
    return false
end

# # Multi-geometry/geometry collections coveredby geometries

#= Multi-geometry or a geometry collection is covered by a geometry if all
elements of the collection are covered by the geometry. =#
function _coveredby(
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g1,
    ::GI.AbstractGeometryTrait, g2,
)
    for sub_g1 in GI.getgeom(g1)
        !coveredby(sub_g1, g2) && return false
    end
    return true
end
