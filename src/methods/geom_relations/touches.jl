# # Touches

export touches

#=
## What is touches?

The touches function checks if one geometry touches another geometry. In other
words, the interiors of the two geometries don't interact, but one of the
geometries must have a boundary point that interacts with either the other
geometry's interior or boundary.


To provide an example, consider these two lines:
```@example touches
import GeometryOps as GO
import GeoInterface as GI
using Makie
using CairoMakie

l1 = GI.Line([(0.0, 0.0), (1.0, 0.0)])
l2 = GI.Line([(1.0, 0.0), (1.0, -1.0)])

f, a, p = lines(GI.getpoint(l1))
lines!(GI.getpoint(l2))
f
```
We can see that these two lines touch only at their endpoints.
```@example touches
GO.touches(l1, l2)  # true
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
const TOUCHES_REQUIRES = (in_require = false, on_require = true, out_require = false)
const TOUCHES_EXACT = (exact = False(),)

"""
    touches(geom1, geom2)::Bool

Return `true` if the first geometry touches the second geometry. In other words,
the two interiors cannot interact, but one of the geometries must have a
boundary point that interacts with either the other geometry's interior or
boundary.

## Examples
```jldoctest setup=:(using GeometryOps, GeometryBasics)
import GeometryOps as GO, GeoInterface as GI

l1 = GI.Line([(0.0, 0.0), (1.0, 0.0)])
l2 = GI.Line([(1.0, 1.0), (1.0, -1.0)])

GO.touches(l1, l2)
# output
true
```
"""
touches(g1, g2)::Bool = _touches(trait(g1), g1, trait(g2), g2)

"""
    touches(g1)

Return a function that checks if its input touches `g1`.
This is equivalent to `x -> touches(x, g1)`.
"""
touches(g1) = Base.Fix2(touches, g1)

# # Convert features to geometries
_touches(::GI.FeatureTrait, g1, ::Any, g2) = touches(GI.geometry(g1), g2)
_touches(::Any, g1, t2::GI.FeatureTrait, g2) = touches(g1, GI.geometry(g2))
_touches(::FeatureTrait, g1, ::FeatureTrait, g2) = touches(GI.geometry(g1), GI.geometry(g2))

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
    TOUCHES_EXACT...,
)

#= Geometry touches a point if the point is on the geometry boundary. =#
_touches(
    trait1::Union{GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    trait2::GI.PointTrait, g2,
) = _touches(trait2, g2, trait1, g1)


# # Lines touching geometries

#= Linestring touches another line if at least one boundary point interacts with
the boundary of interior of the other line, but the interiors don't interact. =#
_touches(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _line_curve_process(
    g1, g2;
    TOUCHES_CURVE_ALLOWED...,
    TOUCHES_REQUIRES...,
    TOUCHES_EXACT...,
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
    TOUCHES_REQUIRES...,
    TOUCHES_EXACT...,
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
    TOUCHES_EXACT...,
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
made up of interior points and no boundary points =#
_touches(
    ::GI.LinearRingTrait, g1,
    ::GI.LinearRingTrait, g2,
) = false

#= Linearring touches a polygon if at least one of the points of the ring
interact with the polygon boundary and non are in the polygon interior. =#
_touches(
    ::GI.LinearRingTrait, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    TOUCHES_POLYGON_ALLOWS...,
    TOUCHES_REQUIRES...,
    TOUCHES_EXACT...,
    closed_line = true,
)


# # Polygons touch geometries

#= Polygon touches a curve if at least one of the curve boundary points interacts
with the polygon's boundary and no curve points interact with the interior.=#
_touches(
    trait1::GI.PolygonTrait, g1,
    trait2::GI.AbstractCurveTrait, g2
) = _touches(trait2, g2, trait1, g1)


#= Polygon touches another polygon if they share at least one boundary point and
no interior points. =#
_touches(
    ::GI.PolygonTrait, g1,
    ::GI.PolygonTrait, g2,
) = _polygon_polygon_process(
    g1, g2;
    TOUCHES_POLYGON_ALLOWS...,
    TOUCHES_REQUIRES...,
    TOUCHES_EXACT...,
)

# # Geometries touch multi-geometry/geometry collections

#= 

A geometry touches a multi-geometry or a collection if the geometry touches at
least one of the elements of the collection.

This is a bit tricky to implement - we have to actually check every geometry, 
and make sure that each geom is either disjoint or touching.

Problem here is that we would end up doing double the work.

Either you check disjointness first, and then check touches - in which case
you have already done the work for the touches check, but can't take advantage of it.

Or you check touches first, and if that is false, you check disjointness.  But if touches failed,
and you don't know _why_ it was false (disjoint or contained / intersecting), you have to iterate
over every point twice -- again!


At this point we actually need a fast return function...or some more detail returned from the process functions.

That's a project for later though.  Right now we need to get this correct, so I'm going to do the dumb thing.

=#
function _touches(
    ::Union{GI.PointTrait, GI.AbstractCurveTrait, GI.PolygonTrait}, g1,
    ::Union{
        GI.MultiPointTrait, GI.AbstractMultiCurveTrait,
        GI.MultiPolygonTrait, GI.GeometryCollectionTrait,
    }, g2,
)
    has_touched = false
    for sub_g2 in GI.getgeom(g2)
        if touches(g1, sub_g2)
            has_touched = true
        else 
            # if not touching, they are either intersecting or disjoint
            # if disjoint, then we can continue
            # else, we can short circuit, since the geoms are not touching and not disjoint
            # i.e. they are intersecting
            disjoint(g1, sub_g2) || return false
        end
    end
    return has_touched
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
    has_touched = false
    for sub_g1 in GI.getgeom(g1)
        if touches(sub_g1, g2)
            has_touched = true
        else 
            # if not touching, they are either intersecting or disjoint
            # if disjoint, then we can continue
            # else, we can short circuit, since the geoms are not touching and not disjoint
            disjoint(sub_g1, g2) || return false
        end
    end
    return has_touched
end

# Extent forwarding


function _touches(t1::GI.AbstractGeometryTrait, g1, t2, e::Extents.Extent)
    return _touches(t1, g1, GI.PolygonTrait(), extent_to_polygon(e))
end
function _touches(t1, e1::Extents.Extent, t2::GI.AbstractGeometryTrait, g2)
    return _touches(GI.PolygonTrait(), extent_to_polygon(e1), t2, g2)
end
function _touches(t1, e1::Extents.Extent, t2, e2::Extents.Extent)
    return Extents.touches(e1, e2)
end


