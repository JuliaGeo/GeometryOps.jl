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
using GeometryOps
using GeometryOps.GeometryBasics
using Makie
using CairoMakie

p1 = Point(0.0, 0.0)
p2 = Point(1.0, 1.0)
l1 = Line(p1, p2)

f, a, p = lines([p1, p2])
scatter!(p1, color = :red)
```
As we can see, `p1` is on the endpoint of l1. This means it is not `within`, but
it does meet the definition of `coveredby`.
```@example coveredby
coveredby(p1, l1)  # true
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait.

Each of these calls a method in the geom_geom_processors file. The methods in
this file determine if the given geometries meet a set of criteria. For the
`coveredby` function and arguments g1 and g2, this criteria is as follows:
    - points of g1 are allowed to be in the interior of g2
    - points of g1 are allowed to be on the boundary of g2
    - points of g1 are not allowed to be in the exterior of g2
    - no points of g1 are required to be in the interior of g2
    - no points of g1 are required to be on the boundary of g2
    - no points of g1 are required to be in the exterior of g2

The code for the specific implementations is in the geom_geom_processors file.
=#

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
coveredby(g1, g2) = coveredby(trait(g1), g1, trait(g2), g2)
coveredby(::GI.FeatureTrait, g1, ::Any, g2) = coveredby(GI.geometry(g1), g2)
coveredby(::Any, g1, t2::GI.FeatureTrait, g2) = coveredby(g1, GI.geometry(g2))


# Points coveredby geometries
"""
    coveredby(::GI.PointTrait, g1, ::GI.PointTrait, g2)::Bool

If a point is coveredby another point, then those points must be equal. If they
are not equal, then they are not coveredby and return false.
"""
coveredby(
    ::GI.PointTrait, g1,
    ::GI.PointTrait, g2,
) = equals(g1, g2)


"""
    coveredby(
        ::GI.PointTrait, g1,
        ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
    )::Bool

A point is coveredby a line or linestring if it is on a vertex or an edge of
that linestring. Return true if those conditions are met, else false.
"""
coveredby(
    ::GI.PointTrait, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _point_curve_process(
    g1, g2;
    in_allow = true, on_allow = true, out_allow = false,
    repeated_last_coord = false,
)

"""
    coveredby(::GI.PointTrait, g1, ::GI.LinearRingTrait, g2)::Bool

A point is coveredby a linear ring if it is on a vertex or an edge of that
linear ring. Return true if those conditions are met, else false.
"""
coveredby(
    ::GI.PointTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _point_curve_process(
    g1, g2;
    in_allow = true, on_allow = true, out_allow = false,
    repeated_last_coord = true,
)

"""
    coveredby(::GI.PointTrait, g1, ::GI.PolygonTrait, g2)::Bool

A point is coveredby a polygon if it is inside of that polygon, including edges 
and vertices. Return true if those conditions are met, else false.
"""
coveredby(
    ::GI.PointTrait, g1,
    ::GI.PolygonTrait, g2,
) = _point_polygon_process(
    g1, g2;
    in_allow = true, on_allow = true, out_allow = false,
)

"""
coveredby(::GI.AbstractTrait, g1, ::GI.PointTrait, g2)::Bool

Points cannot cover any geometry other than points. Return false if not
dispatched to more specific function.
"""
coveredby(
    ::GI.AbstractTrait, g1,
    ::GI.PointTrait, g2,
) = false

# Lines coveredby geometries
"""
    coveredby(
        ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
        ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
    )::Bool

A line or linestring is coveredby another line or linestring if all of the
interior and boundary points of the first line are on the interior and
boundary points of the second line.
"""
coveredby(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _line_curve_process(
    g1, g2;
    in_allow = true, on_allow = true, out_allow = false,
    in_require = false, on_require = false, out_require = false,
    closed_line = false,
    closed_curve = false,
)

"""
    coveredby(
        ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
        ::GI.LinearRingTrait, g2,
    )::Bool

A line or linestring is coveredby a linear ring if all of the interior and
boundary points of the line are on the edges of the ring.
"""
coveredby(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    in_allow = true, on_allow = true, out_allow = false,
    in_require = false, on_require = false, out_require = false,
    closed_line = false,
    closed_curve = true,
)

"""
    coveredby(::GI.LineStringTrait, g1, ::GI.PolygonTrait, g2)::Bool

A line or linestring is coveredby a polygon if all of the interior and boundary
points of the line are in the polygon interior or on its edges. This includes
edges of holes. Return true if those conditions are met, else false.
"""
coveredby(
    ::Union{GI.LineTrait, GI.LineStringTrait}, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    in_allow =  true, on_allow = true, out_allow = false,
    in_require = false, on_require = false, out_require = false,
    closed_line = false,
)

# Rings covered by geometries
"""
    coveredby(
        ::GI.LinearRingTrait, g1,
        ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
    )::Bool

A linear ring is covered by a linestring if all the vertices and edges of the
linear ring are on the edges/vertices of the linear ring. Return true if
those conditions are met, else false.
"""
coveredby(
    ::GI.LinearRingTrait, g1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, g2,
) = _line_curve_process(
    g1, g2;
    in_allow = true, on_allow = true, out_allow = false,
    in_require = false, on_require = false, out_require = false,
    closed_line = true,
    closed_curve = false,
)

"""
    coveredby(::GI.LinearRingTrait, g1, ::GI.LinearRingTrait, g2)::Bool

A linear ring is covered by another linear ring if the vertices and edges of the
first linear ring are on the edges/vertices of the second linear ring. Return
true if those conditions are met, else false.
"""
coveredby(
    ::GI.LinearRingTrait, g1,
    ::GI.LinearRingTrait, g2,
) = _line_curve_process(
    g1, g2;
    in_allow = true, on_allow = true, out_allow = false,
    in_require = false, on_require = false, out_require = false,
    closed_line = true,
    closed_curve = true,
)

"""
    coveredby(::GI.LinearRingTrait, g1, ::GI.PolygonTrait, g2)::Bool

A linear ring is coveredby a polygon if the vertices and edges of the linear
ring are either in the polygon interior or on the polygon edges. This includes
edges of holes. Return true if those conditions are met, else false.
"""
coveredby(
    ::GI.LinearRingTrait, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    g1, g2;
    in_allow =  true, on_allow = true, out_allow = false,
    in_require = true, on_require = false, out_require = false,
    closed_line = true,
)

