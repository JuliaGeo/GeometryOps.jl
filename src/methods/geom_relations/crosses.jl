# # Crossing checks

#=
## What is crosses?

Two geometries cross when they share some, but not all, of their interior
points, and the shared set has a lower dimension than at least one of them. It
is the DE-9IM predicate `crosses`, and which matrix pattern it means depends on
the dimensions of the two arguments:

  - `dim(a) < dim(b)` — `T*T******`: the interior of `a` must meet the interior
    of `b`, and the interior of `a` must also reach the exterior of `b`.
  - `dim(a) > dim(b)` — `T*****T**`: the transpose of the above.
  - both one-dimensional — `0********`: the interiors must meet, and meet only
    in a set of dimension zero. A collinear overlap is one-dimensional, so it
    is not a crossing.
  - any other combination of dimensions — false. Two points cannot cross, and
    neither can two polygons.

A single point never crosses anything: its interior is one point, and one point
cannot be both interior and exterior to the other geometry.

To provide an example, consider these two lines:
```@example crosses
import GeometryOps as GO
import GeoInterface as GI
using Makie
using CairoMakie

l1 = GI.LineString([(0.0, 0.0), (1.0, 0.0)])
l2 = GI.LineString([(0.5, 1.0), (0.5, -1.0)])

f, a, p = lines(GI.getpoint(l1))
lines!(GI.getpoint(l2))
f
```
The two lines meet at a single point that is interior to both of them, so they
cross:
```@example crosses
GO.crosses(l1, l2)  # true
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait.

Each of these calls a method in the geom_geom_processors file, which decides
whether the two geometries meet a set of criteria. `crosses` never constrains
the boundaries of its arguments, so the "allow" flags are permissive and the
work is done by the requirements:

    - the interiors of the two geometries are required to meet
    - the lower-dimensional geometry is required to have a point in the
      exterior of the higher-dimensional one
    - for two curves, a collinear overlap is disallowed instead, since it would
      make the shared interior set one-dimensional

The code for the specific implementations is in the geom_geom_processors file.
=#

#= Interiors must meet, and the lower-dimensional geometry must reach the other
geometry's exterior. Nothing is asked of either boundary. =#
const CROSSES_POLYGON_ALLOWED = (in_allow = true, on_allow = true, out_allow = true)
const CROSSES_POLYGON_REQUIRES = (in_require = true, on_require = false, out_require = true)
const CROSSES_IN_REQUIRES = (in_require = true, on_require = false, out_require = false)
const CROSSES_OUT_REQUIRES = (in_require = false, on_require = false, out_require = true)
const CROSSES_NO_REQUIRES = (in_require = false, on_require = false, out_require = false)
const CROSSES_EXACT = (exact = False(),)

#= Two curves cross when their interiors meet in a zero-dimensional set, so a
collinear overlap is the one interaction that is disallowed. =#
const CROSSES_CURVE_ALLOWED = (over_allow = false, cross_allow = true, on_allow = true, out_allow = true)
#= ...and the companion that allows everything, used to ask only whether the
interiors meet at all. =#
const CROSSES_MEET_ALLOWED = (over_allow = true, cross_allow = true, on_allow = true, out_allow = true)

#= The point processors double as classifiers: allow exactly one of the three
positions and the returned `Bool` says whether the point is in that position. =#
const CROSSES_POINT_IN_ONLY = (in_allow = true, on_allow = false, out_allow = false)
const CROSSES_POINT_OUT_ONLY = (in_allow = false, on_allow = false, out_allow = true)

#= Every geometry trait `crosses` answers for. Deliberately not
`GI.AbstractGeometryTrait`: `GI.RectangleTrait` is one of those, and the extent
forwarding methods in `common.jl` already dispatch on it. =#
const CROSSES_TRAITS = Union{
    GI.PointTrait, GI.MultiPointTrait,
    GI.AbstractCurveTrait, GI.AbstractMultiCurveTrait,
    GI.PolygonTrait, GI.MultiPolygonTrait,
    GI.GeometryCollectionTrait,
}
const CROSSES_LINE_LIKE = Union{GI.AbstractCurveTrait, GI.AbstractMultiCurveTrait}
const CROSSES_AREA_LIKE = Union{GI.PolygonTrait, GI.MultiPolygonTrait}

"""
    crosses([manifold::Manifold], geom1, geom2)::Bool

Return `true` if the two geometries cross: their interiors meet, and neither
geometry is swallowed by the other. Which question that is depends on the
dimensions of the arguments:

  - one geometry lower-dimensional than the other (`T*T******`): the interiors
    must meet, and the lower-dimensional geometry must also reach the other
    geometry's exterior.
  - both one-dimensional (`0********`): the interiors must meet, and meet only
    in a set of dimension zero, so a shared stretch of line is not a crossing.
  - any other pair of dimensions: `false`. Neither two points nor two polygons
    can cross, and a single point never crosses anything.

`manifold` defaults to `Planar()`. `Spherical()` asks the same questions with
great-circle arcs in place of straight segments.

## Examples
```jldoctest
import GeometryOps as GO, GeoInterface as GI

l1 = GI.LineString([(0.0, 0.0), (1.0, 0.0)])
l2 = GI.LineString([(0.5, 1.0), (0.5, -1.0)])

GO.crosses(l1, l2)
# output
true
```
"""
crosses(g1, g2)::Bool = crosses(Planar(), g1, g2)
crosses(::AutoManifold, g1, g2)::Bool = crosses(Planar(), g1, g2)
crosses(m::Manifold, g1, g2)::Bool = _crosses(m, trait(g1), g1, trait(g2), g2)

"""
    crosses(g1)

Return a function that checks if its input crosses `g1`.
This is equivalent to `x -> crosses(x, g1)`.
"""
crosses(g1) = Base.Fix2(crosses, g1)


# # Helpers

#= A curve or multi-curve, seen as its component curves, so that the
multi-geometry methods can reuse the single-curve processors. =#
_crosses_curves(g) = _crosses_curves(GI.trait(g), g)
_crosses_curves(::GI.AbstractCurveTrait, g) = (g,)
_crosses_curves(::GI.AbstractMultiCurveTrait, g) = GI.getgeom(g)

#= As above, for a polygon or multi-polygon. =#
_crosses_polygons(g) = _crosses_polygons(GI.trait(g), g)
_crosses_polygons(::GI.PolygonTrait, g) = (g,)
_crosses_polygons(::GI.MultiPolygonTrait, g) = GI.getgeom(g)

#= The `closed_line`/`closed_curve` flag the curve processors take. A linear
ring is closed whether or not it repeats its first vertex, which the processors
cannot tell from the coordinates alone. =#
_crosses_closed(g) = GI.trait(g) isa GI.LinearRingTrait

#= Is `point` on the interior of at least one of `curves`' component curves?

Not exactly the interior of a multi-curve, which by the mod-2 rule also holds
any vertex shared by an even number of components. Such a vertex is reported
here as a boundary point of each component and so as neither in nor out; the
answer only changes for a multi-curve whose components meet end to end. =#
_crosses_point_in_curve(m::Manifold, point, curves) = any(_crosses_curves(curves)) do curve
    _point_curve_process(m, point, curve; CROSSES_POINT_IN_ONLY..., closed_curve = _crosses_closed(curve))
end

# Is `point` in the exterior of every one of `curves`' component curves?
_crosses_point_out_curve(m::Manifold, point, curves) = all(_crosses_curves(curves)) do curve
    _point_curve_process(m, point, curve; CROSSES_POINT_OUT_ONLY..., closed_curve = _crosses_closed(curve))
end

# Is `point` in the interior of at least one polygon of `polygons`?
_crosses_point_in_polygon(m::Manifold, point, polygons) = any(_crosses_polygons(polygons)) do polygon
    _point_polygon_process(m, point, polygon; CROSSES_POINT_IN_ONLY..., CROSSES_EXACT...)
end

# Is `point` in the exterior of every polygon of `polygons`?
_crosses_point_out_polygon(m::Manifold, point, polygons) = all(_crosses_polygons(polygons)) do polygon
    _point_polygon_process(m, point, polygon; CROSSES_POINT_OUT_ONLY..., CROSSES_EXACT...)
end

@noinline _throw_crosses_unsupported(t1, t2) = throw(ArgumentError(
    "`crosses` has no implementation for $(typeof(t1).name.name) against " *
    "$(typeof(t2).name.name): the lightweight processors work one polygon and " *
    "one curve at a time, and neither the exterior of a multi-polygon nor the " *
    "dimension of a geometry collection is the union of its parts' answers. " *
    "Use `relate_predicate(RelateNG(), pred_crosses(), a, b)` instead."))


# # Points cross nothing

#= A point's interior is a single point, and `crosses` needs the
lower-dimensional interior to reach both the interior and the exterior of the
other geometry. One point cannot do both. =#
_crosses(m::Manifold, ::GI.PointTrait, g1, ::CROSSES_TRAITS, g2) = false
_crosses(
    m::Manifold,
    ::Union{GI.MultiPointTrait, CROSSES_LINE_LIKE, CROSSES_AREA_LIKE}, g1,
    ::GI.PointTrait, g2,
) = false


# # Multipoints cross curves and polygons

# Two zero-dimensional geometries have equal dimension, so they cannot cross.
_crosses(m::Manifold, ::GI.MultiPointTrait, g1, ::GI.MultiPointTrait, g2) = false

#= `T*T******`: at least one of the points must be on the curve's interior, and
at least one must be in its exterior. Points on the curve's boundary are in
neither, so they count for neither half. =#
function _crosses(m::Manifold, ::GI.MultiPointTrait, g1, ::CROSSES_LINE_LIKE, g2)
    in_req_met = false
    out_req_met = false
    for point in GI.getpoint(g1)
        in_req_met || (in_req_met = _crosses_point_in_curve(m, point, g2))
        out_req_met || (out_req_met = _crosses_point_out_curve(m, point, g2))
        in_req_met && out_req_met && return true
    end
    return false
end

#= `T*T******` again, against a polygon or multi-polygon: one point inside, one
point outside, with points on the boundary counting for neither. =#
function _crosses(m::Manifold, ::GI.MultiPointTrait, g1, ::CROSSES_AREA_LIKE, g2)
    in_req_met = false
    out_req_met = false
    for point in GI.getpoint(g1)
        in_req_met || (in_req_met = _crosses_point_in_polygon(m, point, g2))
        out_req_met || (out_req_met = _crosses_point_out_polygon(m, point, g2))
        in_req_met && out_req_met && return true
    end
    return false
end

#= `T*****T**` is the transpose of `T*T******`, so a higher-dimensional
geometry crosses a multipoint exactly when the multipoint crosses it. =#
_crosses(
    m::Manifold,
    trait1::Union{CROSSES_LINE_LIKE, CROSSES_AREA_LIKE}, g1,
    trait2::GI.MultiPointTrait, g2,
) = _crosses(m, trait2, g2, trait1, g1)


# # Curves cross curves

#= `0********`: the two interiors must meet, and only in a zero-dimensional
set. `over_allow = false` rejects the collinear overlaps that would make the
shared set one-dimensional, and `in_require` asks for the meeting. =#
_crosses(
    m::Manifold,
    ::GI.AbstractCurveTrait, g1,
    ::GI.AbstractCurveTrait, g2,
) = _line_curve_process(
    m, g1, g2;
    CROSSES_CURVE_ALLOWED...,
    CROSSES_IN_REQUIRES...,
    CROSSES_EXACT...,
    closed_line = _crosses_closed(g1),
    closed_curve = _crosses_closed(g2),
)

#= The same question for multi-curves, but asked component pair by component
pair: one overlapping pair anywhere makes the shared interior set
one-dimensional, while the meeting may come from any pair at all. =#
function _crosses(m::Manifold, ::CROSSES_LINE_LIKE, g1, ::CROSSES_LINE_LIKE, g2)
    in_req_met = false
    for c1 in _crosses_curves(g1)
        closed_line = _crosses_closed(c1)
        for c2 in _crosses_curves(g2)
            closed_curve = _crosses_closed(c2)
            # a collinear overlap anywhere rules out `0********` everywhere
            _line_curve_process(
                m, c1, c2;
                CROSSES_CURVE_ALLOWED..., CROSSES_NO_REQUIRES..., CROSSES_EXACT...,
                closed_line, closed_curve,
            ) || return false
            in_req_met && continue
            in_req_met = _line_curve_process(
                m, c1, c2;
                CROSSES_MEET_ALLOWED..., CROSSES_IN_REQUIRES..., CROSSES_EXACT...,
                closed_line, closed_curve,
            )
        end
    end
    return in_req_met
end


# # Curves cross polygons

#= `T*T******`: the line's interior must meet the polygon's interior, and the
line must also have a point outside the polygon. =#
_crosses(
    m::Manifold,
    ::GI.AbstractCurveTrait, g1,
    ::GI.PolygonTrait, g2,
) = _line_polygon_process(
    m, g1, g2;
    CROSSES_POLYGON_ALLOWED...,
    CROSSES_POLYGON_REQUIRES...,
    CROSSES_EXACT...,
    closed_line = _crosses_closed(g1),
)

#= For a multi-curve the two halves of `T*T******` can be met by different
components — one component inside the polygon, another outside it — so they are
required separately rather than in a single call. =#
function _crosses(m::Manifold, ::GI.AbstractMultiCurveTrait, g1, ::GI.PolygonTrait, g2)
    in_req_met = false
    out_req_met = false
    for curve in GI.getgeom(g1)
        closed_line = _crosses_closed(curve)
        in_req_met || (in_req_met = _line_polygon_process(
            m, curve, g2;
            CROSSES_POLYGON_ALLOWED..., CROSSES_IN_REQUIRES..., CROSSES_EXACT..., closed_line,
        ))
        out_req_met || (out_req_met = _line_polygon_process(
            m, curve, g2;
            CROSSES_POLYGON_ALLOWED..., CROSSES_OUT_REQUIRES..., CROSSES_EXACT..., closed_line,
        ))
        in_req_met && out_req_met && return true
    end
    return false
end

# A polygon crosses a curve exactly when the curve crosses the polygon.
_crosses(
    m::Manifold,
    trait1::GI.PolygonTrait, g1,
    trait2::CROSSES_LINE_LIKE, g2,
) = _crosses(m, trait2, g2, trait1, g1)


# # Equal dimensions above one never cross

_crosses(m::Manifold, ::CROSSES_AREA_LIKE, g1, ::CROSSES_AREA_LIKE, g2) = false


# # Cases the lightweight processors cannot express

#= A curve's points outside a multi-polygon are not the points outside any one
of its polygons — a line running from one component into a neighbouring one is
outside each of them and outside neither of them — so `out_require` cannot be
decomposed over the components the way `in_require` can. Rather than answer
with the weaker per-component question, refuse and point at the engine that
does compute the union. =#
_crosses(m::Manifold, t1::CROSSES_LINE_LIKE, g1, t2::GI.MultiPolygonTrait, g2) =
    _throw_crosses_unsupported(t1, t2)
_crosses(m::Manifold, t1::GI.MultiPolygonTrait, g1, t2::CROSSES_LINE_LIKE, g2) =
    _throw_crosses_unsupported(t1, t2)

#= A geometry collection's dimension is the largest of its parts', and which of
the three `crosses` patterns applies depends on it, so a collection cannot be
answered part by part either. =#
_crosses(m::Manifold, t1::GI.GeometryCollectionTrait, g1, t2::CROSSES_TRAITS, g2) =
    _throw_crosses_unsupported(t1, t2)
_crosses(
    m::Manifold,
    t1::Union{GI.MultiPointTrait, CROSSES_LINE_LIKE, CROSSES_AREA_LIKE}, g1,
    t2::GI.GeometryCollectionTrait, g2,
) = _throw_crosses_unsupported(t1, t2)
