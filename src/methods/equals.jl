# # Equals

export equals

#=
## What is equals?

The equals function checks if two geometries are equal. They are equal if they
share the same set of points and edges to define the same shape.

To provide an example, consider these two lines:
```@example equals
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
We can see that the two lines do not share a commen set of points and edges in
the plot, so they are not equal:
```@example equals
equals(l1, l2)  # returns false
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait. This is also used in the
implementation, since it's a lot less work! 

Note that while we need the same set of points and edges, they don't need to be
provided in the same order for polygons. For for example, we need the same set
points for two multipoints to be equal, but they don't have to be saved in the
same order. The winding order also doesn't have to be the same to represent the
same geometry. This requires checking every point against every other point in
the two geometries we are comparing. Also, some geometries must be "closed" like
polygons and linear rings. These will be assumed to be closed, even if they
don't have a repeated last point explicity written in the coordinates.
Additionally, geometries and multi-geometries can be equal if the multi-geometry
only includes that single geometry.
=#

"""
    equals(geom1, geom2)::Bool

Compare two Geometries return true if they are the same geometry.

## Examples
```jldoctest
import GeometryOps as GO, GeoInterface as GI
poly1 = GI.Polygon([[(0,0), (0,5), (5,5), (5,0), (0,0)]])
poly2 = GI.Polygon([[(0,0), (0,5), (5,5), (5,0), (0,0)]])

GO.equals(poly1, poly2)
# output
true
```
"""
equals(geom_a, geom_b) = equals(
    GI.trait(geom_a), geom_a,
    GI.trait(geom_b), geom_b,
)

"""
    equals(::T, geom_a, ::T, geom_b)::Bool

Two geometries of the same type, which don't have a equals function to dispatch
off of should throw an error.
"""
equals(::T, geom_a, ::T, geom_b) where T = error("Cant compare $T yet")

"""
    equals(trait_a, geom_a, trait_b, geom_b)

Two geometries which are not of the same type cannot be equal so they always
return false.
"""
equals(trait_a, geom_a, trait_b, geom_b) = false

"""
    equals(::GI.PointTrait, p1, ::GI.PointTrait, p2)::Bool

Two points are the same if they have the same x and y (and z if 3D) coordinates.
"""
function equals(::GI.PointTrait, p1, ::GI.PointTrait, p2)
    GI.ncoord(p1) == GI.ncoord(p2) || return false
    GI.x(p1) == GI.x(p2) || return false
    GI.y(p1) == GI.y(p2) || return false
    if GI.is3d(p1)
        GI.z(p1) == GI.z(p2) || return false 
    end
    return true
end

"""
    equals(::GI.PointTrait, p1, ::GI.MultiPointTrait, mp2)::Bool

A point and a multipoint are equal if the multipoint is composed of a single
point that is equivalent to the given point.
"""
function equals(::GI.PointTrait, p1, ::GI.MultiPointTrait, mp2)
    GI.npoint(mp2) == 1 || return false
    return equals(p1, GI.getpoint(mp2, 1))
end

"""
    equals(::GI.MultiPointTrait, mp1, ::GI.PointTrait, p2)::Bool

A point and a multipoint are equal if the multipoint is composed of a single
point that is equivalent to the given point.
"""
equals(trait1::GI.MultiPointTrait, mp1, trait2::GI.PointTrait, p2) =
    equals(trait2, p2, trait1, mp1)

"""
    equals(::GI.MultiPointTrait, mp1, ::GI.MultiPointTrait, mp2)::Bool

Two multipoints are equal if they share the same set of points.
"""
function equals(::GI.MultiPointTrait, mp1, ::GI.MultiPointTrait, mp2)
    GI.npoint(mp1) == GI.npoint(mp2) || return false
    for p1 in GI.getpoint(mp1)
        has_match = false  # if point has a matching point in other multipoint
        for p2 in GI.getpoint(mp2)
            if equals(p1, p2)
                has_match = true
                break
            end
        end
        has_match || return false  # if no matching point, can't be equal
    end
    return true  # all points had a match
end

"""
    _equals_curves(c1, c2, closed_type1, closed_type2)::Bool

Two curves are equal if they share the same set of point, representing the same
geometry. Both curves must must be composed of the same set of points, however,
they do not have to wind in the same direction, or start on the same point to be
equivalent.
Inputs:
    c1 first geometry
    c2 second geometry
    closed_type1::Bool true if c1 is closed by definition (polygon, linear ring)
    closed_type2::Bool true if c2 is closed by definition (polygon, linear ring)
"""
function _equals_curves(c1, c2, closed_type1, closed_type2)
    # Check if both curves are closed or not
    n1 = GI.npoint(c1)
    n2 = GI.npoint(c2)
    c1_repeat_point = GI.getpoint(c1, 1) == GI.getpoint(c1, n1)
    n2 = GI.npoint(c2)
    c2_repeat_point = GI.getpoint(c2, 1) == GI.getpoint(c2, n2)
    closed1 = closed_type1 || c1_repeat_point
    closed2 = closed_type2 || c2_repeat_point
    closed1 == closed2 || return false
    # How many points in each curve
    n1 -= c1_repeat_point ? 1 : 0
    n2 -= c2_repeat_point ? 1 : 0
    n1 == n2 || return false
    n1 == 0 && return true
    # Find offset between curves
    jstart = nothing
    p1 = GI.getpoint(c1, 1)
    for i in 1:n2
        if equals(p1, GI.getpoint(c2, i))
            jstart = i
            break
        end
    end
    # no point matches the first point
    isnothing(jstart) && return false
    # found match for only point
    n1 == 1 && return true
    # if isn't closed and first or last point don't match, not same curve
    !closed_type1 && (jstart != 1 && jstart != n1) && return false
    # Check if curves are going in same direction
    i = 2
    j = jstart + 1
    j -= j > n2 ? n2 : 0
    same_direction = equals(GI.getpoint(c1, i), GI.getpoint(c2, j))
    # if only 2 points, we have already compared both
    n1 == 2 && return same_direction
    # Check all remaining points are the same wrapping around line
    jstep = same_direction ? 1 : -1
    for i in 2:n1
        ip = GI.getpoint(c1, i)
        j = jstart + (i - 1) * jstep
        j += (0 < j <= n2) ? 0 : (n2 * -jstep)
        jp = GI.getpoint(c2, j)
        equals(ip, jp) || return false
    end
    return true
end

"""
    equals(
        ::Union{GI.LineTrait, GI.LineStringTrait}, l1,
        ::Union{GI.LineTrait, GI.LineStringTrait}, l2,
    )::Bool

Two lines/linestrings are equal if they share the same set of points going
along the curve. Note that lines/linestrings aren't closed by defintion.
"""
equals(
    ::Union{GI.LineTrait, GI.LineStringTrait}, l1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, l2,
) = _equals_curves(l1, l2, false, false)

"""
    equals(
        ::Union{GI.LineTrait, GI.LineStringTrait}, l1,
        ::GI.LinearRingTrait, l2,
    )::Bool

A line/linestring and a linear ring are equal if they share the same set of
points going along the curve. Note that lines aren't closed by defintion, but
rings are, so the line must have a repeated last point to be equal
"""
equals(
    ::Union{GI.LineTrait, GI.LineStringTrait}, l1,
    ::GI.LinearRingTrait, l2,
) = _equals_curves(l1, l2, false, true)

"""
    equals(
        ::GI.LinearRingTrait, l1,
        ::Union{GI.LineTrait, GI.LineStringTrait}, l2,
    )::Bool

A linear ring and a line/linestring are equal if they share the same set of
points going along the curve. Note that lines aren't closed by defintion, but
rings are, so the line must have a repeated last point to be equal
"""
equals(
    ::GI.LinearRingTrait, l1,
    ::Union{GI.LineTrait, GI.LineStringTrait}, l2,
) = _equals_curves(l1, l2, true, false)

"""
    equals(
        ::GI.LinearRingTrait, l1,
        ::GI.LinearRingTrait, l2,
    )::Bool

Two linear rings are equal if they share the same set of points going along the
curve. Note that rings are closed by definition, so they can have, but don't
need, a repeated last point to be equal.
"""
equals(
    ::GI.LinearRingTrait, l1,
    ::GI.LinearRingTrait, l2,
) = _equals_curves(l1, l2, true, true)

"""
    equals(::GI.PolygonTrait, geom_a, ::GI.PolygonTrait, geom_b)::Bool

Two polygons are equal if they share the same exterior edge and holes.
"""
function equals(::GI.PolygonTrait, geom_a, ::GI.PolygonTrait, geom_b)
    # Check if exterior is equal
    _equals_curves(
        GI.getexterior(geom_a), GI.getexterior(geom_b),
        true, true,  # linear rings are closed by definition
    ) || return false
    # Check if number of holes are equal
    GI.nhole(geom_a) == GI.nhole(geom_b) || return false
    # Check if holes are equal
    for ihole in GI.gethole(geom_a)
        has_match = false
        for jhole in GI.gethole(geom_b)
            if _equals_curves(
                ihole, jhole,
                true, true,  # linear rings are closed by definition
            )
                has_match = true
                break
            end
        end
        has_match || return false
    end
    return true
end

"""
    equals(::GI.PolygonTrait, geom_a, ::GI.MultiPolygonTrait, geom_b)::Bool

A polygon and a multipolygon are equal if the multipolygon is composed of a
single polygon that is equivalent to the given polygon.
"""
function equals(::GI.PolygonTrait, geom_a, ::MultiPolygonTrait, geom_b)
    GI.npolygon(geom_b) == 1 || return false
    return equals(geom_a, GI.getpolygon(geom_b, 1))
end

"""
    equals(::GI.MultiPolygonTrait, geom_a, ::GI.PolygonTrait, geom_b)::Bool

A polygon and a multipolygon are equal if the multipolygon is composed of a
single polygon that is equivalent to the given polygon.
"""
equals(trait_a::GI.MultiPolygonTrait, geom_a, trait_b::PolygonTrait, geom_b) = 
    equals(trait_b, geom_b, trait_a, geom_a)

"""
    equals(::GI.PolygonTrait, geom_a, ::GI.PolygonTrait, geom_b)::Bool

Two multipolygons are equal if they share the same set of polygons.
"""
function equals(::GI.MultiPolygonTrait, geom_a, ::GI.MultiPolygonTrait, geom_b)
    # Check if same number of polygons
    GI.npolygon(geom_a) == GI.npolygon(geom_b) || return false
    # Check if each polygon has a matching polygon
    for poly_a in GI.getpolygon(geom_a)
        has_match = false
        for poly_b in GI.getpolygon(geom_b)
            if equals(poly_a, poly_b)
                has_match = true
                break
            end
        end
        has_match || return false
    end
    return true
end