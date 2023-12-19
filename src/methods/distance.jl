# # Distance and signed distance

export distance, signed_distance

#=
## What is distance? What is signed distance?

Distance is the distance of a point to another geometry. This is always a
positive number. If a point is inside of geometry, so on a curve or inside of a
polygon, the distance will be zero. Signed distance is mainly used for polygons
and multipolygons. If a point is outside of a geometry, signed distance has the
same value as distance. However, points within the geometry have a negative
distance representing the distance of a point to the closest boundary.
Therefore, for all "non-filled" geometries, like curves, the distance will
either be postitive or 0.

To provide an example, consider this rectangle:
```@example rect
using GeometryOps
using GeometryOps.GeometryBasics
using Makie

rect = Polygon([Point(0,0), Point(0,1), Point(1,1), Point(1,0), Point(0, 0)])
point_in = Point(0.5, 0.5) 
point_out = Point(0.5, 1.5)
f, a, p = poly(rect; axis = (; aspect = DataAspect()))
scatter!(f, point_in)
scatter!(f, point_out)
f
```
This is clearly a rectangle with one point inside and one point outside. The
points are both an equal distance to the polygon. The distance to point_in is
negative while the distance to point_out is positive.
```@example rect
distance(point_in, poly)  # == 0
signed_distance(point_in, poly)  # < 0
signed_distance(point_out, poly)  # > 0
```

## Implementation

This is the GeoInterface-compatible implementation. First, we implement a
wrapper method that dispatches to the correct implementation based on the
geometry trait. This is also used in the implementation, since it's a lot less
work!

Distance and signed distance are only implemented for points to other geometries
right now. This could be extended to include distance from other geometries in
the future.

The distance calculated is the Euclidean distance using the Pythagorean theorem.
Also note that singed_distance only makes sense for "filled-in" shapes, like
polygons, so it isn't implemented for curves.
=#

"""
    distance(g1, g2)::Real

Calculates the  ditance from the geometry `g1` to the `point`. The distance
will always be positive or zero.
"""
distance(point, geom) = distance(
    GI.trait(point), point,
    GI.trait(geom), geom,
)

"""
    distance(::GI.AbstractTrait, geom, ::GI.PointTrait, point)::Real

All distance functions below are defined with the point trait and point as the
first two arguments. If the geometry trait and geometry are first, swap the
argument order.
"""
distance(gtrait::GI.AbstractTrait, geom, ptrait::GI.PointTrait, point) = 
    distance(ptrait, point, gtrait, geom)

"""
    signed_distance(point, geom)::Real

Calculates the signed distance from the geometry `geom` to the point
defined by `(x, y)`.  Points within `geom` have a negative distance,
and points outside of `geom` have a positive distance.

If `geom` is a MultiPolygon, then this function returns the maximum distance 
to any of the polygons in `geom`.
"""
signed_distance(point, geom) = signed_distance(
    GI.trait(point), point,
    GI.trait(geom), geom,
)

"""
    signed_distance(::GI.AbstractTrait, geom, ::GI.PointTrait, point)::Real

All signed distance functions below are defined with the point trait and point
as the first two arguments. If the geometry trait and geometry are first, swap
the argument order.
"""
signed_distance(gtrait::GI.AbstractTrait, geom, ptrait::GI.PointTrait, point) = 
    signed_distance(ptrait, point, gtrait, geom)

"""
    signed_distance(::GI.PointTrait, point, ::GI.AbstractTrait, geom)::Real

The signed distance from a point to a geometry that isn't defined below (polygon
and multipolygon) is simply equal to the distance between those two points
"""
signed_distance(ptrait::GI.PointTrait, point, gtrait::GI.AbstractTrait, geom) =
    distance(ptrait, point, gtrait, geom)

"""
    distance(::GI.PointTrait, point, ::GI.PointTrait, geom)::Real

The distance from a point to a point is just the Euclidean distance between the
points.
"""
distance(::GI.PointTrait, point, ::GI.PointTrait, geom) =
    euclid_distance(point, geom)

"""
    distance(::GI.PointTrait, point, ::GI.MultiPointTrait, geom)::Real

The distance from a point to a multipolygon is the shortest distance from a the
given point to any point within the multipoint object.
"""
function distance(::GI.PointTrait, point, ::GI.MultiPointTrait, geom)
    T = typeof(GI.x(point))
    min_dist = typemax(T)
    for p in GI.getpoint(geom)
        dist = euclid_distance(point, p)
        min_dist = dist < min_dist ? dist : min_dist
    end
    return min_dist
end

"""
    distance(::GI.PointTrait, point, ::GI.LineTrait, geom)::Real

The distance from a point to a line is the minimum distance from the point to
the closest point on the given line.
"""
distance(::GI.PointTrait, point, ::GI.LineTrait, geom) = 
    _distance_line(point, GI.getpoint(geom, 1), GI.getpoint(geom, 2))

"""
    distance(::GI.PointTrait, point, ::GI.LineStringTrait, geom)::Real

The distance from a point to a linestring is the minimum distance from the point
to the closest segment of the linestring.
"""
distance(::GI.PointTrait, point, ::GI.LineStringTrait, geom) =
    _distance_curve(point, geom, close_curve = false)

"""
    distance(::GI.PointTrait, point, ::GI.LinearRingTrait, geom)::Real

The distance from a point to a linear ring is the minimum distance from the
point to the closest segment of the linear ring. Note that the linear ring is
closed by definition, but is not filled in, so the signed distance will always
be positive or zero.
"""
distance(::GI.PointTrait, point, ::GI.LinearRingTrait, geom) =
    _distance_curve(point, geom, close_curve = true)

"""
    distance(::GI.PointTrait, point, ::GI.PolygonTrait, geom)::Real

The distance from a point to a polygon is zero if the point is within the
polygon and otherwise is the minimum distance from the point to an edge of the
polygon. This includes edges created by holes.
"""
function distance(::GI.PointTrait, point, ::GI.PolygonTrait, geom)
    T = typeof(GI.x(point))
    GI.within(point, geom) && return T(0)
    return _distance_polygon(point, geom)
end

"""
    signed_distance(::GI.PointTrait, point, ::GI.PolygonTrait, geom)::Real

The signed distance from a point to a polygon is negative if the point is within
the polygon and is positive otherwise. The value of the distance is the minimum
distance from the point to an edge of the polygon. This includes edges created
by holes.
"""
function signed_distance(::GI.PointTrait, point, ::GI.PolygonTrait, geom)
    min_dist = _distance_polygon(point, geom)
    # should be negative if point is inside polygon
    return GI.within(point, geom) ? -min_dist : min_dist
end

"""
    distance(::GI.PointTrait, point, ::GI.MultiPolygonTrait, geom)

The distance from a point to a multipolygon is zero if the point is within the
multipolygon and otherwise is the minimum distance from the point to the closest
edge of any of the polygons within the multipolygon. This includes edges created
by holes of the polygons as well.
"""
function distance(::GI.PointTrait, point, ::GI.MultiPolygonTrait, geom)
    min_dist = distance(point, GI.getpolygon(geom, 1))
    for i in 2:GI.npolygon(geom)
        min_dist == 0 && return min_dist  # point inside of last polygon checked
        dist = distance(point, GI.getpolygon(geom, i))
        min_dist = dist < min_dist ? dist : min_dist
    end
    return min_dist
end

"""
    signed_distance(::GI.PointTrait, point, ::GI.MultiPolygonTrait, geom)

The signed distance from a point to a mulitpolygon is negative if the point is
within one of the polygons that make up the multipolygon and is positive
otherwise. The value of the distance is the minimum distance from the point to
an edge of the multipolygon. This includes edges created by holes of the
polygons as well.
"""
function signed_distance(::GI.PointTrait, point, ::GI.MultiPolygonTrait, geom)
    min_dist = signed_distance(point, GI.getpolygon(geom, 1))
    for i in 2:GI.npolygon(geom)
        dist = signed_distance(point, GI.getpolygon(geom, i))
        min_dist = dist < min_dist ? dist : min_dist
    end
    return min_dist
end

"""
    euclid_distance(x1::Real, y1::Real, x2::Real, y2::Real)::Real

Returns the Euclidean distance between two points given their x and y values.
"""
Base.@propagate_inbounds _euclid_distance(x1, y1, x2, y2) =
    sqrt((x2 - x1)^2 + (y2 - y1)^2)

"""
    euclid_distance(p1::Point, p2::Point)::Real

Returns the Euclidean distance between two points.
"""
Base.@propagate_inbounds euclid_distance(p1, p2) = _euclid_distance(
    GeoInterface.x(p1), GeoInterface.y(p1),
    GeoInterface.x(p2), GeoInterface.y(p2),
)

"""
    _distance_line(p0, p1, p2)::Real

Returns the minimum distance from point p0 to the line defined by endpoints p1
and p2.
"""
function _distance_line(p0, p1, p2)
    x0, y0 = GeoInterface.x(p0), GeoInterface.y(p0)
    x1, y1 = GeoInterface.x(p1), GeoInterface.y(p1)
    x2, y2 = GeoInterface.x(p2), GeoInterface.y(p2)

    xfirst, yfirst, xlast, ylast = x1 < x2 ?
        (x1, y1, x2, y2) : (x2, y2, x1, y1)
    
    #=
    Vectors from first point to last point (v) and from first point to point of
    interest (w) to find the projection of w onto v to find closest point
    =#
    v = (xlast - xfirst, ylast - yfirst)
    w = (x0 - xfirst, y0 - yfirst)

    c1 = sum(w .* v)
    if c1 <= 0  # p0 is closest to first endpoint
        return _euclid_distance(x0, y0, xfirst, yfirst)
    end

    c2 = sum(v .* v)
    if c2 <= c1 # p0 is closest to last endpoint
        return _euclid_distance(x0, y0, xlast, ylast)
    end

    b2 = c1 / c2  # projection fraction
    return _euclid_distance(x0, y0, xfirst + (b2 * v[1]), yfirst + (b2 * v[2]))
end

"""
    _distance_curve(point, curve; close_curve = false)

Returns the minimum distance from the given point to the given curve. If
close_curve is true, make sure to include the edge from the first to last point
of the curve, even if it isn't explicitly repeated.
"""
function _distance_curve(point, curve; close_curve = false)
    # See if linear ring has explicitly repeated last point in coordinates
    np = GI.npoint(curve)
    first_last_equal = equals(GI.getpoint(curve, 1), GI.getpoint(curve, np))
    close_curve &= first_last_equal
    np -= first_last_equal ? 1 : 0 
    # Find minimum distance
    T = typeof(GI.x(point))
    min_dist = typemax(T)
    p1 = GI.getpoint(curve, close_curve ? np : 1)
    for i in (close_curve ? 1 : 2):np
        p2 = GI.getpoint(curve, i)
        dist = _distance_line(point, p1, p2)
        min_dist = dist < min_dist ? dist : min_dist
        p1 = p2
    end
    return min_dist
end

"""
    _distance_polygon(point, poly)

Returns the minimum distance from the given point to an edge of the given
polygon, including from edges created by holes. Assumes polygon isn't filled and
treats the exterior and each hole as a linear ring. 
"""
function _distance_polygon(point, poly)
    min_dist = _distance_curve(point, GI.getexterior(poly); close_curve = true)
    @inbounds for hole in GI.gethole(poly)
        dist = _distance_curve(point, hole; close_curve = true)
        min_dist = dist < min_dist ? dist : min_dist
    end
    return min_dist
end


