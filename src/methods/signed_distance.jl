# # Signed distance

export signed_distance

#=
## What is signed distance?

Signed distance is the distance of a point to a given geometry. Points within
the geometry have a negative distance and points outside of the geometry have a
positive distance.

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


## Implementation

This is the GeoInterface-compatible implementation. First, we implement a
wrapper method that dispatches to the correct implementation based on the
geometry trait. This is also used in the implementation, since it's a lot less
work!
=#

Base.@propagate_inbounds euclid_distance(x1, y1, x2, y2) =
    sqrt((x2 - x1)^2 + (y2 - y1)^2)

Base.@propagate_inbounds euclid_distance(p1, p2) = euclid_distance(
    GeoInterface.x(p1), GeoInterface.y(p1),
    GeoInterface.x(p2), GeoInterface.y(p2),
)

"""
    signed_distance(geom, x::Real, y::Real)::Float64

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

signed_distance(::GI.PointTrait, p0, ::GI.PointTrait, point) =
    euclid_distance(p0, point)

function signed_distance(::GI.PointTrait, p0, ::GI.MultiPointTrait, mpoint)
    T = typeof(GI.x(p0))
    min_dist = typemax(T)
    for p1 in GI.getpoint(mpoint)
        dist = euclid_distance(p0, p1)
        min_dist = dist < min_dist ? dist : min_dist
    end
    return min_dist
end

signed_distance(::GI.PointTrait, p0, ::GI.LineTrait, line) = 
    _distance(p0, GI.getpoint(line, 1), GI.getpoint(line, 2))

function signed_distance(::GI.PointTrait, p0, ::GI.LineStringTrait, linestring)
    T = typeof(GI.x(p0))
    min_dist = typemax(T)
    p1 = GI.getpoint(linestring, 1)
    for i in 2:GI.npoint(linestring)
        p2 = GI.getpoint(linestring, i)
        dist = _distance(p0, p1, p2)
        min_dist = dist < min_dist ? dist : min_dist
        p1 = p2
    end
    return min_dist
end

function signed_distance(::GI.PointTrait, p0, ::GI.LinearRingTrait, ring)
    # See if linear ring has explicitly repeated last point in coordinates
    np = GI.npoint(ring)
    closed = equals(GI.getpoint(ring, 1), GI.getpoint(ring, np))
    np -= closed ? 1 : 0
    # Find minimum distance
    T = typeof(GI.x(p0))
    min_dist = typemax(T)
    p1 = GI.getpoint(ring, 1)
    for i in 2:np
        p2 = GI.getpoint(ring, i)
        dist = _distance(p0, p1, p2)
        min_dist = dist < min_dist ? dist : min_dist
        p1 = p2
    end
    # Make sure to check closing edge
    min_dist = min(min_dist, _distance(p0, p1, GI.getpoint(ring, 1)))
    return min_dist
end

function signed_distance(::GI.PointTrait, p0, ::GeoInterface.PolygonTrait, poly)
    min_dist = signed_distance(p0, GeoInterface.getexterior(poly))
    @inbounds for hole in GeoInterface.gethole(poly)
        dist = signed_distance(p0, hole)
        min_dist = dist < min_dist ? dist : min_dist
    end
    # should be negative if point is inside polygon
    return GI.contains(poly, p0) ? min_dist : -min_dist
end

function signed_distance(::GI.PointTrait, p0, ::GI.MultiPolygonTrait, mpoly)
    max_min_dist = signed_distance(p0, GI.getpolygon(mpoly, 1))
    for i in 2:GI.npolygon(mpoly)
        dist = signed_distance(p0, GI.getpolygon(mpoly, i))
        max_min_dist = dist > min_dist ? dist : max_min_dist
    end
    return max_min_dist
end

function _distance(p0, p1, p2)
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
        return euclid_distance(x0, y0, xfirst, yfirst)
    end

    c2 = sum(v .* v)
    if c2 <= c1 # p0 is closest to last endpoint
        return euclid_distance(x0, y0, xlast, ylast)
    end

    b2 = c1 / c2  # projection fraction
    return euclid_distance(x0, y0, xfirst + (b2 * v[1]), yfirst + (b2 * v[2]))
end


