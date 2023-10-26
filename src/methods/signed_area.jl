# # Area and signed area

export area, signed_area

#=
## What is area? What is signed area?

Area is the amount of space occupied by a two-dimensional figure. It is always
a positive value. Signed area is simply the integral over the exterior path of
a polygon, minus the sum of integrals over its interior holes. It is signed such
that a clockwise path has a positive area, and a counterclockwise path has a
negative area. The area is the absolute value of the signed area.

To provide an example, consider this rectangle:
```@example rect
using GeometryOps
using GeometryOps.GeometryBasics
using Makie

rect = Polygon([Point(0,0), Point(0,1), Point(1,1), Point(1,0), Point(0, 0)])
f, a, p = poly(rect; axis = (; aspect = DataAspect()))
```
This is clearly a rectangle, etc.  But now let's look at how the points look:
```@example rect
lines!(a, rect; color = 1:length(coordinates(rect))+1)
f
```
The points are ordered in a clockwise fashion, which means that the signed area
is negative.  If we reverse the order of the points, we get a postive area.

## Implementation

This is the GeoInterface-compatible implementation. First, we implement a
wrapper method that dispatches to the correct implementation based on the
geometry trait. This is also used in the implementation, since it's a lot less
work!

Note that area (and signed area) are zero for all points and curves, even
if the curves are closed like with a linear ring. Also note that signed area
really only makes sense for polygons, given with a multipolygon can have several
polygons each with a different orientation and thus the absolute value of the
signed area might not be the area. Caution when using this function!
=#

"""
    area(geom)::Real

Returns the area of the geometry.
"""
area(geom) = area(GI.trait(geom), geom)

"""
    signed_area(geom)::Real

Returns the signed area of the geometry, based on winding order.
"""
signed_area(geom) = signed_area(GI.trait(geom), geom)

"""
    area(::GI.PointTrait, point)::Real

The area of a point is always zero. 
"""
function area(::GI.PointTrait, point)
    T = typeof(GI.x(point))
    return T(0)
end

"""
    signed_area(::GI.PointTrait, point)::Real

The signed area of a point is always zero. 
"""
signed_area(trait::GI.PointTrait, point) = signed_area(trait, point)

"""
    area(::GI.AbstractCurveTrait, curve)::Real

The area of a curve is always zero. 
"""
function area(::CT, curve) where CT <: GI.AbstractCurveTrait
    T = typeof(GI.x(GI.getpoint(curve, 1)))
    return T(0)
end

"""
    signed_area(::GI.AbstractCurveTrait, curve)::Real

The signed area of a curve is always zero. 
"""
signed_area(trait::CT, curve) where CT <: GI.AbstractCurveTrait =
    area(trait, curve)

"""
    area(::GI.PolygonTrait, curve)::Real

Finds the area of a polygon, which is the absolute value of the signed area.
"""
area(trait::GI.PolygonTrait, geom) = abs(signed_area(trait, geom))

"""
    signed_area(::GI.PolygonTrait, curve)::Real

Finds the signed area of a polygon. This is positive if the polygon is clockwise
and negative if it is a counterclockwise path.
"""
function signed_area(::GI.PolygonTrait, poly)
    s_area = _signed_area(GI.getexterior(poly))
    area = abs(s_area)
    for hole in GI.gethole(poly)
        area -= abs(_signed_area(hole))
    end
    return area * sign(s_area)
end

"""
    area(::GI.MultiPolygonTrait, curve)::Real

Finds the area of a multi-polygon, which is the sum of the areas of all of the
sub-polygons.
"""
area(::GI.MultiPolygonTrait, geom) =
    sum((area(poly) for poly in GI.getpolygon(geom)))

"""
    signed_area(::GI.MultiPolygonTrait, curve)::Real

Finds the signed area of a multi-polygon. This value doesn't really have an
inuitive meaning given each sub-polygon can be clockwise or couterclockwise.
"""
signed_area(::GI.MultiPolygonTrait, geom) =
    sum((signed_area(poly) for poly in GI.getpolygon(geom)))

"""
    _signed_area(geom)::Real

Calculates the signed area of a given curve. This is equivalent to integrating
to find the area under the curve.
"""
function _signed_area(geom)
    # Integrate the area under the curve
    point₁ = GI.getpoint(geom, 1)
    point₂ = GI.getpoint(geom, 2)
    area = GI.x(point₁) * GI.y(point₂) - GI.y(point₁) * GI.x(point₂)
    for point in GI.getpoint(geom)
        # Advance the point buffers by 1 point
        point₁ = point₂
        point₂ = point
        # Accumulate the area into `area`
        area += GI.x(point₁) * GI.y(point₂) - GI.y(point₁) * GI.x(point₂)
    end
    area /= 2
    return area
end