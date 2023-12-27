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
signed area might not be the area. This is why signed area is only implemented
for polygons.
=#

"""
    area(geom)::Real

Returns the area of the geometry. This is computed slighly differently for
different geometries:
    - The area of a point is always zero.
    - The area of a curve is always zero.
    - The area of a polygon is the absolute value of the signed area.
    - The area multi-polygon is the sum of the areas of all of the sub-polygons.
"""
area(geom) = area(GI.trait(geom), geom)

"""
    signed_area(geom)::Real

Returns the signed area of the geometry, based on winding order. This is
computed slighly differently for different geometries:
    - The signed area of a point is always zero.
    - The signed area of a curve is always zero.
    - The signed area of a polygon is computed with the shoelace formula and is
    positive if the polygon coordinates wind clockwise and negative if
    counterclockwise.
    - You cannot compute the signed area of a multipolygon as it doesn't have a
    meaning as each sub-polygon could have a different winding order.
"""
signed_area(geom) = signed_area(GI.trait(geom), geom)

# Points
area(::GI.PointTrait, point) = GI.isempty(point) ?
    0 : zero(typeof(GI.x(point)))

signed_area(trait::GI.PointTrait, point) = area(trait, point)

# MultiPoints
function area(::GI.MultiPointTrait, multipoint)
    GI.isempty(multipoint) && return 0
    np = GI.npoint(multipoint)
    np == 0 && return 0
    return zero(typeof(GI.x(GI.getpoint(multipoint, np))))
end

signed_area(trait::GI.MultiPointTrait, multipoint) = area(trait, multipoint)

# Curves
function area(::CT, curve) where CT <: GI.AbstractCurveTrait
    GI.isempty(curve) && return 0
    np = GI.npoint(curve)
    np == 0 && return 0
    return zero(typeof(GI.x(GI.getpoint(curve, np))))
end

signed_area(trait::CT, curve) where CT <: GI.AbstractCurveTrait =
    area(trait, curve)

# MultiCurves
function area(::MCT, multicurve) where MCT <: GI.AbstractMultiCurveTrait
    GI.isempty(multicurve) && return 0
    ng = GI.ngeom(multicurve)
    ng == 0 && return 0
    np = GI.npoint(GI.getgeom(multicurve, ng))
    np == 0 && return 0
    return zero(typeof(GI.x(GI.getpoint(GI.getgeom(multicurve, ng), np))))
end

signed_area(trait::MCT, curve) where MCT <: GI.AbstractMultiCurveTrait =
    area(trait, curve)

# Polygons
area(trait::GI.PolygonTrait, geom) = abs(signed_area(trait, geom))

function signed_area(::GI.PolygonTrait, poly)
    GI.isempty(poly) && return 0
    s_area = _signed_area(GI.getexterior(poly))
    area = abs(s_area)
    area == 0 && return area
    # Remove hole areas from total
    for hole in GI.gethole(poly)
        area -= abs(_signed_area(hole))
    end
    # Winding of exterior ring determines sign
    return area * sign(s_area)
end

# MultiPolygons
area(::GI.MultiPolygonTrait, multipoly) =
    sum((area(poly) for poly in GI.getpolygon(multipoly)), init = 0)

# GeometryCollections
area(::GI.GeometryCollectionTrait, collection) = 
    sum((area(geom) for geom in GI.getgeom(collection)), init = 0)
#=
Helper function:

Calculates the signed area of a given curve. This is equivalent to integrating
to find the area under the curve. Even if curve isn't explicitly closed by
repeating the first point at the end of the coordinates, curve is still assumed
to be closed.
=#
function _signed_area(geom)
    # Close curve, even if last point isn't explicitly repeated 
    np = GI.npoint(geom)
    np == 0 && return 0
    first_last_equal = equals(GI.getpoint(geom, 1), GI.getpoint(geom, np))
    np -= first_last_equal ? 1 : 0 
    # Integrate the area under the curve
    p1 = GI.getpoint(geom, np)
    T = typeof(GI.x(p1))
    area = zero(T)
    for i in 1:np
        p2 = GI.getpoint(geom, i)
        # Accumulate the area into `area`
        area += GI.x(p1) * GI.y(p2) - GI.y(p1) * GI.x(p2)
        p1 = p2
    end
    return area / 2
end