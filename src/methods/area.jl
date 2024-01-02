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

const _AREA_TARGETS = Union{GI.PolygonTrait,GI.AbstractCurveTrait,GI.MultiPontTrait,GI.PointTrait}

"""
    area(geom, ::Type{T} = Float64)::T

Returns the area of the geometry. This is computed slighly differently for
different geometries:
    - The area of a point/multipoint is always zero.
    - The area of a curve/multicurve is always zero.
    - The area of a polygon is the absolute value of the signed area.
    - The area multi-polygon is the sum of the areas of all of the sub-polygons.
    - The area of a geometry collection is the sum of the areas of all of the
    sub-geometries. 

Result will be of type T, where T is an optional argument with a default value
of Float64.
"""
function area(geom, ::Type{T} = Float64; threaded=false) where T <: AbstractFloat
    applyreduce(+, _AREA_TARGETS, geom; threaded) do g
        _area(T, GI.trait(g), g)
    end
end


"""
    signed_area(geom, ::Type{T} = Float64)::T

Returns the signed area of the geometry, based on winding order. This is
computed slighly differently for different geometries:
    - The signed area of a point is always zero.
    - The signed area of a curve is always zero.
    - The signed area of a polygon is computed with the shoelace formula and is
    positive if the polygon coordinates wind clockwise and negative if
    counterclockwise.
    - You cannot compute the signed area of a multipolygon as it doesn't have a
    meaning as each sub-polygon could have a different winding order.

Result will be of type T, where T is an optional argument with a default value
of Float64.
"""
function signed_area(geom, ::Type{T} = Float64) where T <: AbstractFloat
    applyreduce(+, _AREA_TARGETS, geom) do g
        _signed_area(T, GI.trait(g), g)
    end
end

# Points, MultiPoints, Curves, MultiCurves
_area(::Type{T}, ::GI.AbstractGeometryTrait, geom) where T = zero(T)

_signed_area(::Type{T}, ::GI.AbstractGeometryTrait, geom) where T = zero(T)

# Polygons
_area(::Type{T}, trait::GI.PolygonTrait, poly) where T =
    abs(_signed_area(T, trait, poly))

function _signed_area(::Type{T}, ::GI.PolygonTrait, poly) where T
    GI.isempty(poly) && return zero(T)
    s_area = _signed_area(T, GI.getexterior(poly))
    area = abs(s_area)
    area == 0 && return area
    # Remove hole areas from total
    for hole in GI.gethole(poly)
        area -= abs(_signed_area(T, hole))
    end
    # Winding of exterior ring determines sign
    return area * sign(s_area)
end

#=
Helper function:

Calculates the signed area of a given curve. This is equivalent to integrating
to find the area under the curve. Even if curve isn't explicitly closed by
repeating the first point at the end of the coordinates, curve is still assumed
to be closed.
=#
function _signed_area(::Type{T}, geom) where T
    area = zero(T)
    np = GI.npoint(geom)
    np == 0 && return area

    first = true
    local pfirst, p1
    # Integrate the area under the curve
    for p2 in GI.getpoint(geom)
        # Skip the first and do it later 
        # This lets us work within one iteration over geom, 
        # which means on C call when using points from external libraries.
        if first
            p1 = pfirst = p2
            first = false
            continue
        end
        # Accumulate the area into `area`
        area += GI.x(p1) * GI.y(p2) - GI.y(p1) * GI.x(p2)
        p1 = p2
    end
    # Complete the last edge.
    # If the first and last where the same this will be zero
    p2 = pfirst
    area += GI.x(p1) * GI.y(p2) - GI.y(p1) * GI.x(p2)
    return T(area / 2)
end
