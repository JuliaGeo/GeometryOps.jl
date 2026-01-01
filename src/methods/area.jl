# # Area

export area, signed_area

#=
## What is area? What is signed area?

Area is the amount of space occupied by a two-dimensional figure. It is always a positive
value. Signed area is simply the integral over the exterior path of a polygon, minus the sum
of integrals over its interior holes. It is signed such that a clockwise path has a positive
area, and a counterclockwise path has a negative area. The area is the absolute value of the
signed area.

To provide an example, consider this rectangle:
```@example rect
import GeometryOps as GO
import GeoInterface as GI
using Makie
using CairoMakie

rect = GI.Polygon([[(0,0), (0,1), (1,1), (1,0), (0, 0)]])
f, a, p = poly(collect(GI.getpoint(rect)); axis = (; aspect = DataAspect()))
```
This is clearly a rectangle, etc.  But now let's look at how the points look:
```@example rect
lines!(
    collect(GI.getpoint(rect));
    color = 1:GI.npoint(rect), linewidth = 10.0)
f
```
The points are ordered in a counterclockwise fashion, which means that the signed area
is negative.  If we reverse the order of the points, we get a positive area.
```@example rect
GO.signed_area(rect)  # -1.0
```

## Implementation

This is the GeoInterface-compatible implementation. First, we implement a wrapper method
that dispatches to the correct implementation based on the geometry trait. This is also used
in the implementation, since it's a lot less work!

Note that area and signed area are zero for all points and curves, even if the
curves are closed like with a linear ring. Also note that signed area really only makes
sense for polygons, given with a multipolygon can have several polygons each with a
different orientation and thus the absolute value of the signed area might not be the area.
This is why signed area is only implemented for polygons.
=#

# Targets for applys functions
const _AREA_TARGETS = TraitTarget{Union{GI.PolygonTrait,GI.AbstractCurveTrait,GI.MultiPointTrait,GI.PointTrait}}()

"""
    area(geom, [T = Float64])::T
    area(manifold::Manifold, geom, [T = Float64])::T
    area(algorithm::Algorithm, geom, [T = Float64])::T

Returns the area of a geometry or collection of geometries.
This is computed slightly differently for different geometries:

    - The area of a point/multipoint is always zero.
    - The area of a curve/multicurve is always zero.
    - The area of a polygon is the absolute value of the signed area.
    - The area multi-polygon is the sum of the areas of all of the sub-polygons.
    - The area of a geometry collection, feature collection of array/iterable
        is the sum of the areas of all of the sub-geometries.

## Manifold support

- `Planar()`: Uses the shoelace formula for 2D Cartesian coordinates (default).
- `Spherical()`: Uses Girard's theorem for spherical polygons. Coordinates
   are interpreted as (longitude, latitude) in degrees. Returns area in
   square units of the sphere's radius (default: Earth's mean radius in meters).
- `Geodesic()`: Uses geodesic calculations (requires Proj extension).

## Examples

```julia
import GeometryOps as GO
import GeoInterface as GI

# Planar area (default)
rect = GI.Polygon([[(0,0), (1,0), (1,1), (0,1), (0,0)]])
GO.area(rect)  # 1.0

# Spherical area (1/8 of Earth's surface)
octant = GI.Polygon([[(0.0, 0.0), (90.0, 0.0), (0.0, 90.0), (0.0, 0.0)]])
GO.area(GO.Spherical(), octant)  # ≈ 6.38e13 m²

# Spherical area with custom radius (unit sphere)
GO.area(GO.Spherical(radius=1.0), octant)  # ≈ π/2
```

Result will be of type T, where T is an optional argument with a default value
of Float64.
"""
function area(geom, ::Type{T} = Float64; threaded=false, kwargs...) where T <: AbstractFloat
    area(Planar(), geom, T; threaded, kwargs...)
end

function area(::Planar, geom, ::Type{T} = Float64; threaded=false, kwargs...) where T <: AbstractFloat
    applyreduce(WithTrait((trait, g) -> _area(T, trait, g)), +, _AREA_TARGETS, geom; threaded, init=zero(T), kwargs...)
end

"""
    signed_area(geom, [T = Float64])::T

Returns the signed area of a single geometry, based on winding order. 
This is computed slightly differently for different geometries:

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
signed_area(geom, ::Type{T} = Float64) where T <: AbstractFloat =
    _signed_area(T, GI.trait(geom), geom)

# Points, MultiPoints, Curves, MultiCurves
_area(::Type{T}, ::GI.AbstractGeometryTrait, geom) where T = zero(T)

_signed_area(::Type{T}, ::GI.AbstractGeometryTrait, geom) where T = zero(T)

# LibGEOS treats linear rings as zero area.   I disagree with that but we should probably maintain compatibility...

_area(::Type{T}, tr::GI.LinearRingTrait, geom) where T = 0 # could be abs(_signed_area(T, tr, geom))

_signed_area(::Type{T}, ::GI.LinearRingTrait, geom) where T = 0 # could be _signed_area(T, tr, geom)
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

# One term of the shoelace area formula
_area_component(p1, p2) = GI.x(p1) * GI.y(p2) - GI.y(p1) * GI.x(p2)

#= Calculates the signed area of a given curve. This is equivalent to integrating
to find the area under the curve. Even if curve isn't explicitly closed by
repeating the first point at the end of the coordinates, curve is still assumed
to be closed. =#
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
        area += _area_component(p1, p2)
        p1 = p2
    end
    # Complete the last edge.
    # If the first and last where the same this will be zero
    p2 = pfirst
    area += _area_component(p1, p2)
    return T(area / 2)
end

# ## Spherical Area
# The first implementation here is a naive triangulated implementation.
# The second cut implementation that is planned will use the algorithm that Google's s2 uses
# to get numerically stable triangles from a spherical polygon.

export NaiveTriangulatedSphericalArea

abstract type SphericalTriangleAreaMethod end

struct Girard <: SphericalTriangleAreaMethod end
struct Eriksson <: SphericalTriangleAreaMethod end
struct NaiveTriangulatedSphericalArea{S <: Spherical, T <: SphericalTriangleAreaMethod} <: SingleManifoldAlgorithm{S}
    manifold::S
    method::T
end
NaiveTriangulatedSphericalArea(; radius = Spherical().radius, method = Eriksson()) = NaiveTriangulatedSphericalArea(Spherical(; radius), method)
NaiveTriangulatedSphericalArea(manifold::Spherical) = NaiveTriangulatedSphericalArea(manifold, Eriksson())
GeometryOpsCore.manifold(alg::NaiveTriangulatedSphericalArea) = alg.manifold

using .UnitSpherical: UnitSphericalPoint

# Compute signed area of a spherical triangle on the unit sphere using the half-angle formula.
# Returns the spherical excess E, which equals the area on the unit sphere.
function _spherical_triangle_area(::Girard, p1::UnitSphericalPoint, p2::UnitSphericalPoint, p3::UnitSphericalPoint)
    cross_23 = p2 × p3
    triple = p1 ⋅ cross_23
    d12 = p1 ⋅ p2
    d23 = p2 ⋅ p3
    d31 = p3 ⋅ p1
    denom = 1 + d12 + d23 + d31
    abs(denom) < eps(Float64) && return zero(Float64)
    return 2 * atan(triple, denom)
end

# Using Eriksson's formula for the area of spherical triangles: https://www.jstor.org/stable/2691141
function _spherical_triangle_area(::Eriksson, a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    #t = abs(dot(a, cross(b, c)))
    #t /= 1 + dot(b,c) + dot(c, a) + dot(a, b)
    t = abs(dot(a, (cross(b - a, c - a))) / dot(b + a, c + a))
    return 2*atan(t)
end



# Compute signed area of a ring using streaming iteration (no allocation)
function _naive_triangulated_spherical_ring_area(method::SphericalTriangleAreaMethod, trait::GI.AbstractCurveTrait, ring, T)
    GI.npoint(trait, ring) < 3 && return zero(T)
    # Get first point and remaining points
    p1_geo, rest = Iterators.peel(GI.getpoint(trait, ring))
    p1 = UnitSphericalPoint(GI.PointTrait(), p1_geo)
    pfirst = p1
    # Collect remaining points, converting to unit sphere
    points = collect(Iterators.map(p -> UnitSphericalPoint(GI.PointTrait(), p), rest))
    isempty(points) && return zero(T)
    # Skip closing point if it matches first
    if points[end] ≈ pfirst
        pop!(points)
    end
    length(points) < 2 && return zero(T)
    # Triangulate from first vertex
    area = zero(T)
    for i in 1:(length(points)-1)
        area += _spherical_triangle_area(method, pfirst, points[i], points[i+1])
    end
    return area
end
# Dispatch area(::Spherical, ...) to use NaiveTriangulatedSphericalArea with Eriksson's formula for triangles
function area(m::Spherical, geom, ::Type{T} = Float64; threaded=false, kwargs...) where T <: AbstractFloat
    area(NaiveTriangulatedSphericalArea(m), geom, T; threaded, kwargs...)
end

# Main implementation for NaiveTriangulatedSphericalArea
function area(alg::NaiveTriangulatedSphericalArea, geom, ::Type{T} = Float64; threaded=false, kwargs...) where T <: AbstractFloat

    function _polygon_area(trait::GI.PolygonTrait, alg::SphericalTriangleAreaMethod, poly)
        GI.isempty(poly) && return zero(T)
        ext = GI.getexterior(poly)
        ext_area = abs(_naive_triangulated_spherical_ring_area(alg, GI.trait(ext), ext, T))
        for hole in GI.gethole(poly)
            hole_trait = GI.trait(hole)
            ext_area -= abs(_naive_triangulated_spherical_ring_area(alg, hole_trait, hole, T))
        end
        return ext_area
    end
    _polygon_area(::GI.PointTrait, alg, geom) = zero(T)

    unit_area = applyreduce(
        WithTrait((trait, g) -> _polygon_area(trait, alg.method, g)),
        +,
        TraitTarget{Union{GI.PolygonTrait, GI.PointTrait}}(),
        geom;
        threaded,
        init=zero(T),
        kwargs...
    )
    return T(unit_area * manifold(alg).radius^2)
end
