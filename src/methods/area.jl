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

- `AutoManifold()` (default): When the Proj extension is loaded, recognized
  geographic CRSs use a geodesic calculation on the CRS ellipsoid and
  recognized projected CRSs use native-unit `Planar()` calculations. Without
  Proj, geographic geometries use degree-based `Spherical()` calculations,
  while projected and unknown geometries use native-unit `Planar()` calculations.
  The manifold is selected once from the top-level input's CRS and applies to
  all geometries contained in that input.
- `Planar()`: Uses the shoelace formula in native coordinate units squared,
  regardless of CRS.
- `Spherical()`: Uses Girard's theorem for spherical polygons. Coordinates
   are interpreted as (longitude, latitude) in degrees. Returns area in
   square units of the sphere's radius (default: Earth's mean radius in meters).
- `Geodesic()`: Uses geodesic calculations (requires Proj extension).

Projected map area is a planar grid measurement and can differ from surface
area, particularly for projections with distortion.

## Examples

```julia
import GeometryOps as GO
import GeoInterface as GI

# CRS-free planar area (the AutoManifold default)
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
    area(AutoManifold(), geom, T; threaded, kwargs...)
end

function area(::AutoManifold, geom, ::Type{T} = Float64; threaded=false, kwargs...) where T <: AbstractFloat
    _area_auto(GI.crstrait(geom), GI.crs(geom), geom, T; threaded, kwargs...)
end

function _area_auto(trait::GI.AbstractCRSTrait, crs, geom, ::Type{T}; threaded=false, kwargs...) where T
    isnothing(crs) && return area(Planar(), geom, T; threaded, kwargs...)
    _area_auto_with_crs(trait, crs, geom, T; threaded, kwargs...)
end

# Without Proj, GeoInterface's geographic trait is the only signal to interpret coordinates as lon/lat.
_area_auto_with_crs(::GI.AbstractGeographicTrait, crs, geom, ::Type{T}; threaded=false, kwargs...) where T =
    area(Spherical(), geom, T; threaded, kwargs...)

# A projected CRS has map-plane coordinates; preserve their native square units without Proj.
_area_auto_with_crs(::GI.AbstractProjectedTrait, crs, geom, ::Type{T}; threaded=false, kwargs...) where T =
    area(Planar(), geom, T; threaded, kwargs...)

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
# This is the Van Oosterom–Strackee tangent-half-angle form, kept natively SIGNED so that
# reflex triangles in a fan triangulation subtract (concave rings) instead of adding.
# `numerator` is the raw signed triple product a ⋅ (b × c); `denominator` is 1 + a⋅b + b⋅c + c⋅a.
# The `(b - a, c - a)` / `(b + a, c + a)` rewrites are algebraically identical to the commented
# lines below but better conditioned for tiny triangles (Eriksson 1990), so small-polygon
# accuracy is preserved. `atan2` also branches correctly when `denominator ≤ 0` (large fan
# triangles whose true |area| exceeds π), which `atan(abs(t))` cannot represent.
function _spherical_triangle_area(::Eriksson, a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    #numerator = dot(a, cross(b, c))
    #denominator = 1 + dot(b,c) + dot(c, a) + dot(a, b)
    numerator = dot(a, (cross(b - a, c - a)))
    denominator = dot(b + a, c + a)
    return 2 * atan(numerator, denominator)
end

# ## Fan triangulation of a ring
#
# A ring's signed area is the sum of the signed excesses of the triangles `(apex, vᵢ, vᵢ₊₁)`
# over its edges, with the apex one of its own vertices. Two things about the apex matter.
#
# The apex must not have a vertex near its antipode. Both triangle formulas above lose
# their conditioning there — for `Eriksson`, a vertex `v` with `v ≈ -apex` makes
# `denominator = (v + apex) ⋅ (…)` cancel to rounding noise, and the two triangles sharing
# that fan side come out with the wrong magnitude (or, past ~1e-9, the wrong quadrant).
# With `v = -apex` exactly the fan side is not even a well-defined arc. So the fan starts
# from the first vertex, and moves to another one only if the ring puts a vertex within
# `_FAN_APEX_ANTIPODE_TOL` (chord length, about 0.06°) of the first vertex's antipode: that
# can only happen to a ring spanning nearly 180°, so ordinary rings sum exactly as before.
#
# The apex also fixes which of the two regions a ring bounds is "outside": the fan never
# covers the apex's antipode, so the sum is `A` when that antipode is outside the ring and
# `A - 4π` when it is inside. Under `oriented = false` the interior is by definition the
# smaller region, and `_fold_to_hemisphere` picks the representative in `[-2π, 2π]`, so the
# result does not depend on which vertex the ring starts from. Under `oriented = true` the
# sum is returned as is: choosing between `A` and `A - 4π` there needs the ring's
# orientation, which this kernel does not see.
const _FAN_APEX_ANTIPODE_TOL = 1e-3
_near_antipodal(a::UnitSphericalPoint, b::UnitSphericalPoint) = dot(a, b) < -1 + _FAN_APEX_ANTIPODE_TOL^2 / 2

function _fold_to_hemisphere(s::T) where T
    s > T(2π) && return s - T(4π)
    s < -T(2π) && return s + T(4π)
    return s
end

# Fan of the first `n` entries of `pts` (any GeoInterface points) from `pts[k]`, in ring
# order. `ok` is false if some vertex disqualifies the apex; with `bail` the sum is then
# abandoned early.
function _spherical_fan_from(method::SphericalTriangleAreaMethod, pts, n, k, ::Type{T}, bail::Bool) where T
    apex = UnitSphericalPoint(GI.PointTrait(), pts[k])
    ok = true
    area = zero(T)
    b = UnitSphericalPoint(GI.PointTrait(), pts[mod1(k + 1, n)])
    for j in 1:(n - 2)
        a = b
        b = UnitSphericalPoint(GI.PointTrait(), pts[mod1(k + j + 1, n)])
        if _near_antipodal(apex, a)
            ok = false
            bail && return (area, false)
        end
        area += _spherical_triangle_area(method, apex, a, b)
    end
    _near_antipodal(apex, b) && (ok = false)
    return (area, ok)
end

# Signed unit-sphere area of the ring `pts[1:n]` (open: no repeated closing point).
function _spherical_fan_area(method::SphericalTriangleAreaMethod, oriented::Bool, pts, n, ::Type{T}) where T
    area, ok = _spherical_fan_from(method, pts, n, 1, T, false)
    if !ok
        #-- Candidates spread around the ring: a band symmetric about the equator is
        #-- antipodally symmetric along a stretch, so consecutive vertices would all fail.
        #-- A ring where every candidate fails is essentially a great circle, and keeps the
        #-- first-vertex sum.
        for j in (4, 2, 6, 1, 3, 5, 7)
            k = 1 + (j * n) ÷ 8
            k == 1 && continue
            alt, ok = _spherical_fan_from(method, pts, n, k, T, true)
            if ok
                area = alt
                break
            end
        end
    end
    return oriented ? area : _fold_to_hemisphere(area)
end

# Compute signed area of a ring given as a geometry
function _naive_triangulated_spherical_ring_area(alg::NaiveTriangulatedSphericalArea, trait::GI.AbstractCurveTrait, ring, T)
    GI.npoint(trait, ring) < 3 && return zero(T)
    points = collect(Iterators.map(p -> UnitSphericalPoint(GI.PointTrait(), p), GI.getpoint(trait, ring)))
    n = length(points)
    # Skip closing point if it matches first
    points[n] ≈ points[1] && (n -= 1)
    n < 3 && return zero(T)
    return _spherical_fan_area(alg.method, manifold(alg).oriented, points, n, T)
end
# Dispatch area(::Spherical, ...) to use NaiveTriangulatedSphericalArea with Eriksson's formula for triangles
function area(m::Spherical, geom, ::Type{T} = Float64; threaded=false, kwargs...) where T <: AbstractFloat
    area(NaiveTriangulatedSphericalArea(m), geom, T; threaded, kwargs...)
end

# Compute the area of a single polygon (exterior minus holes) on the unit sphere.
# These must be top-level functions: a multi-method local function captured by a
# closure gets lowered into a `Core.Box`, making `area` infer as `Any`.
# See https://github.com/JuliaGeo/GeometryOps.jl/issues/407.
function _naive_triangulated_spherical_polygon_area(alg::NaiveTriangulatedSphericalArea, ::Type{T}, ::GI.PolygonTrait, poly) where T
    GI.isempty(poly) && return zero(T)
    ext = GI.getexterior(poly)
    ext_area = abs(_naive_triangulated_spherical_ring_area(alg, GI.trait(ext), ext, T))
    for hole in GI.gethole(poly)
        hole_trait = GI.trait(hole)
        ext_area -= abs(_naive_triangulated_spherical_ring_area(alg, hole_trait, hole, T))
    end
    return ext_area
end
_naive_triangulated_spherical_polygon_area(::NaiveTriangulatedSphericalArea, ::Type{T}, ::GI.PointTrait, point) where T = zero(T)

# ## Ring area over a plain vector of points
#
# `intersection_area` measures the rings a clipping engine produced without wrapping them
# in a geometry first, so these take a vector of points. `closed` says whether that vector
# repeats its first point at the end — the caller knows, and it cannot be re-derived here:
# the two ends of a sliver ring are legitimately close together, and no tolerance can tell
# that apart from a closing point.
#
# They must stay term-for-term equivalent to `_signed_area` and
# `_naive_triangulated_spherical_ring_area` above: that equality is exactly what makes
# `intersection_area(alg, a, b)` and `area(manifold(alg), intersection(alg, a, b))` agree
# to the last bit. Change a formula there, change it here.

# Shoelace, wrapping from the last point to the first. A repeated closing point contributes
# a zero term, so closed and open rings both work.
function _ring_area(::Planar, pts::AbstractVector, ::Type{T}; closed::Bool = true) where T
    n = length(pts)
    n < 3 && return zero(T)
    area = zero(T)
    for i in 1:n
        area += _area_component(pts[i], pts[mod1(i + 1, n)])
    end
    return T(area / 2)
end

# Signed unit-sphere area, by the same fan triangulation.
function _ring_area(m::Spherical, pts::AbstractVector, ::Type{T}; closed::Bool = true) where T
    n = length(pts)
    n < 3 && return zero(T)
    #-- Drop the closing point. Only a ring the caller called closed has one: the `≈` is
    #-- `isapprox`'s default `rtol` (~1.5e-8), which on an OPEN ring would swallow the last
    #-- vertex of any sliver whose ends fall within that of each other — halving its area,
    #-- or zeroing it outright once the remaining fan degenerates.
    closed && UnitSphericalPoint(GI.PointTrait(), pts[n]) ≈ UnitSphericalPoint(GI.PointTrait(), pts[1]) && (n -= 1)
    n < 3 && return zero(T)
    return T(_spherical_fan_area(Eriksson(), m.oriented, pts, n, T))
end

# The factor an area on the unit sphere is scaled by to reach the manifold's own units.
_area_scale(::Planar) = 1
_area_scale(m::Spherical) = m.radius^2

# Main implementation for NaiveTriangulatedSphericalArea
function area(alg::NaiveTriangulatedSphericalArea, geom, ::Type{T} = Float64; threaded=false, kwargs...) where T <: AbstractFloat
    unit_area = applyreduce(
        WithTrait((trait, g) -> _naive_triangulated_spherical_polygon_area(alg, T, trait, g)),
        +,
        TraitTarget{Union{GI.PolygonTrait, GI.PointTrait}}(),
        geom;
        threaded,
        init=zero(T),
        kwargs...
    )
    return T(unit_area * manifold(alg).radius^2)
end
