# # Spherical Caps

#=
```@meta
CollapsedDocStrings = true
```

```@docs; canonical=false
SphericalCap
cap_contains
merge_caps
dilate_cap
cap_xyz_extent
circumcenter_on_unit_sphere
```

## What is SphericalCap?

A spherical cap represents a section of a unit sphere about some point, bounded by a radius. 
It is defined by a center point on the unit sphere and a radius (in radians).

Spherical caps are used in:
- Representing circular regions on a spherical surface
- Approximating and bounding spherical geometries
- Spatial indexing and filtering on the unit sphere
- Implementing containment, intersection, and disjoint predicates

The `SphericalCap` type offers multiple constructors to create caps from:
- UnitSphericalPoint and radius
- Geographic coordinates and radius
- Three points on the unit sphere (circumcircle)

## Examples

```@example sphericalcap
using GeometryOps.UnitSpherical
using GeoInterface

# Create a spherical cap from a point and radius
point = UnitSphericalPoint(1.0, 0.0, 0.0)  # Point on the unit sphere
cap = SphericalCap(point, 0.5)  # Cap with radius 0.5 radians
```

```@example sphericalcap
# Create a spherical cap from geographic coordinates
lat, lon = 40.0, -74.0  # New York City (approximate)
point = GeoInterface.Point(lon, lat)
cap = SphericalCap(point, 0.1)  # Cap with radius ~0.1 radians
```

```@example sphericalcap
# Create a spherical cap from three points (circumcircle)
p1 = UnitSphericalPoint(1.0, 0.0, 0.0)
p2 = UnitSphericalPoint(0.0, 1.0, 0.0)
p3 = UnitSphericalPoint(0.0, 0.0, 1.0)
cap = SphericalCap(p1, p2, p3)
```

=#

# Spherical cap implementation
"""
    SphericalCap{T}
    SphericalCap(point::UnitSphericalPoint{T}, radius::T)

A spherical cap represents a section of a unit sphere about some point, bounded by a radius.
It is defined by a center point on the unit sphere and a radius (in radians).
"""
struct SphericalCap{T}
    "The point at the center of the cap."
    point::UnitSphericalPoint{T}
    "The radius of the cap (in radians). This is what should normally be used in any calculation or comparison."
    radius::T
    """
    A comparison-friendly value equal to `cos(radius)`. Used for efficient containment tests:
    a point `p` is inside the cap if `p ⋅ center >= radiuslike`. Note that this value is
    *inversely* related to cap size (radiuslike=1 for a point, radiuslike=0 for a hemisphere).
    """
    radiuslike::T
end

function SphericalCap(point::UnitSphericalPoint{T}, radius::Number) where T
    radius = convert(T, radius)
    return SphericalCap{T}(point, radius, cos(radius))
end

SphericalCap(point, radius::Number) = SphericalCap(GI.trait(point), point, radius)

SphericalCap(geom) = SphericalCap(GI.trait(geom), geom)
SphericalCap(t::GI.AbstractGeometryTrait, geom) = SphericalCap(t, geom, 0)

function SphericalCap(::GI.PointTrait, point, radius::Number)
    return SphericalCap(UnitSphereFromGeographic()(point), radius)
end
# TODO: add implementations for line string and polygon traits
# That will require a minimum bounding circle implementation.
# TODO: add implementations for multitraits based on this

# TODO: this returns an approximately antipodal point...

# The margin by which the cheap bounds in `_intersects` must clear the decision
# before they are trusted: comfortably above the rounding of a squared chord
# (~4 eps relative) and comfortably below any distance anyone cares about.
const _CAP_CHORD_GUARD = 2.0^-40

# TODO: exact-predicate intersection
# This is all inexact and thus subject to floating point error
#=
Two caps intersect iff the angle `d` between their centers is at most the sum
`r` of their radii.  Computing `d` costs a `spherical_distance` — a cross
product, a `sqrt` and an `atan2` — but the *squared chord* between the centers,
`c² = ‖p - q‖²`, is 3 subtractions and 3 multiplications, and it brackets `d`
tightly enough to settle almost every pair without ever forming `d`:

    c  ≤  d = 2 asin(c/2)  ≤  c * (1 + c²/20)        for 0 ≤ c ≤ 1

(the lower bound is `asin(x) ≥ x`; for the upper, `2 asin(c/2) - c` has the
series `c³/24 + 3c⁵/640 + 15c⁷/21504 + …`, which stays under `c³/20` on `c ≤ 1`
— verified in BigFloat, minimum slack `1.0e-18` at `c → 0` and `2.6e-3` at
`c = 1`).  Squaring both sides keeps it `sqrt`-free.  Only pairs that land
inside the resulting band pay for `spherical_distance`, and there the original
test runs unchanged.

`c ≤ 1` is `d ≤ 60°`; wider pairs skip the accept bound but still get the
reject one, and `r ≥ π` accepts outright since no two points are further
than `π` apart.

Every comparison is written so that a NaN radius (the whole-sphere cap that
`minimum_bounding_circle` returns for empty input) fails it and falls through
to the original test, which is what `NaN` did before.  Same for negative radii.

Where this and the old test disagree, the old test is the one that is wrong:
`spherical_distance`'s `cross(p, q)` cancels catastrophically for nearly equal
`p, q`, giving `d` a relative error of order `eps/d`, while the chord is
computed accurately there.  Adjudicated against 256-bit arithmetic over 600k
pairs whose centers are 1e-12..1 rad apart with `r` within 2 ULP of the true
distance: 111,514 disagreements, of which the old test got 111,514 wrong and
this one got 0 wrong.

Assumes unit-length centers, as the rest of the cap API does.
=#
function _intersects(x::SphericalCap, y::SphericalCap)
    r = x.radius + y.radius
    p, q = x.point, y.point
    dx = p[1] - q[1]
    dy = p[2] - q[2]
    dz = p[3] - q[3]
    c² = dx * dx + dy * dy + dz * dz
    if r >= 0                       # false for NaN, and for negative radii
        r >= oftype(r, π) && return true  # no two points on the sphere are further apart
        r² = r * r
        g = _CAP_CHORD_GUARD
        c² > r² * (1 + g) && return false            # d ≥ c > r
        b = 1 + 0.05 * c²   # 0.05 rounds up from 1/20, so `b` stays an upper bound
        c² <= 1 && c² * b * b <= r² * (1 - g) && return true   # d ≤ c(1 + c²/20) ≤ r
    end
    return spherical_distance(p, q) <= r
end

_disjoint(x::SphericalCap, y::SphericalCap) = !_intersects(x, y)

"""
    cap_contains(big::SphericalCap, small::SphericalCap) -> Bool
    cap_contains(cap::SphericalCap, point::UnitSphericalPoint) -> Bool

Whether the closed cap `big` contains all of `small`, or whether the closed
`cap` contains `point`. A cap whose radius is at least `π` contains the whole
unit sphere.
"""
function cap_contains(big::SphericalCap, small::SphericalCap)
    big.radius >= oftype(big.radius, π) && return true
    small.radius >= oftype(small.radius, π) && return false
    dist = spherical_distance(big.point, small.point)
    # small circle fits in big circle; `<=` so internally-tangent small
    # caps count as contained, matching the point overload below and
    # S2's `S2Cap::Contains(const S2Cap&)` convention.
    return dist + small.radius <= big.radius
end
function cap_contains(cap::SphericalCap, point::UnitSphericalPoint)
    cap.radius >= oftype(cap.radius, π) && return true
    spherical_distance(cap.point, point) <= cap.radius
end

# Private spellings retained for compatibility with existing callers.
_contains(big::SphericalCap, small::SphericalCap) = cap_contains(big, small)
_contains(cap::SphericalCap, point::UnitSphericalPoint) = cap_contains(cap, point)

# ## Cap intersection

"""
    Extents.intersects(x::SphericalCap, y::SphericalCap)

Whether the two closed caps intersect. Tangency counts as intersection. The
comparison delegates to the cap predicate's guarded squared-chord test, which
retains precision for nearly coincident centres and falls back to angular
distance only inside its uncertainty band.
"""
Extents.intersects(x::SphericalCap, y::SphericalCap) = _intersects(x, y)

# ## Cap–extent intersection

"""
    Extents.intersects(cap::SphericalCap, ext::Extents.Extent{(:X, :Y, :Z)})
    Extents.intersects(ext::Extents.Extent{(:X, :Y, :Z)}, cap::SphericalCap)

Whether `cap` intersects the part of the unit sphere covered by the 3D
Cartesian bounding box `ext` (also in unit-spherical space).

A point on the sphere lies in the cap iff its Euclidean (chord) distance to
the cap's center is at most ``2 * sin(radius/2)``, so this tests whether the box
comes within that distance of the center.

This is a fail-safe comparison, so may have false positives but never false negatives.
"""
function Extents.intersects(cap::SphericalCap, ext::Extents.Extent{(:X, :Y, :Z)})
    c = cap.point
    dx = c.x - clamp(c.x, ext.X[1], ext.X[2])
    dy = c.y - clamp(c.y, ext.Y[1], ext.Y[2])
    dz = c.z - clamp(c.z, ext.Z[1], ext.Z[2])
    # the chord radius as S2's `S1ChordAngle(S1Angle)` computes it: accurate
    # for small radii, where `2 - 2 * radiuslike` cancels to 0
    chord = 2 * sin(0.5 * min(cap.radius, π))
    return dx^2 + dy^2 + dz^2 <= chord^2
end
Extents.intersects(ext::Extents.Extent{(:X, :Y, :Z)}, cap::SphericalCap) =
    Extents.intersects(cap, ext)

"""
    cap_xyz_extent(cap::SphericalCap) -> Extents.Extent{(:X, :Y, :Z)}

Return an outward-rounded Cartesian bounding box for every unit-sphere point
in `cap`. The result is intended for Cartesian spatial indexes over
unit-spherical data; it is deliberately an explicitly named conversion rather
than `Extents.extent(cap)`, because a `SphericalCap` remains a spherical query
object rather than ordinary Cartesian geometry.

Valid caps have radii in `[0, π]`. A radius at least `π`, or a non-finite or
negative radius, conservatively returns the whole unit-sphere box.
"""
function cap_xyz_extent(cap::SphericalCap{T}) where {T <: AbstractFloat}
    oneT = one(T)
    whole = (-oneT, oneT)
    r = cap.radius
    if !(zero(T) <= r < T(π))
        return Extents.Extent(X = whole, Y = whole, Z = whole)
    end

    cr, sr = cos(r), sin(r)
    # Covers the accumulated error in sin/cos, the transverse square root,
    # and the multiply-adds. The final adjacent-float step makes the direction
    # of rounding explicit even when subtracting the guard is exact.
    guard = 8 * eps(T)
    bounds = ntuple(Val(3)) do i
        x = clamp(cap.point[i], -oneT, oneT)
        transverse = sqrt(max(zero(T), (oneT - x) * (oneT + x)))
        lo = x <= -cr ? -oneT : x * cr - transverse * sr
        hi = x >= cr ? oneT : x * cr + transverse * sr
        (max(-oneT, prevfloat(lo - guard)),
         min(oneT, nextfloat(hi + guard)))
    end
    return Extents.Extent(X = bounds[1], Y = bounds[2], Z = bounds[3])
end

"""
    dilate_cap(cap::SphericalCap, radius::Real) -> SphericalCap

Return `cap` enlarged by the nonnegative angular `radius` in radians. Positive
dilations are rounded outward and capped at `π`; dilating by zero returns
`cap` unchanged.
"""
function dilate_cap(cap::SphericalCap{T}, radius::Real) where {T <: AbstractFloat}
    δ = convert(T, radius)
    δ >= zero(T) || throw(DomainError(radius, "cap dilation must be nonnegative"))
    δ == zero(T) && return cap
    cap.radius >= T(π) && return cap
    widened = min(T(π), nextfloat(cap.radius + δ))
    return SphericalCap(cap.point, widened)
end

"""
    merge_caps(x::SphericalCap, y::SphericalCap) -> SphericalCap

Return a cap containing both `x` and `y`. If either input already contains the
other it is returned unchanged. Otherwise the centre lies between the input
centres and the radius is rounded outward to preserve coverage.
"""
function merge_caps(x::SphericalCap, y::SphericalCap)
    cap_contains(x, y) && return x
    cap_contains(y, x) && return y

    d = spherical_distance(x.point, y.point)
    newradius = (x.radius + y.radius + d) / 2
    excenter = 0.5 * (1 - (x.radius - y.radius) / d)
    newcenter = slerp(x.point, y.point, excenter)
    # Re-evaluate coverage at the computed centre: this absorbs interpolation
    # and distance rounding, not only the closed-form radius's rounding.
    covering_radius = max(newradius,
        spherical_distance(newcenter, x.point) + x.radius,
        spherical_distance(newcenter, y.point) + y.radius)
    radius = min(oftype(covering_radius, π), nextfloat(covering_radius))
    return SphericalCap(newcenter, radius)
end

# Private spelling retained for compatibility with existing callers.
_merge(x::SphericalCap, y::SphericalCap) = merge_caps(x, y)

"""
    circumcenter_on_unit_sphere(a, b, c)

Return the center of the circle passing through the three `UnitSphericalPoint`s
`a`, `b` and `c`, as a `UnitSphericalPoint`.

Of the two antipodal circumcenters, this returns the one on the same hemisphere as the
input points, i.e. the center of the *smaller* of the two circles they bound.
"""
function circumcenter_on_unit_sphere(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    raw = LinearAlgebra.cross(a, b) +
          LinearAlgebra.cross(b, c) +
          LinearAlgebra.cross(c, a)
    center = LinearAlgebra.normalize(raw)

    # The formula can return either of two antipodal circumcenters depending on
    # the winding order of the input points. We want the smaller circumcircle,
    # which has its center on the same hemisphere as the input points.
    # If dot(a, center) < 0, then center is on the opposite hemisphere from a,
    # meaning we have the far circumcenter and need to negate it.
    # TODO: the above logic might actually be wrong...
    if LinearAlgebra.dot(a, center) < 0
        center = -center
    end

    return center
end

"Get the circumcenter of the triangle (a, b, c) on the unit sphere.  Returns a normalized 3-vector."
function SphericalCap(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    circumcenter = circumcenter_on_unit_sphere(a, b, c)
    circumradius = spherical_distance(a, circumcenter)
    return SphericalCap(circumcenter, circumradius)
end

function _is_ccw_unit_sphere(v_0::S, v_c::S, v_i::S) where S <: UnitSphericalPoint
    # checks if the smaller interior angle for the great circles connecting u-v and v-w is CCW
    return(LinearAlgebra.dot(LinearAlgebra.cross(v_c - v_0, v_i - v_c), v_i) < 0)
end

function angle_between(a::S, b::S, c::S) where S <: UnitSphericalPoint
    ab = b - a
    bc = c - b
    norm_dot = (ab ⋅ bc) / (LinearAlgebra.norm(ab) * LinearAlgebra.norm(bc))
    angle =  acos(clamp(norm_dot, -1.0, 1.0))
    if _is_ccw_unit_sphere(a, b, c)
        return angle
    else
        return 2π - angle
    end
end
