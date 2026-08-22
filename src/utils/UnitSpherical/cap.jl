# # Spherical Caps

#=
```@meta
CollapsedDocStrings = true
```

```@docs; canonical=false
SphericalCap
Extents.contains
Extents.within
Extents.union
Extents.grow
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
`r` of their radii.  We can shortcut by computing the *squared chord* between the centers,
`c² = ‖p - q‖²`, which is 3 subtractions and 3 multiplications. It brackets `d`
tightly enough to settle almost every pair without ever forming `d`:

    c  ≤  d = 2 asin(c/2)  ≤  c * (1 + c²/20)        for 0 ≤ c ≤ 1

(the lower bound is `asin(x) ≥ x`; for the upper, `2 asin(c/2) - c` has the
series `c³/24 + 3c⁵/640 + 15c⁷/21504 + …`, which stays under `c³/20` on `c ≤ 1`
, minimum slack `1.0e-18` at `c → 0` and `2.6e-3` at `c = 1`).
Squaring both sides keeps it `sqrt`-free.

This acts as an adaptive filter, similar to how the predicates work.
Anything which is sufficiently far away gets caught by this, closer things
go to the old comparison for accuracy.

The accept bound requires `c ≤ 1` (`d ≤ 60°`), but the reject bound does not.
Since no two points are more than `π` apart, `r ≥ π` accepts immediately.
NaN radii (used for empty-input whole-sphere caps) and negative radii fail the
filter comparisons and reach the original test, preserving its semantics.

For nearly equal centers, `spherical_distance` suffers catastrophic cancellation
in `cross(p, q)` and relative error of order `eps/d`; the chord remains accurate.
Against 256-bit arithmetic on 600k pairs 1e-12..1 rad apart, with `r` within
2 ULP of the true distance, all 111,514 disagreements were old-test errors;
the chord test had none.

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
    Extents.contains(big::SphericalCap, small::SphericalCap; strict=false)

Whether the closed cap `big` contains all of `small`.
"""
function Extents.contains(big::SphericalCap, small::SphericalCap; strict=false)
    big.radius >= oftype(big.radius, π) && return true
    small.radius >= oftype(small.radius, π) && return false
    dist = spherical_distance(big.point, small.point)
    # small circle fits in big circle; `<=` so internally-tangent small
    # caps count as contained, matching the point overload below and
    # S2's `S2Cap::Contains(const S2Cap&)` convention.
    return dist + small.radius <= big.radius
end
"""
    Extents.within(small::SphericalCap, big::SphericalCap; strict=false)

Whether the closed cap `small` is within `big`.
"""
Extents.within(small::SphericalCap, big::SphericalCap; strict=false) =
    Extents.contains(big, small; strict)

# Point membership stays internal: a point is not an `Extents.Extent`, so the
# Extents DE-9IM containment generic does not describe this mixed relation.
function _contains(cap::SphericalCap, point::UnitSphericalPoint)
    cap.radius >= oftype(cap.radius, π) && return true
    spherical_distance(cap.point, point) <= cap.radius
end

# ## Cap intersection

"""
    Extents.intersects(x::SphericalCap, y::SphericalCap)

Whether the two closed caps intersect, including at tangency.
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
    convert(::Type{<:Extents.Extent}, cap::SphericalCap) -> Extents.Extent{(:X, :Y, :Z)}

Convert `cap` to an outward-rounded Cartesian bounding box containing every
unit-sphere point in the cap. The result is intended for Cartesian spatial
indexes over unit-spherical data.

Extent extraction conventionally means that an object *has* a Cartesian
extent. A `SphericalCap` is instead a spherical query object, so conversion
makes this boundary crossing explicit.

Valid caps have radii in `[0, π]`. A radius at least `π`, or a non-finite or
negative radius, conservatively returns the whole unit-sphere box.

The abstract target `Extents.Extent` and compatible concrete targets are
accepted. Incompatible concrete targets throw an `ArgumentError`.
"""
function Base.convert(target::Type{<:Extents.Extent}, cap::SphericalCap{T}) where {
        T <: AbstractFloat}
    oneT = one(T)
    whole = (-oneT, oneT)
    r = cap.radius
    if !(zero(T) <= r < T(π))
        result = Extents.Extent(X = whole, Y = whole, Z = whole)
        result isa target || throw(ArgumentError(
            "cannot convert SphericalCap{$T} to incompatible extent type $target"))
        return result
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
    result = Extents.Extent(X = bounds[1], Y = bounds[2], Z = bounds[3])
    result isa target || throw(ArgumentError(
        "cannot convert SphericalCap{$T} to incompatible extent type $target"))
    return result
end

"""
    Extents.grow(cap::SphericalCap, factor::Real) -> SphericalCap

Grow the cap's angular diameter by `factor` on each side.
"""
function Extents.grow(cap::SphericalCap{T}, factor::Real) where {T <: AbstractFloat}
    δ = convert(T, factor)
    δ == zero(T) && return cap
    δ > zero(T) && cap.radius >= T(π) && return cap
    grown = cap.radius + 2 * cap.radius * δ
    grown >= zero(T) || throw(DomainError(factor, "cap growth factor produces a negative radius"))
    radius = grown > cap.radius ? min(T(π), nextfloat(grown)) : grown
    return SphericalCap(cap.point, radius)
end

"""
    Extents.union(x::SphericalCap, y::SphericalCap; strict=false) -> SphericalCap

Return a cap containing both `x` and `y`.

If either input already contains the other it is returned unchanged.
Otherwise, the centre lies between the input caps' centres, and the radius is
rounded outward to preserve coverage.

`strict` is accepted for consistency with the Extents API and has no effect
because both inputs use the same domain.
"""
function Extents.union(x::SphericalCap, y::SphericalCap; strict=false)
    Extents.contains(x, y; strict) && return x
    Extents.contains(y, x; strict) && return y

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
