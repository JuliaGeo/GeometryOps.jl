# # Spherical Caps

#=
```@meta
CollapsedDocStrings = true
```

```@docs; canonical=false
SphericalCap
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
using GeometryOps
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

#=
## Cap predicates

The cap is the closed region `{p : ‖p‖ = 1, p ⋅ c ≥ k}` for the stored
center `c` and `k = radiuslike` — a geodesic ball of angular radius
`arccos(k/‖c‖)` around `c/‖c‖`.  The predicates decide relations between
these regions exactly: a float evaluation with a `1e-7` screen (two
orders above the worst rounding error of these expressions, which the
`√` of cancellation-prone products caps at ~4e-8), then exact
`Rational{BigInt}` case analysis on the same closed-form conditions.
Everything is closed-set — tangency counts as intersecting/contained,
matching S2's `S2Cap` conventions.

With `Sx = cx ⋅ cx`, `D = cx ⋅ cy`, and radii `rᵢ = arccos(kᵢ/‖cᵢ‖)`:
- empty(X) ⟺ kx > ‖cx‖;  full(X) ⟺ kx ≤ −‖cx‖
- intersects(X, Y) ⟺ d ≤ rx + ry: always true when `rx + ry ≥ π`
  (⟺ kx‖cy‖ + ky‖cx‖ ≤ 0), else
  `(D − kx ky) + √((Sx − kx²)(Sy − ky²)) ≥ 0`.
- contains(B, S) ⟺ d + rs ≤ rb: needs `D + ks‖cb‖ ≥ 0` (else
  `d + rs > π`) and `D ks − kb Ss ≥ √((Sb Ss − D²)(Ss − ks²))`.
=#

const _CAP_SCREEN_TOL = 1e-7

function _intersects(x::SphericalCap, y::SphericalCap)
    cx, cy = Float64.(Tuple(x.point)), Float64.(Tuple(y.point))
    kx, ky = Float64(x.radiuslike), Float64(y.radiuslike)
    Sx, Sy = sum(abs2, cx), sum(abs2, cy)
    nx, ny = sqrt(Sx), sqrt(Sy)
    τ = _CAP_SCREEN_TOL
    if kx < nx - τ && ky < ny - τ               # both certainly nonempty
        g = kx * ny + ky * nx
        g < -τ && return true                   # rx + ry certainly > π
        if g > τ
            D = sum(cx .* cy)
            e = (D - kx * ky) + sqrt(max((Sx - kx^2) * (Sy - ky^2), 0.0))
            abs(e) > τ && return e > 0
        end
    elseif kx > nx + τ || ky > ny + τ
        return false                            # certainly empty
    end
    return _intersects_exact(x, y)
end

_disjoint(x::SphericalCap, y::SphericalCap) = !_intersects(x, y)

function _intersects_exact(x::SphericalCap, y::SphericalCap)
    cx, cy = _rat.(Tuple(x.point)), _rat.(Tuple(y.point))
    kx, ky = _rat(x.radiuslike), _rat(y.radiuslike)
    Sx, Sy = sum(abs2, cx), sum(abs2, cy)
    (kx > 0 && kx^2 > Sx) && return false       # X empty
    (ky > 0 && ky^2 > Sy) && return false       # Y empty
    _sign_sqrtpair(kx, Sy, ky, Sx) <= 0 && return true   # rx + ry ≥ π
    G = sum(cx .* cy) - kx * ky
    G >= 0 && return true
    return (Sx - kx^2) * (Sy - ky^2) >= G^2
end

function _contains(big::SphericalCap, small::SphericalCap)
    cb, cs = Float64.(Tuple(big.point)), Float64.(Tuple(small.point))
    kb, ks = Float64(big.radiuslike), Float64(small.radiuslike)
    Sb, Ss = sum(abs2, cb), sum(abs2, cs)
    nb, ns = sqrt(Sb), sqrt(Ss)
    τ = _CAP_SCREEN_TOL
    if ks > ns + τ
        return true                             # small certainly empty
    elseif ks < ns - τ && abs(kb) < nb - τ      # small nonempty, big neither empty nor full
        D = sum(cb .* cs)
        g = D + ks * nb
        g < -τ && return false                  # d + rs certainly > π
        if g > τ
            H = D * ks - kb * Ss
            e = H - sqrt(max((Sb * Ss - D^2) * (Ss - ks^2), 0.0))
            abs(e) > τ && return e > 0
        end
    end
    return _contains_exact(big, small)
end

function _contains_exact(big::SphericalCap, small::SphericalCap)
    cb, cs = _rat.(Tuple(big.point)), _rat.(Tuple(small.point))
    kb, ks = _rat(big.radiuslike), _rat(small.radiuslike)
    Sb, Ss = sum(abs2, cb), sum(abs2, cs)
    (ks > 0 && ks^2 > Ss) && return true        # small empty
    (kb > 0 && kb^2 > Sb) && return false       # big empty
    (kb <= 0 && kb^2 >= Sb) && return true      # big is the whole sphere
    D = sum(cb .* cs)
    _sign_a_plus_bsqrt(D, ks, Sb) < 0 && return false    # d + rs > π
    H = D * ks - kb * Ss
    H < 0 && return false
    return H^2 >= (Sb * Ss - D^2) * (Ss - ks^2)
end

# Half-space membership of the point taken as-is (its norm is not
# checked), exact w.r.t. the stored floats.
_contains(cap::SphericalCap, point::UnitSphericalPoint) =
    _sign_dot3mk(Float64(point[1]), Float64(point[2]), Float64(point[3]),
        Float64(cap.point[1]), Float64(cap.point[2]), Float64(cap.point[3]),
        Float64(cap.radiuslike)) >= 0

# sign(u√A + v√B) for A, B ≥ 0 — exact on rationals.
function _sign_sqrtpair(u, A, v, B)
    su = (u == 0 || A == 0) ? 0 : Int(sign(u))
    sv = (v == 0 || B == 0) ? 0 : Int(sign(v))
    su == 0 && return sv
    sv == 0 && return su
    su == sv && return su
    c = u^2 * A - v^2 * B
    return c == 0 ? 0 : (c > 0 ? su : sv)
end

# sign(a + b√S) for S ≥ 0 — exact on rationals.
function _sign_a_plus_bsqrt(a, b, S)
    sb = (b == 0 || S == 0) ? 0 : Int(sign(b))
    sa = Int(sign(a))
    sb == 0 && return sa
    sa == 0 && return sb
    sa == sb && return sa
    c = a^2 - b^2 * S
    return c == 0 ? 0 : (c > 0 ? sa : sb)
end

#Comment by asinghvi: this could be transformed to GO.union
function _merge(x::SphericalCap, y::SphericalCap)

    d = spherical_distance(x.point, y.point)
    newradius = (x.radius + y.radius + d) / 2
    if newradius < x.radius
        #x contains y
        x
    elseif newradius < y.radius
        #y contains x
        y
    else
        excenter = 0.5 * (1 - (x.radius - y.radius) / d)
        newcenter = slerp(x.point, y.point, excenter)
        SphericalCap(newcenter, newradius)
    end
end

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
