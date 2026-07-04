# # Cap and arc extents

#=
3D `Extents.Extent`s of spherical objects — caps and great-circle arcs —
and conservative intersection between a [`SphericalCap`](@ref) and a 3D
box in unit-spherical ℝ³ space, so caps can drive `SpatialTreeInterface`
queries over trees whose leaves are boxes around spherical geometry.

The cap is the closed region `{p : ‖p‖ = 1, p ⋅ center ≥ radiuslike}`,
defined by the *stored* floats.  `Extents.intersects(cap, ext)` is a
filter: `false` proves the box contains no cap point, `true` means
"possibly intersecting".  Every separation test is monotone under box
enlargement, so pruning an internal tree node (whose extent is the union
of its children) can never discard a leaf the leaf-level test would keep.
=#

# The X/Y/Z bound pairs of an extent, erroring clearly on non-3D extents.
function _xyz_bounds(ext::Extents.Extent)
    K = keys(ext)
    (:X in K && :Y in K && :Z in K) || throw(ArgumentError(
        "cap–extent intersection works in unit-spherical ℝ³ and needs an extent with X, Y and Z dimensions, got $K"))
    return ext.X, ext.Y, ext.Z
end

"""
    Extents.extent(cap::SphericalCap)

The 3D axis-aligned bounding box of the cap, as an
`Extents.Extent{(:X, :Y, :Z)}`.  Conservative: bounds are widened outward
by ~`2√eps` (and clamped to `[-1, 1]`), so every cap point lies inside at
the cost of a loose last few digits.
"""
function Extents.extent(cap::SphericalCap)
    F = float(typeof(cap.radiuslike))
    c = cap.point
    n = sqrt(F(c[1])^2 + F(c[2])^2 + F(c[3])^2)
    k = F(cap.radiuslike) / n
    return Extents.Extent(
        X = _cap_axis_bounds(F(c[1]) / n, k),
        Y = _cap_axis_bounds(F(c[2]) / n, k),
        Z = _cap_axis_bounds(F(c[3]) / n, k),
    )
end

#=
Bounds of `pᵢ` over the cap: `cᵢ k ∓ √((1 − cᵢ²)(1 − k²))`, i.e.
`cos(θᵢ ± r)` for `cos θᵢ = cᵢ`, clamped to `±1` when the axis pole lies
inside the cap.  The outward pad dominates every rounding error here,
including the `√` of a cancellation-prone difference (absolute error up to
`~√eps`) and the not-exactly-unit stored center.
=#
function _cap_axis_bounds(ci::F, k::F) where F
    s = sqrt(max(one(F) - ci * ci, zero(F)) * max(one(F) - k * k, zero(F)))
    pad = 2 * sqrt(eps(F))
    hi = ci >= k ? one(F) : min(muladd(ci, k, s) + pad, one(F))
    lo = -ci >= k ? -one(F) : max(muladd(ci, k, -s) - pad, -one(F))
    return (lo, hi)
end

"""
    arc_extent(a, b)

The 3D axis-aligned bounding box of the minor great-circle arc from `a`
to `b` (unit vectors), as an `Extents.Extent{(:X, :Y, :Z)}`.
Conservative: the endpoints' box is padded by the arc's sagitta
`1 − cos(θ/2)` — the arc's maximum deviation from its chord — plus
rounding slack.  Exactly antipodal endpoints (where the minor arc is
ambiguous) get a pad of 1, covering every candidate arc.
"""
function arc_extent(a::UnitSphericalPoint, b::UnitSphericalPoint)
    F = float(promote_type(eltype(a), eltype(b)))
    d = clamp(F(a ⋅ b), -one(F), one(F))                        # cos θ
    pad = one(F) - sqrt((one(F) + d) / 2) + 4 * sqrt(eps(F))    # sagitta + slack
    return Extents.Extent(
        X = (min(a[1], b[1]) - pad, max(a[1], b[1]) + pad),
        Y = (min(a[2], b[2]) - pad, max(a[2], b[2]) + pad),
        Z = (min(a[3], b[3]) - pad, max(a[3], b[3]) + pad),
    )
end
arc_extent(a, b) = arc_extent(UnitSphericalPoint(a), UnitSphericalPoint(b))

"""
    Extents.intersects(cap::SphericalCap, ext::Extents.Extent)
    Extents.intersects(ext::Extents.Extent, cap::SphericalCap)

Whether the cap — the closed region `{p : ‖p‖ = 1, p ⋅ center ≥
radiuslike}` — may intersect the 3D box `ext`.  Conservative: `false` is a
proof of disjointness, established by exact-arithmetic separations against
the cap's half-space and the unit-sphere shell plus the cap's bounding
box; `true` means possibly intersecting.  Suitable as a spatial tree query
predicate (see `SpatialTreeInterface.sanitize_predicate`).
"""
Extents.intersects(cap::SphericalCap, ext::Extents.Extent) =
    _cap_intersects_box(cap, Extents.extent(cap), ext)
Extents.intersects(ext::Extents.Extent, cap::SphericalCap) =
    Extents.intersects(cap, ext)

function _cap_intersects_box(cap::SphericalCap, capbox::Extents.Extent, ext::Extents.Extent)
    bx, by, bz = _xyz_bounds(ext)
    xlo, xhi = Float64(bx[1]), Float64(bx[2])
    ylo, yhi = Float64(by[1]), Float64(by[2])
    zlo, zhi = Float64(bz[1]), Float64(bz[2])
    cx, cy, cz = Float64(cap.point[1]), Float64(cap.point[2]), Float64(cap.point[3])
    k = Float64(cap.radiuslike)
    # The cap lies in the half-space `p ⋅ c ≥ k`: disjoint when even the box
    # corner most aligned with `c` falls short.
    sx = cx >= 0 ? xhi : xlo
    sy = cy >= 0 ? yhi : ylo
    sz = cz >= 0 ? zhi : zlo
    _sign_dot3mk(sx, sy, sz, cx, cy, cz, k) < 0 && return false
    # The cap lies on the unit sphere: disjoint when the box misses the
    # shell (lies strictly inside or strictly outside it).
    nx = xlo > 0 ? xlo : (xhi < 0 ? xhi : 0.0)
    ny = ylo > 0 ? ylo : (yhi < 0 ? yhi : 0.0)
    nz = zlo > 0 ? zlo : (zhi < 0 ? zhi : 0.0)
    _sign_sqnorm3m1(nx, ny, nz) > 0 && return false
    fx = abs(xlo) >= abs(xhi) ? xlo : xhi
    fy = abs(ylo) >= abs(yhi) ? ylo : yhi
    fz = abs(zlo) >= abs(zhi) ? zlo : zhi
    _sign_sqnorm3m1(fx, fy, fz) < 0 && return false
    # The cap lies in its own bounding box.
    return Extents.intersects(capbox, ext)
end

# A tree-query predicate with the cap's bounding box hoisted out of the
# traversal, so no `sqrt` is paid per visited node.
struct _CapIntersects{C <: SphericalCap, E <: Extents.Extent}
    cap::C
    capbox::E
end
(p::_CapIntersects)(ext::Extents.Extent) = _cap_intersects_box(p.cap, p.capbox, ext)

SpatialTreeInterface.sanitize_predicate(cap::SphericalCap) =
    _CapIntersects(cap, Extents.extent(cap))

# ## Exact sign kernels

#=
ExactPredicates' `@genpredicate` filter assumes multihomogeneous
polynomials, which `s ⋅ c − k` is not (the `k` term has degree 0 in the
coordinates), so these use the same two-stage structure by hand: a float
evaluation against a forward error bound (`γ₄ ≈ 2 eps` relative to the sum
of term magnitudes, doubled for slack), then exact `Rational{BigInt}`
arithmetic when inconclusive.  Floats convert to rationals exactly, so the
fallback — and therefore the sign — is exact.
=#

# Exact rational lift; `Float64(x)` is exact for the float types we take.
_rat(x::Real) = Rational{BigInt}(Float64(x))

# sign(s ⋅ c − k), exact.
function _sign_dot3mk(sx, sy, sz, cx, cy, cz, k)
    tx, ty, tz = sx * cx, sy * cy, sz * cz
    val = (tx + ty + tz) - k
    err = 4 * eps(Float64) * ((abs(tx) + abs(ty) + abs(tz)) + abs(k))
    abs(val) > err && return val > 0 ? 1 : -1
    return Int(sign((_rat(sx) * _rat(cx) + _rat(sy) * _rat(cy) + _rat(sz) * _rat(cz)) - _rat(k)))
end

# sign(‖v‖² − 1), exact.
function _sign_sqnorm3m1(vx, vy, vz)
    tx, ty, tz = vx * vx, vy * vy, vz * vz
    val = (tx + ty + tz) - 1.0
    err = 4 * eps(Float64) * (tx + ty + tz + 1.0)
    abs(val) > err && return val > 0 ? 1 : -1
    return Int(sign(_rat(vx)^2 + _rat(vy)^2 + _rat(vz)^2 - 1))
end
