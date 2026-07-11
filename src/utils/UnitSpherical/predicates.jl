# # Spherical Predicates
#=
This file contains geometric predicates for spherical geometry on the unit sphere.
These predicates determine spatial relationships between points and arcs on the sphere.
=#

"""
    spherical_orient(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint) -> Int

Determine the orientation of point `c` with respect to the great circle arc from `a` to `b`.

Returns:
- `1` if `c` is to the left of the arc (counter-clockwise)
- `-1` if `c` is to the right of the arc (clockwise)
- `0` if `c` is on the great circle (collinear)

Uses [`robust_cross_product`](@ref) for numerical stability with nearly identical or antipodal points.

# Examples
```jldoctest
using GeometryOps.UnitSpherical: UnitSphericalPoint, spherical_orient
a = UnitSphericalPoint(1.0, 0.0, 0.0)
b = UnitSphericalPoint(0.0, 1.0, 0.0)
c = UnitSphericalPoint(0.0, 0.0, 1.0)
spherical_orient(a, b, c)
# output
1
```
"""
function spherical_orient(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    # The orientation is determined by sign((a × b) · c)
    # Use robust_cross_product for numerical stability
    n = robust_cross_product(a, b)
    dot_product = n ⋅ c

    # Use a tolerance for near-zero values
    tol = eps(Float64) * 16  # Same tolerance as S2 geometry
    if abs(dot_product) < tol
        return 0
    end
    return dot_product > 0 ? 1 : -1
    # return ExactPredicates.orient(a, b, UnitSphericalPoint((0., 0., 0.)), c)
end

# Convenience method for raw vectors
function spherical_orient(a::AbstractVector, b::AbstractVector, c::AbstractVector)
    return spherical_orient(
        UnitSphericalPoint(a),
        UnitSphericalPoint(b),
        UnitSphericalPoint(c)
    )
end

"""
    point_on_spherical_arc(p::UnitSphericalPoint, a::UnitSphericalPoint, b::UnitSphericalPoint) -> Bool

Check if point `p` lies on the great circle arc from `a` to `b`.

The arc is the shorter path along the great circle connecting `a` and `b`.
Returns `true` if `p` is on the arc (including endpoints), `false` otherwise.

# Examples
```jldoctest
using GeometryOps.UnitSpherical: UnitSphericalPoint, point_on_spherical_arc
a = UnitSphericalPoint(1.0, 0.0, 0.0)
b = UnitSphericalPoint(0.0, 1.0, 0.0)
mid = UnitSphericalPoint(1/√2, 1/√2, 0.0)
point_on_spherical_arc(mid, a, b)
# output
true
```
"""
function point_on_spherical_arc(p::UnitSphericalPoint, a::UnitSphericalPoint, b::UnitSphericalPoint)
    # First check: is p on the great circle through a and b?
    if spherical_orient(a, b, p) != 0
        return false
    end

    # Second check: is p between a and b on the arc?
    # For the shorter arc, p is between a and b if:
    # (a · p) ≥ (a · b) and (b · p) ≥ (a · b)
    # This works because dot product on unit sphere = cos(angle)
    # If p is between a and b, the angles a-p and b-p are both ≤ angle a-b

    ab = a ⋅ b  # cos(angle between a and b)
    ap = a ⋅ p  # cos(angle between a and p)
    bp = b ⋅ p  # cos(angle between b and p)

    tol = eps(Float64) * 16

    # p is on arc if it's "closer" to both endpoints than they are to each other
    # (in terms of angle, so larger dot product)
    return (ap ≥ ab - tol) && (bp ≥ ab - tol)
end

"""
    spherical_ring_contains(pts, n, q; orient, on_arc, proper_crossing) -> Union{Bool, Nothing}

Whether `q` lies in the closed region on the left of the ring `pts[1:n]`
(S2 loop convention: counterclockwise winding, interior on the left, so a
clockwise ring contains the complement).  The closing edge `pts[n] → pts[1]`
is implied; boundary points count as contained.  Returns `nothing` when
every anchor edge is degenerate with respect to `q` — callers must treat
that conservatively.

Containment is decided by crossing parity, the way `S2Loop::Contains` /
`InitBound` decide pole containment: which side of an anchor edge `q` falls
on, flipped once per transversal crossing of the arc from the anchor's
midpoint to `q` with the other edges; degenerate anchors are skipped and
the next edge tried.

The geometric predicates are injectable, for callers with stricter
requirements.  They receive the input points untouched (which may be
non-unit for scale-invariant predicates — the defaults assume unit input);
only the constructed reference midpoint is normalized.

- `orient(a, b, c)`: sign-valued orientation of `c` against the oriented
  great circle through `a, b`; default [`spherical_orient`](@ref).
- `on_arc(q, a, b)::Bool`: boundary membership; default
  [`point_on_spherical_arc`](@ref).  Pass `Returns(false)` when boundary
  points are already classified.
- `proper_crossing(q, m, a, b)::Int`: `1` if the minor arcs `(q, m)` and
  `(a, b)` cross transversally in both interiors, `0` if not, `-1` for too
  close to call; consulted once `orient` places both endpoint pairs
  strictly transversally.  The default uses `robust_cross_product` with a
  small tolerance band.
"""
function spherical_ring_contains(pts, n, q;
        orient = spherical_orient,
        on_arc = point_on_spherical_arc,
        proper_crossing = _hemisphere_proper_crossing)
    for j in 1:n
        on_arc(q, pts[j], pts[mod1(j + 1, n)]) && return true
    end
    nq = norm(q)
    for j in 1:n
        a, b = pts[j], pts[mod1(j + 1, n)]
        a == b && continue
        side = orient(a, b, q)
        side == 0 && continue
        mid = a + b
        # near-antipodal edge: the midpoint direction is unstable
        norm(mid) < 1e-9 * (norm(a) + norm(b)) && continue
        m = UnitSphericalPoint(normalize(mid))
        # test arc q → m would span a half turn
        dot(q, m) < (-1 + 1e-9) * nq && continue
        crossings = 0
        ok = true
        for k in 1:n
            k == j && continue
            c = _arc_crossing_parity(q, m, pts[k], pts[mod1(k + 1, n)]; orient, proper_crossing)
            if c == -1
                ok = false
                break
            end
            crossings += c
        end
        ok || continue
        # walking from `m` toward `q` departs onto `q`'s side of the anchor
        # edge (the arc meets that great circle again only at `-m`); positive
        # side is the interior, and each crossing flips it
        return isodd(crossings) ? side < 0 : side > 0
    end
    return nothing
end

# Crossing parity of the test arc q → m against ring edge a → b: 1 for a
# transversal crossing, 0 for none, -1 for too close to degenerate to call
# (with an exact `orient`, only exact incidences return -1).
function _arc_crossing_parity(q, m, a, b; orient, proper_crossing)
    # a vertex at `-q` lies on every great circle through `q`; its edges can
    # reach the test arc only at `q` itself, excluded by the on-boundary check
    (a == -q || b == -q) && return 0
    sa = orient(q, m, a)
    sb = orient(q, m, b)
    (sa == 0 || sb == 0) && return -1
    (sa > 0) == (sb > 0) && return 0
    # `q` on this edge's great circle but off the edge (checked upfront): the
    # circles meet only at `±q`, out of the test arc's reach — no crossing.
    # Anchor-independent (lonlat meridian edges hold `±eₓ`/`±e_y` exactly),
    # so resolve instead of returning -1.
    sq = orient(a, b, q)
    sq == 0 && return 0
    sm = orient(a, b, m)
    sm == 0 && return -1
    (sq > 0) == (sm > 0) && return 0
    return proper_crossing(q, m, a, b)
end

# Default transversality decision: the circles' intersection direction `x`
# must point into both arcs' hemispheres (each arc holds exactly one of `±x`
# once the endpoint sides are strict). Tolerance-banded; assumes unit input.
function _hemisphere_proper_crossing(q, m, a, b)
    x = cross(normalize(robust_cross_product(q, m)),
              normalize(robust_cross_product(a, b)))
    d1 = dot(x, q + m)
    d2 = dot(x, a + b)
    tol = 16 * eps(Float64) * norm(x)
    (abs(d1) <= tol || abs(d2) <= tol) && return -1
    return (d1 > 0) == (d2 > 0) ? 1 : 0
end
