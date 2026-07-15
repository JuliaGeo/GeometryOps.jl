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

"""
    spherical_exterior_anchor(pts, n) -> Union{UnitSphericalPoint{Float64}, Nothing}

A reference point exterior BY DEFINITION of the enclosed-region semantics
of the ring `pts[1:n]`: the antipode of the ring's normalized vertex mass
(the sum of the unit vertex directions). For any ring whose enclosed region
is meaningfully smaller than a hemisphere, the vertex mass points into the
cap the vertices bound, so its antipode lies in the larger — exterior —
region.

Returns `nothing` when the mass norm is tiny (below `1e-6` per vertex):
near-hemisphere or vertex-symmetric rings, whose vertices spread over a
near-great circle. There the enclosed/complement distinction is itself
near-degenerate (the turning-angle winding tolerance already treats exact
hemispheres permissively — see `_ring_is_ccw`), so callers fall back to the
winding-consistent wedge bootstrap of [`spherical_ring_contains`](@ref).
"""
function spherical_exterior_anchor(pts, n)
    n == 0 && return nothing
    mass = normalize(SVector{3, Float64}(pts[1]))
    for i in 2:n
        mass += normalize(SVector{3, Float64}(pts[i]))
    end
    norm(mass) <= 1e-6 * n && return nothing
    return UnitSphericalPoint(-normalize(mass))
end

"""
    spherical_ring_encloses(pts, n, q;
        anchor, orient, on_arc, proper_crossing) -> Union{Bool, Nothing}

Whether `q` lies in the region ENCLOSED by the ring `pts[1:n]` (the closing
edge `pts[n] → pts[1]` is implied; boundary points count as enclosed):
even-odd crossing parity of the arc from `q` to a reference point that is
exterior *by definition* of the enclosed-region semantics — `anchor`, by
default the antipode of the normalized vertex mass
([`spherical_exterior_anchor`](@ref)).

Winding-independent, like [`spherical_ring_contains`](@ref) composed with a
winding test — but where that composition bootstraps the interior from a
local wedge at one edge and a global turning-angle sum, both of which a
ring that self-intersects *on the sphere* defeats (a figure-eight's lobes
cancel the turning angle, and the wedge answer is anchored to whichever
lobe hosts the edge — S2's forced-through behavior, globally inverted on
real data), parity from a definitionally exterior point degrades to
even-odd semantics: both lobes enclosed, the far side out.

Returns `nothing` — callers fall back conservatively — when:

- `anchor === nothing` (degenerate vertex mass, see
  [`spherical_exterior_anchor`](@ref));
- `q` is (nearly) antipodal to the anchor (the test arc is ill-defined:
  `q` sits at the center of the vertex mass);
- the anchor lies exactly ON a ring edge (the test arc ends on the ring);
  or
- `proper_crossing` reports a crossing as too close to call (`-1`; never
  with exact injected predicates).

The `orient`/`on_arc`/`proper_crossing` predicates are injectable exactly
as in [`spherical_ring_contains`](@ref); `on_test_arc(v, a, b)` decides
whether a point already known to lie on the great circle of `(a, b)` lies
on the closed minor arc (the vertex-grazing resolution below — exact
callers inject their span test).
"""
function spherical_ring_encloses(pts, n, q;
        anchor = spherical_exterior_anchor(pts, n),
        orient = spherical_orient,
        on_arc = point_on_spherical_arc,
        on_test_arc = point_on_spherical_arc,
        proper_crossing = _hemisphere_proper_crossing)
    for j in 1:n
        on_arc(q, pts[j], pts[mod1(j + 1, n)]) && return true
    end
    anchor === nothing && return nothing
    # test arc q → anchor would span (nearly) a half turn
    dot(q, anchor) < (-1 + 1e-9) * norm(q) && return nothing
    crossings = 0
    for k in 1:n
        c = _anchor_crossing_parity(q, anchor, pts[k], pts[mod1(k + 1, n)];
            orient, on_test_arc, proper_crossing)
        c == -1 && return nothing
        crossings += c
    end
    return isodd(crossings)
end

#=
Crossing parity of the closed test arc q → z (z the definitional exterior
anchor) against ring edge a → b: `_arc_crossing_parity` with the two
exactly-degenerate configurations that helper refuses (-1) resolved the way
the indexed locator's `count_arc_segment!` resolves them — symbolically, S2
`VertexCrossing` style — so a symmetric ring (whose vertex mass can point
exactly at a crossing point, putting the anchor on an edge's great circle)
cannot force every query back onto the wedge bootstrap:

- an edge endpoint exactly on the test arc's great circle (`sa == 0` /
  `sb == 0`): two distinct great circles meet only at one antipodal pair,
  so the edge can touch the test arc only at that endpoint — count iff the
  endpoint lies ON the closed test arc and the other endpoint is strictly
  on the positive side, so a crossing pair of incident edges counts once
  and a same-side pair counts zero or twice (parity-equal);
- the anchor exactly on the edge's great circle (`sm == 0`): the circles
  meet only at ±z, and the minor test arc reaches z but never −z — no
  crossing, unless the edge itself contains z (the anchor ON the ring:
  refuse with -1, the caller falls back).

`sq == 0` (q on the edge's circle but not on the edge — boundary is
excluded upfront) stays 0, and edges with a vertex at −q stay 0, exactly
as in `_arc_crossing_parity`.
=#
function _anchor_crossing_parity(q, z, a, b; orient, on_test_arc, proper_crossing)
    (a == -q || b == -q) && return 0
    a == b && return 0
    sa = orient(q, z, a)
    sb = orient(q, z, b)
    if sa == 0 || sb == 0
        if sa == 0 && sb == 0
            # edge collinear with the test circle: its neighbors decide the
            # parity — unless it holds the anchor itself
            return on_test_arc(z, a, b) ? -1 : 0
        end
        von, s_off = sa == 0 ? (a, sb) : (b, sa)
        return (s_off > 0 && on_test_arc(von, q, z)) ? 1 : 0
    end
    (sa > 0) == (sb > 0) && return 0
    sq = orient(a, b, q)
    sq == 0 && return 0
    sm = orient(a, b, z)
    if sm == 0
        return on_test_arc(z, a, b) ? -1 : 0
    end
    (sq > 0) == (sm > 0) && return 0
    return proper_crossing(q, z, a, b)
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
