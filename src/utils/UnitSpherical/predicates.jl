# # Spherical Predicates
#=
This file contains geometric predicates for spherical geometry on the unit sphere.
These predicates determine spatial relationships between points and arcs on the sphere.
=#

"""
    spherical_orient(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint) -> Int

Classify `c` against the directed great-circle arc `a → b`.

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

# Extended help

## Why this does not simply call `robust_cross_product`

Use unnormalized `cross(a - b, a + b)` and a squared degeneracy test to avoid normalization
and square roots. Boundary rounding may change zero to a sign, but cannot reverse a sign.

Nearly equal or antipodal endpoints use [`robust_cross_product`](@ref).
"""
function spherical_orient(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    # The orientation is determined by sign((a × b) · c).
    #
    # Use the stable unnormalized `cross(a - b, a + b) = 2(a × b)`, expanded componentwise.
    d1 = a[1] - b[1]; d2 = a[2] - b[2]; d3 = a[3] - b[3]
    s1 = a[1] + b[1]; s2 = a[2] + b[2]; s3 = a[3] + b[3]
    n1 = d2 * s3 - d3 * s2
    n2 = d3 * s1 - d1 * s3
    n3 = d1 * s2 - d2 * s1
    nsqr = n1 * n1 + n2 * n2 + n3 * n3
    # Same stability criterion `robust_cross_product` applies internally.
    kmin = min_stable_norm(promote_type(eltype(a), eltype(b)))
    if nsqr >= kmin * kmin
        dot_product = n1 * c[1] + n2 * c[2] + n3 * c[3]
        tol = eps(Float64) * 16  # Same tolerance as S2 geometry
        # `abs(dot_product) / sqrt(nsqr) < tol`, without the sqrt
        dot_product * dot_product < (tol * tol) * nsqr && return 0
        return dot_product > 0 ? 1 : -1
    end

    # Recover the normal of nearly equal or antipodal endpoints with exact arithmetic and
    # symbolic perturbation.
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

# ## exact_spherical_orient

#= Expand the determinant because `ExactPredicates.det` is internal.
Pass `a` through two formal parameters to keep the accumulator multihomogeneous:
each input gate can belong to only one group.

The bound holds for all four input tuples, including the equal-parameter case.
The exact polynomial remains `det[a; b-a; c-a]`. =#
ExactPredicates.Codegen.@genpredicate function _exact_spherical_orient(
        a1 :: 3, a2 :: 3, b :: 3, c :: 3)
    u = b - a2
    v = c - a2
    ExactPredicates.Codegen.group!(a1...)
    ExactPredicates.Codegen.group!(u...)
    ExactPredicates.Codegen.group!(v...)
    a1[1] * (u[2]*v[3] - u[3]*v[2]) -
    a1[2] * (u[1]*v[3] - u[3]*v[1]) +
    a1[3] * (u[1]*v[2] - u[2]*v[1])
end

#= ExactPredicates requires `Float64`. Widening `Float16` or `Float32` is lossless,
so the result remains the exact sign of the input determinant. =#
@inline _ep_widen(p) = (Float64(p[1]), Float64(p[2]), Float64(p[3]))
@inline _ep_widen(p::UnitSphericalPoint{Float64}) = p
@inline _ep_widen(p::NTuple{3, Float64}) = p

"""
    exact_spherical_orient(a, b, c) -> Int

Return the exact sign of `(a × b) · c`: 1 for left, -1 for right, and 0 for coplanarity with
the origin. Inputs need not be unit vectors.

# Extended help

## Why this exists alongside `spherical_orient`

The `eps*16` band in [`spherical_orient`](@ref) can classify cell-scale crossings as
collinear. Exact signs preserve their topology.

## Why not `ExactPredicates.orient(a, b, c, (0, 0, 0))`

Grouping unit vectors separately gives a filter bound near `5e-15`, even for nearby points.
Instead use `det[a; b; c] == det[a; b-a; c-a]`, whose bound scales with separation squared.

Exact orientation costs about twice a plain triple product and 1.2 times
[`spherical_orient`](@ref).

Unresolved signs use exact arithmetic. Non-finite coordinates are outside the contract and
throw.

## Narrower float widths

Widen `Float16` and `Float32` losslessly to `Float64`; preserve `Float64` inputs unchanged.
"""
@inline function exact_spherical_orient(a, b, c)
    a64 = _ep_widen(a)
    return _exact_spherical_orient(a64, a64, _ep_widen(b), _ep_widen(c))
end

"""
    point_on_spherical_arc(p::UnitSphericalPoint, a::UnitSphericalPoint, b::UnitSphericalPoint) -> Bool

Return whether `p` lies on the shorter great-circle arc from `a` to `b`, including endpoints.

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

    # Test the cosine span conditions `(a · p) ≥ (a · b)` and `(b · p) ≥ (a · b)`. Unit-vector
    # dot products are angle cosines.

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

Test the closed region left of `pts[1:n]`, with implicit closure. Boundary points count as
contained. Clockwise rings contain the complement under the S2 convention.

Use anchor-edge side and crossing parity, as in `S2Loop::Contains` / `InitBound`. Skip
degenerate anchors; return `nothing` if none works, requiring conservative handling.

Injected predicates receive unchanged input points; only the reference midpoint is normalized.
Defaults require unit input; scale-invariant replacements may accept non-unit directions.

- `orient(a, b, c)`: sign-valued orientation of `c` against the oriented
  great circle through `a, b`; default [`spherical_orient`](@ref).
- `on_arc(q, a, b)::Bool`: boundary membership; default
  [`point_on_spherical_arc`](@ref).  Pass `Returns(false)` when boundary
  points are already classified.
- `proper_crossing(q, m, a, b)::Int`: 1 for a transversal interior crossing, 0 for none, or -1
  if undecidable. Called after strict orientation straddling; defaults to tolerance-banded
  `robust_cross_product`.
"""
function spherical_ring_contains(pts, n, q;
        orient = spherical_orient,
        on_arc = point_on_spherical_arc,
        proper_crossing = _hemisphere_proper_crossing)
    return _ring_contains(pts, n, q, orient, on_arc, proper_crossing)
end

#= Bind predicate types explicitly so forwarded callbacks specialize without
dynamic dispatch or per-edge allocation. =#
function _ring_contains(pts, n, q, orient::O, on_arc::OA, proper_crossing::PC) where {O, OA, PC}
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

Return the antipode of normalized vertex mass, the exterior reference for enclosed-region
semantics. For rings well below a hemisphere, it lies outside the vertex cap.

Return `nothing` when mass norm is below `1e-6` per vertex. Near-hemisphere or symmetric rings
then require the winding-based fallback in [`spherical_ring_contains`](@ref).
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

Test winding-independent even-odd containment in `pts[1:n]`, including the boundary and
implicit closing edge. Count crossings from `q` to the exterior `anchor`.

The default anchor is [`spherical_exterior_anchor`](@ref). Parity gives even-odd semantics for
self-intersections, including both lobes of a figure-eight.

Return `nothing` in these cases; callers must fall back conservatively:

- `anchor === nothing` (degenerate vertex mass, see
  [`spherical_exterior_anchor`](@ref));
- `q` is (nearly) antipodal to the anchor (the test arc is ill-defined:
  `q` sits at the center of the vertex mass);
- the anchor lies exactly ON a ring edge (the test arc ends on the ring);
  or
- `proper_crossing` reports a crossing as too close to call (`-1`; never
  with exact injected predicates).

Inject predicates as in [`spherical_ring_contains`](@ref). `on_test_arc(v, a, b)` tests closed
minor-arc span for a point already known to be on the great circle.
"""
function spherical_ring_encloses(pts, n, q;
        anchor = spherical_exterior_anchor(pts, n),
        orient = spherical_orient,
        on_arc = point_on_spherical_arc,
        on_test_arc = point_on_spherical_arc,
        proper_crossing = _hemisphere_proper_crossing)
    _on_ring_boundary(pts, n, q, on_arc) && return true
    anchor === nothing && return nothing
    return _ring_encloses_parity(pts, n, q, anchor, orient, on_test_arc, proper_crossing)
end

# Boundary scan, with `on_arc` bound to a type parameter so it specializes.
function _on_ring_boundary(pts, n, q, on_arc::OA) where {OA}
    for j in 1:n
        on_arc(q, pts[j], pts[mod1(j + 1, n)]) && return true
    end
    return false
end

#= Pass the anchor positionally and bind predicate types to specialize the parity walk.
This avoids union-typed keyword values and dynamically dispatched callbacks per edge. =#
function _ring_encloses_parity(pts, n, q, z, orient::O, on_test_arc::OT,
        proper_crossing::PC) where {O, OT, PC}
    # test arc q → z would span (nearly) a half turn
    dot(q, z) < (-1 + 1e-9) * norm(q) && return nothing
    crossings = 0
    for k in 1:n
        c = _anchor_crossing_parity(q, z, pts[k], pts[mod1(k + 1, n)];
            orient, on_test_arc, proper_crossing)
        c == -1 && return nothing
        crossings += c
    end
    return isodd(crossings)
end

#=
Count crossings of the test arc `q → z` against edge `a → b`, resolving exact degeneracies
with S2-style symbolic vertex crossing.

For `sa == 0` or `sb == 0`, count the endpoint only if it lies on the closed test arc and its
neighbor lies strictly positive. Incident edges then preserve parity.

For `sm == 0`, the circles meet at `±z`; the test arc reaches only `z`. Return -1 if the edge
contains the anchor, otherwise 0.

Return 0 when `q` lies on the edge's circle but outside the edge, or an edge vertex
equals `−q`, as in `_arc_crossing_parity`.
=#
function _anchor_crossing_parity(q, z, a, b; orient::O, on_test_arc::OT,
        proper_crossing::PC) where {O, OT, PC}
    (a == -q || b == -q) && return 0
    a == b && return 0
    sa = orient(q, z, a)
    sb = orient(q, z, b)
    if sa == 0 || sb == 0
        if sa == 0 && sb == 0
            # Collinear edge: neighbors determine parity unless it contains the anchor.
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
    return _proper_crossing_with_orients(proper_crossing, orient, q, z, a, b, sa, sb, sq, sm)
end

# Return crossing parity of `q → m` against `a → b`: 1 for crossing, 0 for none, -1 if
# undecidable. Exact orientation makes only exact incidences undecidable.
function _arc_crossing_parity(q, m, a, b; orient::O, proper_crossing::PC) where {O, PC}
    # a vertex at `-q` lies on every great circle through `q`; its edges can
    # reach the test arc only at `q` itself, excluded by the on-boundary check
    (a == -q || b == -q) && return 0
    sa = orient(q, m, a)
    sb = orient(q, m, b)
    (sa == 0 || sb == 0) && return -1
    (sa > 0) == (sb > 0) && return 0
    # If `q` lies on the edge's circle but outside the edge, the circles meet only at `±q` and
    # cannot cross inside the test arc.
    sq = orient(a, b, q)
    sq == 0 && return 0
    sm = orient(a, b, m)
    sm == 0 && return -1
    (sq > 0) == (sm > 0) && return 0
    return _proper_crossing_with_orients(proper_crossing, orient, q, m, a, b, sa, sb, sq, sm)
end

# Allow opt-in callbacks to reuse four orientation signs. Other callbacks retain the
# four-point protocol and its undecidable result.
@inline function _proper_crossing_with_orients(proper_crossing::PC, orient::O,
        q, m, a, b, sa, sb, sq, sm) where {PC, O}
    return proper_crossing(q, m, a, b)
end

# Require the circle intersection direction to lie in both arc hemispheres. Strict endpoint
# signs select one candidate per arc. Uses a tolerance band and unit input.
function _hemisphere_proper_crossing(q, m, a, b)
    x = cross(normalize(robust_cross_product(q, m)),
              normalize(robust_cross_product(a, b)))
    d1 = dot(x, q + m)
    d2 = dot(x, a + b)
    tol = 16 * eps(Float64) * norm(x)
    (abs(d1) <= tol || abs(d2) <= tol) && return -1
    return (d1 > 0) == (d2 > 0) ? 1 : 0
end
