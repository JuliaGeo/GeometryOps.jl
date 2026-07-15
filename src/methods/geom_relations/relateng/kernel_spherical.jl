# # Spherical RelateKernel
#
#=
Spherical implementation of the RelateKernel contract declared in `kernel.jl`,
over `UnitSphericalPoint{Float64}`. Every predicate is a sign of
det(u, v, w) = (u×v)·w, so the exact path mirrors planar: a float filter then an
exact fallback (`ExactPredicates.orient` for the plain orient, `Rational{BigInt}`
on the xyz components for composites). No intersection coordinate is ever
constructed. See the design doc 2026-06-15.

All of `cross`, `⋅`, `normalize` (LinearAlgebra), `ExactPredicates`, and the
`UnitSpherical` names are already in scope here — this file is `include`d into
`GeometryOps`, which `using`s them at the top of the module. The
`_rebuild_point(::UnitSphericalPoint, …)` hook that keeps node points typed
lives next to the generic `_rebuild_point` in `kernel.jl`.
=#

# xyz tuple of a 3D point, for the ExactPredicates / Rational{BigInt} paths.
@inline _tup3(u) = (GI.x(u), GI.y(u), GI.z(u))

# ## rk_orient

# Orientation of `c` relative to the great-circle arc `(a, b)`: the sign of the
# scalar triple product (a×b)·c. Exact path: `ExactPredicates.orient` over the
# xyz tuples about the origin (the spike measured 3.5 ns); NOT
# `UnitSpherical.spherical_orient`, whose eps*16 tolerance is unfit for the
# exact contract. Float path: the plain triple product.
rk_orient(::Spherical, a, b, c; exact) = _rk_orient(booltype(exact), a, b, c)
@inline _rk_orient(::True, a, b, c) =
    ExactPredicates.orient(_tup3(a), _tup3(b), _tup3(c), (0.0, 0.0, 0.0))
@inline _rk_orient(::False, a, b, c) = cross(a, b) ⋅ c

# ## Exact-aware 3-vector arithmetic
#
# Composite predicates (arc membership, proper crossing, node coincidence)
# reduce to signs of polynomials in the xyz components. With `exact = True()` we
# evaluate over `Rational{BigInt}` (Float64 are dyadic rationals → exact); with
# `False()`, Float64. `_vec3(bt, p)` lifts a point to the chosen number type; the
# rest are plain tuple cross/dot, so one code path serves both — exactly how the
# planar kernel threads `exact`.
@inline _vec3(::True, u) = (Rational{BigInt}(GI.x(u)), Rational{BigInt}(GI.y(u)), Rational{BigInt}(GI.z(u)))
@inline _vec3(::False, u) = (Float64(GI.x(u)), Float64(GI.y(u)), Float64(GI.z(u)))
@inline _cross3(a, b) = (a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1])
@inline _dot3(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
@inline _iszero3(a) = iszero(a[1]) && iszero(a[2]) && iszero(a[3])
@inline _neg3(a) = (-a[1], -a[2], -a[3])
# `w` strictly interior to the minor arc (a, b) with normal n = a×b.
@inline _strictly_in_arc3(w, a, b, n) = _dot3(_cross3(a, w), n) > 0 && _dot3(_cross3(w, b), n) > 0
_usp_eq(p, q) = GI.x(p) == GI.x(q) && GI.y(p) == GI.y(q) && GI.z(p) == GI.z(q)

# ## rk_point_on_segment

# Whether `p` lies on the closed minor arc `[q0, q1]`. Two conditions: `p` is on
# the arc's great circle (coplanar with `q0`, `q1`, origin — an exact orient ==
# 0), and within the minor-arc span. The span test is scale-invariant: writing
# the coplanar `p` as `α q0 + β q1`, `p` is on the closed minor arc iff α, β ≥ 0,
# and `sign(β) = sign((q0×p)·n)`, `sign(α) = sign((p×q1)·n)` with `n = q0×q1` — a
# pure determinant sign, correct for unit and non-unit inputs alike.
function rk_point_on_segment(m::Spherical, p, q0, q1; exact)
    rk_orient(m, q0, q1, p; exact) == 0 || return false
    return _on_arc_span(booltype(exact), p, q0, q1)
end
@inline function _on_arc_span(bt, p, q0, q1)
    P = _vec3(bt, p); Q0 = _vec3(bt, q0); Q1 = _vec3(bt, q1)
    n = _cross3(Q0, Q1)
    if _iszero3(n)
        # parallel endpoints — a zero-length arc (real rings carry repeated
        # vertices; NE 110m North Korea has an `[A, A, B, A]` sliver ring) or
        # an ill-defined antipodal pair: the closed arc holds only its
        # endpoints, but with `n == 0` the span tests below are `0 >= 0` and
        # would accept every `p` on the great circle (which the orient gate
        # already reduced to every `p`, since orient against a zero normal is
        # identically 0). Membership is direction coincidence with an endpoint.
        return (_iszero3(_cross3(P, Q0)) && _dot3(P, Q0) > 0) ||
               (_iszero3(_cross3(P, Q1)) && _dot3(P, Q1) > 0)
    end
    return _dot3(_cross3(Q0, P), n) >= 0 && _dot3(_cross3(P, Q1), n) >= 0
end

# ## Ingest and interaction bounds

# Renormalize to unit length (Float32-sourced data — e.g. Natural Earth GeoJSON
# converted to Float64 — is ~1e-8 off unit and trips `robust_cross_product`).
@inline rk_normalize_usp(u) = UnitSphericalPoint(normalize(u))

# Canonical kernel point of a GeoInterface point: lon/lat (2D) → unit xyz, or an
# already-3D point treated as xyz; renormalized and signed-zero normalized so
# the same vertex always produces identical bits (NodeKey equality). The vertex
# ingest (Phase 3 `_to_kernel_point`) and the extent computation below share
# this, so a vertex and its extent agree exactly.
@inline function _spherical_kernel_point(p)
    u = GI.is3d(p) ?
        UnitSphericalPoint(Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p))) :
        UnitSphereFromGeographic()((Float64(GI.x(p)), Float64(GI.y(p))))
    return _node_point(rk_normalize_usp(u))
end

# Phase 3 ingest hooks (the planar methods live in kernel.jl). The spherical
# kernel point type is the unit-sphere xyz point; conversion is the canonical
# `_spherical_kernel_point` (lon/lat → unit xyz, or an already-xyz point
# renormalized), so an ingested vertex agrees bit-for-bit with its extent.
_kernel_point_type(::Spherical) = UnitSphericalPoint{Float64}
@inline _to_kernel_point(::Spherical, p) = _spherical_kernel_point(p)

@inline _widen(lo, hi) = (prevfloat(lo, 4), nextfloat(hi, 4))

@noinline _throw_antipodal_edge(a, b) = throw(ArgumentError(
    "spherical edge between antipodal vertices $(_tup3(a)) and $(_tup3(b)) has no " *
    "unique great-circle arc; densify it first with the `AntipodalEdgeSplit` " *
    "correction (it inserts the lon/lat midpoint)"))

# Exactly antipodal pair: vanishing cross product, opposed directions. A
# vanishing cross with `u ⋅ v > 0` is a zero-length/repeated vertex, fine.
_exactly_antipodal(u, v) = iszero(cross(u, v)) && (u ⋅ v) < 0.0

# Ingest validation, once per curve at `RelateGeometry` construction: an
# exactly-antipodal edge has no unique great-circle arc, so throw rather
# than pick one (`spherical_arc_extent` picks a stable plane, never throws).
function _validate_relate_edges(::Spherical, curve)
    n = GI.npoint(curve)
    n < 2 && return nothing
    prev = _spherical_kernel_point(GI.getpoint(curve, 1))
    for i in 2:n
        cur = _spherical_kernel_point(GI.getpoint(curve, i))
        _exactly_antipodal(prev, cur) && _throw_antipodal_edge(prev, cur)
        prev = cur
    end
    return nothing
end

# ## rk_classify_intersection
#
# The two great circles meet at ±d, d = (a0×a1)×(b0×b1). `SS_PROPER` iff one of
# ±d is strictly interior to both minor arcs (the candidate-direct formulation —
# planar straddle tests are NOT sufficient on the sphere, where arcs can straddle
# each other's great circle while meeting only at the antipodal point). Endpoint
# incidences are exact arc-membership; collinear = the arcs share a great circle
# (d == 0). No intersection coordinate is constructed.
function rk_classify_intersection(m::Spherical, a0, a1, b0, b1; exact)
    a0_on_b = rk_point_on_segment(m, a0, b0, b1; exact)
    a1_on_b = rk_point_on_segment(m, a1, b0, b1; exact)
    b0_on_a = rk_point_on_segment(m, b0, a0, a1; exact)
    b1_on_a = rk_point_on_segment(m, b1, a0, a1; exact)
    return _sph_classify(booltype(exact), a0, a1, b0, b1, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
end

function _sph_classify(bt, a0, a1, b0, b1, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    A0 = _vec3(bt, a0); A1 = _vec3(bt, a1); B0 = _vec3(bt, b0); B1 = _vec3(bt, b1)
    na = _cross3(A0, A1); nb = _cross3(B0, B1)
    d = _cross3(na, nb)
    n_inc = a0_on_b + a1_on_b + b0_on_a + b1_on_a
    if _iszero3(d)   # same great circle (or a degenerate, zero-length arc)
        n_inc == 0 && return SegSegClass(SS_DISJOINT, false, false, false, false)
        # a degenerate (zero-length) arc on the other is a touch, not an overlap
        zero_len = _iszero3(na) || _iszero3(nb)
        shared_only = n_inc == 2 && (a0_on_b || a1_on_b) && (b0_on_a || b1_on_a) &&
            (_usp_eq(a0, b0) || _usp_eq(a0, b1) || _usp_eq(a1, b0) || _usp_eq(a1, b1))
        kind = (shared_only || zero_len) ? SS_TOUCH : SS_COLLINEAR
        return SegSegClass(kind, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    end
    if a0_on_b || a1_on_b || b0_on_a || b1_on_a
        return SegSegClass(SS_TOUCH, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    end
    nd = _neg3(d)
    if (_strictly_in_arc3(d, A0, A1, na) && _strictly_in_arc3(d, B0, B1, nb)) ||
       (_strictly_in_arc3(nd, A0, A1, na) && _strictly_in_arc3(nd, B0, B1, nb))
        return SegSegClass(SS_PROPER, false, false, false, false)
    end
    return SegSegClass(SS_DISJOINT, false, false, false, false)
end

# ## Angle ordering at nodes (tangent-plane port of PolygonNodeTopology)
#
# Directions around an apex `n` live in the tangent plane at `n`. Pick a
# reference axis `r` (the coordinate axis least aligned with `n`, so `r ≁ ±n`);
# the tangent frame is `u = r - (r·n̂)n̂`, `v = n × r`, right-handed with `u×v =
# n̂`. A direction toward `p` has tangent coordinates `(p·u, p·v)`, and we only
# need their *signs* — both are determinant signs of `n, r, p`, exact for
# integer inputs and scale-corrected so the apex need not be unit:
#   sign(p·u) = sign((p·r)(n·n) - (r·n)(p·n)),   sign(p·v) = sign((n×r)·p).
# Feeding these to the planar quadrant scheme, with the same-quadrant tiebreak
# `rk_orient(m, n, q, p) = sign((n×q)·p)` (already the tangent-plane CCW sign),
# reproduces PolygonNodeTopology exactly.

# Coordinate axis least aligned with `n3` (smallest |component|, first-index
# tiebreak — matches `argmin`), as a unit vector of `n3`'s element type.
@inline function _ref_axis(n3)
    ax, ay, az = abs(n3[1]), abs(n3[2]), abs(n3[3])
    o = one(ax); z = zero(ax)
    if ax <= ay && ax <= az
        return (o, z, z)
    elseif ay <= az
        return (z, o, z)
    else
        return (z, z, o)
    end
end

# JTS quadrant of the direction toward `P3` around apex `n3` with reference
# `r3`: NE=0, NW=1, SW=2, SE=3, axis directions on the `>= 0` side.
@inline function _sph_quadrant3(n3, r3, P3)
    nn = _dot3(n3, n3); nr = _dot3(n3, r3); pn = _dot3(P3, n3); pr = _dot3(P3, r3)
    su = pr * nn - nr * pn               # sign of P·u
    sv = _dot3(_cross3(n3, r3), P3)      # sign of P·v
    (su == 0 && sv == 0) &&
        throw(ArgumentError("cannot compute the quadrant of a zero-length direction"))
    if su >= 0
        return sv >= 0 ? 0 : 3
    else
        return sv >= 0 ? 1 : 2
    end
end

function rk_quadrant(::Spherical, origin, p)
    n3 = _tup3(origin)
    return _sph_quadrant3(n3, _ref_axis(n3), _tup3(p))
end

# compareAngle around an explicit apex direction `n3` (a vec3 tuple): the
# crossing-apex slow path, where `n3` is the *constructed* crossing direction
# and so must be compared with explicit determinant signs (not ExactPredicates,
# which needs Float64 vertices). Mirrors `_compare_angle`: quadrant first, then
# the orient tiebreak `sign((n×q)·p)`.
function _sph_compare_around(bt, n3, p, q)
    P = _vec3(bt, p); Q = _vec3(bt, q)
    r3 = _ref_axis(n3)
    qp = _sph_quadrant3(n3, r3, P)
    qq = _sph_quadrant3(n3, r3, Q)
    qp > qq && return 1
    qp < qq && return -1
    o = _dot3(_cross3(n3, Q), P)
    return o > 0 ? 1 : (o < 0 ? -1 : 0)
end

# The crossing direction (the sphere point where the two arcs of a crossing
# node meet): ±(na×nb), the candidate strictly interior to both minor arcs.
function _sph_crossing_dir(bt, node::NodeKey)
    A0 = _vec3(bt, node.pt); A1 = _vec3(bt, node.a1)
    B0 = _vec3(bt, node.b0); B1 = _vec3(bt, node.b1)
    na = _cross3(A0, A1); nb = _cross3(B0, B1)
    d = _cross3(na, nb)
    (_strictly_in_arc3(d, A0, A1, na) && _strictly_in_arc3(d, B0, B1, nb)) && return d
    return _neg3(d)
end

function rk_compare_edge_dir(m::Spherical, node::NodeKey, p, q; exact)
    node.is_crossing || return _compare_angle(m, node.pt, p, q; exact)
    # Crossing apex: unlike the plane, the tangent direction apex→x is not
    # parallel to opp(x)→x, so the planar endpoint substitution does not carry
    # over. Compare around the (exact, on-arc) crossing direction instead —
    # the slow path, only on crossing-node edge ordering.
    bt = booltype(exact)
    return _sph_compare_around(bt, _sph_crossing_dir(bt, node), p, q)
end

# ## rk_nodes_coincide (exact slow path)
#
# Whether two node keys denote the same sphere point. The point of a vertex node
# is its coordinate direction; of a crossing node, the on-arc crossing direction
# `±(na×nb)`. Two directions denote the same sphere point iff they are parallel
# (cross product zero) and point into the same hemisphere (positive dot) — `-d`
# is the antipodal point, a different node. Exact via `Rational{BigInt}` (the
# `True()` branch of `_vec3`), mirroring the planar D3 rational slow path.
@inline _exact_node_dir(bt, k::NodeKey) =
    k.is_crossing ? _sph_crossing_dir(bt, k) : _vec3(bt, k.pt)

function rk_nodes_coincide(::Spherical, k1::NodeKey, k2::NodeKey; exact)
    k1 == k2 && return true
    bt = booltype(exact)
    d1 = _exact_node_dir(bt, k1); d2 = _exact_node_dir(bt, k2)
    return _iszero3(_cross3(d1, d2)) && _dot3(d1, d2) > 0
end

# ## Ring orientation

#=
Spherical method of `_ring_is_ccw` (relate_geometry.jl — the port of JTS
`Orientation.isCCW` used by `_orient_ring`). The planar extreme-vertex cap
algorithm assumes a coordinate plane: its y-extreme vertex pick and flat-cap
`del_x` tiebreak are meaningless on xyz points (a ring symmetric about the
equator has two exactly-equal extreme-y vertices and reads CW in *both*
windings). On the sphere the ring is CCW iff its signed area (Girard fan sum,
area.jl) is positive — iff the region on its left is the enclosed one, the
one smaller than a hemisphere. This is the sole place the engine resolves
which of the two ring-bounded regions a ring means; `_orient_ring` (edge-side
topology), `rk_point_in_ring`, and `rk_interaction_bounds` all inherit it, so
they agree by construction.

`exact` is accepted for signature parity but unused: the sign only selects
the region convention, and a near-zero sum means the ring splits the sphere
into equal halves — intrinsically ambiguous, not a float artifact.
=#
function _ring_is_ccw(::Spherical, ring::Vector; exact)
    n = length(ring)
    n > 1 && ring[end] == ring[1] && (n -= 1)
    n < 3 && return false
    # renormalize: the Girard quadrant split assumes unit vectors, and the
    # conformance suite feeds exact-integer non-unit rings
    apex = _girard_fan_apex(ring, n)
    total = 0.0
    prev = rk_normalize_usp(ring[n])
    for i in 1:n
        cur = rk_normalize_usp(ring[i])
        total += _spherical_triangle_area(Girard(), apex, prev, cur)
        prev = cur
    end
    return total > 0
end

#=
Fan apex for the Girard sum. The signed-area fan is apex-invariant in exact
math, but a fan triangle with an (even nearly) antipodal apex–vertex pair
has no well-defined connecting geodesic — its Girard excess degenerates to
`atan(≈0, ≈0)` — and corrupts the sum. Ingest only rejects antipodal
*edges*: a fan chord to a non-adjacent vertex can still be antipodal, and
post-`AntipodalEdgeSplit` rings carry exactly such vertex pairs (the split
keeps both endpoints of the offending edge). Pick the first vertex — then
the first edge midpoint — with no vertex within ~1e-9 of its antipode. A
ring defeating every candidate pairs its whole vertex set with antipodes;
fall back to the first vertex, which is no worse than the sum being taken
at all (such a ring splits the sphere near-evenly, where the sign is
intrinsically ambiguous).
=#
function _girard_fan_apex(ring, n)
    for i in 1:n
        apex = rk_normalize_usp(ring[i])
        _clean_fan_apex(apex, ring, n) && return apex
    end
    for i in 1:n
        mid = rk_normalize_usp(ring[i]) + rk_normalize_usp(ring[mod1(i + 1, n)])
        norm(mid) < 1e-9 && continue   # near-antipodal edge: unstable midpoint
        apex = UnitSphericalPoint(normalize(mid))
        _clean_fan_apex(apex, ring, n) && return apex
    end
    return rk_normalize_usp(ring[1])
end

_clean_fan_apex(apex, ring, n) =
    all(i -> (rk_normalize_usp(ring[i]) ⋅ apex) > -1 + 1e-9, 1:n)

# ## rk_point_in_ring (anchor-retry crossing parity, winding-independent)

# Whether the two minor arcs (p0,p1) and (q0,q1) cross properly (interior to
# both). The great circles meet at ±d, d = (p0×p1)×(q0×q1); a proper crossing is
# one of ±d strictly interior to both arcs (the spike's `arcs_cross_properly`).
function _arcs_cross_properly(bt, p0, p1, q0, q1)
    P0 = _vec3(bt, p0); P1 = _vec3(bt, p1); Q0 = _vec3(bt, q0); Q1 = _vec3(bt, q1)
    na = _cross3(P0, P1); nb = _cross3(Q0, Q1)
    d = _cross3(na, nb)
    _iszero3(d) && return false
    (_strictly_in_arc3(d, P0, P1, na) && _strictly_in_arc3(d, Q0, Q1, nb)) && return true
    nd = _neg3(d)
    return _strictly_in_arc3(nd, P0, P1, na) && _strictly_in_arc3(nd, Q0, Q1, nb)
end

# Ring vertices as spherical kernel points. A 3D ring is already in kernel
# coordinates (e.g. the conformance suite's exact integer USP rings) and is read
# verbatim — renormalizing would perturb the exact orient the boundary test
# relies on. A 2D (lon/lat) ring — the engine's ingested polygon — is converted
# to unit xyz.
_ring_kernel_pts(ring) = _ring_kernel_pts(booltype(GI.is3d(GI.getpoint(ring, 1))), ring)
_ring_kernel_pts(::True, ring) = _node_points(ring)
_ring_kernel_pts(::False, ring) = _ring_usp(ring)

# Location of `p` relative to the area enclosed by `ring` — the kernel
# contract (kernel.jl) is winding-independent, like the planar ray-crossing
# parity: real-world rings arrive in either winding (Natural Earth ships
# shapefile-convention CW shells), and `_locate_point_in_polygonal` passes
# them unoriented. Boundary first (exact arc membership), then the shared
# `spherical_ring_contains` (which reports the region on the ring's *left*)
# with this kernel's predicates injected — `rk_orient` for sides,
# `_arcs_cross_properly` for transversality — so the parity decision is as
# exact as the predicates; the left region is the interior iff the ring is
# canonically CCW (`_ring_is_ccw` above, the same bit `_orient_ring` feeds
# the edge-side topology). All anchors degenerate (unreachable for a
# non-degenerate ring and an off-boundary point) is refused, not answered
# wrong.
rk_point_in_ring(m::Spherical, p, ring; exact) =
    rk_point_in_ring(m, p, SphericalKernelRing(m, ring; exact); exact)

"""
    SphericalKernelRing(m::Spherical, ring; exact)

The cached kernel-space form of one ring: the converted
`UnitSphericalPoint` vertex vector (`pts` — the boundary edge walk), its
deduped open form (`ded`/`n` — the parity walk; aliases `pts` when the
ring has no repeated vertices), and the ring's Girard orientation bit
(`_ring_is_ccw`, the same bit edge topology and interaction bounds use).

`rk_point_in_ring` re-derived all of this from lon/lat on every query —
vertex conversion alone was ~60% of a prepared spherical point query. The
point-in-area locators (indexed_point_in_area.jl) convert each ring once
and query on this form (Layer 1 of the 2026-07-14 spherical-indexed-locator
design).

Repeated consecutive vertices are dropped from the parity walk (real rings
carry them — NE 110m North Korea's sliver is `[A, A, B, A]`; JTS removes
them at ingest, but this path receives the raw ring): a retraced edge lies
exactly under the anchor midpoint and breaks the parity count. After dedup
a ring with fewer than 3 distinct vertices bounds no area.
"""
struct SphericalKernelRing
    pts::Vector{UnitSphericalPoint{Float64}}
    ded::Vector{UnitSphericalPoint{Float64}}
    n::Int
    is_ccw::Bool
end

function SphericalKernelRing(m::Spherical, ring; exact)
    pts = _ring_kernel_pts(ring)
    n = length(pts)
    n > 1 && pts[end] == pts[1] && (n -= 1)
    ded, n = _drop_repeated_ring_pts(pts, n)
    is_ccw = n >= 3 && _ring_is_ccw(m, ded; exact)
    return SphericalKernelRing(pts, ded, n, is_ccw)
end

# Type-stable functors for the predicates injected into
# `spherical_ring_contains` (Layer 3 of the spherical-indexed-locator
# design): the anonymous closures they replace were rebuilt per call and
# heap-boxed, costing an allocation and dynamic dispatch per predicate call
# on the point-in-area hot path. The injectable-predicate design of
# `spherical_ring_contains` is unchanged.
struct _RKOrient{M <: Spherical, E} <: Function
    m::M
    exact::E
end
(f::_RKOrient)(a, b, c) = rk_orient(f.m, a, b, c; exact = f.exact)

struct _RKProperCrossing{BT} <: Function
    bt::BT
end
(f::_RKProperCrossing)(q, mid, a, b) = _arcs_cross_properly(f.bt, q, mid, a, b) ? 1 : 0

function rk_point_in_ring(m::Spherical, p, kr::SphericalKernelRing; exact)
    pts = kr.pts
    @inbounds for i in 1:length(pts)-1
        rk_point_on_segment(m, p, pts[i], pts[i+1]; exact) && return LOC_BOUNDARY
    end
    #-- rings are closed regardless of a repeated last point (the kernel
    #-- contract), so an implicitly closed ring's closing edge is boundary
    #-- too — the same edge set the longitude-interval index walks
    if length(pts) > 1 && pts[end] != pts[1] &&
            rk_point_on_segment(m, p, pts[end], pts[1]; exact)
        return LOC_BOUNDARY
    end
    kr.n < 3 && return LOC_EXTERIOR
    inside = spherical_ring_contains(kr.ded, kr.n, p;
        orient = _RKOrient(m, exact),
        on_arc = Returns(false),   # boundary classified exactly above
        proper_crossing = _RKProperCrossing(booltype(exact)))
    inside === nothing && _throw_degenerate_point_in_ring(p)
    return inside == kr.is_ccw ? LOC_INTERIOR : LOC_EXTERIOR
end

@noinline _throw_degenerate_point_in_ring(p) = throw(ArgumentError(
    "rk_point_in_ring: every anchor edge of the ring is degenerate with " *
    "respect to the query point $(_tup3(p)) — the ring is degenerate at " *
    "this point"))

# Interaction bounds on the sphere: the shared substrate
# (`spherical_arc_extent` per edge, `_spherical_region_extent` for area
# interiors) over kernel-converted points, so box and ingested vertices
# agree bit-for-bit. Rings are dim-1 linework here (JTS semantics), not S2
# regions — a CW hole must not become a complement region. Boxes get a few
# ulps of padding so a vertex from another conversion path still prunes as
# interacting.
rk_interaction_bounds(m::Spherical, geom) =
    _pad_bounds(_sph_interaction_extent(m, GI.trait(geom), geom))

_sph_interaction_extent(m::Spherical, ::GI.AbstractPointTrait, geom) =
    GI.extent(_spherical_kernel_point(geom))
function _sph_interaction_extent(m::Spherical, ::GI.AbstractCurveTrait, geom)
    n = GI.npoint(geom)
    prev = _spherical_kernel_point(GI.getpoint(geom, 1))
    # seeding with pts[1]'s box covers the degenerate n == 1 curve; it is
    # absorbed by the first edge box otherwise
    ext = spherical_arc_extent(prev, prev)
    for i in 2:n
        cur = _spherical_kernel_point(GI.getpoint(geom, i))
        ext = Extents.union(ext, spherical_arc_extent(prev, cur))
        prev = cur
    end
    return ext
end
function _sph_interaction_extent(m::Spherical, ::GI.AbstractPolygonTrait, geom)
    # region box of the exterior ring: edge arc extents plus enclosed-axis
    # widening, on the same converted points the engine ingests.
    # `_spherical_region_extent` bounds the region on the ring's left, so
    # orient the shell canonically CCW first — a CW-wound input (shapefile
    # convention) would otherwise bound the complement, under-covering an
    # enclosed pole (`exact` is unused by the spherical `_ring_is_ccw`)
    pts = _orient_ring(m, _ring_usp(GI.getexterior(geom)), false; exact = True())
    ext = _spherical_region_extent(pts)
    # a valid polygon's holes lie inside that region — but JTS's element
    # envelope also covers a stray hole outside the shell, and extraction
    # relies on that to keep the element alive
    # (see `_extract_segment_strings_from_atomic!`)
    for hole in GI.gethole(geom)
        GI.isempty(hole) && continue
        ext = Extents.union(ext, _sph_interaction_extent(m, GI.trait(hole), hole))
    end
    return ext
end
function _sph_interaction_extent(m::Spherical, ::GI.AbstractGeometryTrait, geom)
    ext = nothing
    for g in GI.getgeom(geom)
        GI.isempty(g) && continue
        e = _sph_interaction_extent(m, GI.trait(g), g)
        ext = ext === nothing ? e : Extents.union(ext, e)
    end
    return ext
end

# Converted (kernel-ingest: unit, signed-zero) vertices of a ring/curve.
_ring_usp(ring) = [_spherical_kernel_point(p) for p in GI.getpoint(ring)]

# `pts[1:n]` (implied closure) with repeated consecutive vertices removed,
# copying only when one exists; wraparound repeats included.
function _drop_repeated_ring_pts(pts, n)
    has_dup = false
    for i in 1:n
        if pts[i] == pts[mod1(i + 1, n)]
            has_dup = true
            break
        end
    end
    has_dup || return pts, n
    ded = empty(pts)
    sizehint!(ded, n)
    for i in 1:n
        (isempty(ded) || ded[end] != pts[i]) && push!(ded, pts[i])
    end
    length(ded) > 1 && ded[end] == ded[1] && pop!(ded)
    return ded, length(ded)
end

_pad_bounds(::Nothing) = nothing
_pad_bounds(ext) = Extents.Extent(
    X = _widen(ext.X...), Y = _widen(ext.Y...), Z = _widen(ext.Z...))
