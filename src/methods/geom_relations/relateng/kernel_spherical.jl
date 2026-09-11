# NOTE: This functionality is experimental and may change at any time.

# # Spherical RelateKernel
#
#=
Implement `RelateKernel` on `UnitSphericalPoint{Float64}`. Exact predicates use floating-point
filters followed by `exact_spherical_orient` or rational arithmetic on xyz components.

This file shares the `GeometryOps` module imports. The `UnitSphericalPoint` rebuild hook is
defined in `kernel.jl`.
=#

# xyz tuple of a 3D point, for the ExactPredicates / Rational{BigInt} paths.
@inline _tup3(u) = (GI.x(u), GI.y(u), GI.z(u))

# ## rk_orient

# Return the sign of `(a×b)·c`. The exact path uses `exact_spherical_orient`, whose filter
# bound scales with point separation. The approximate path uses the plain triple product.
rk_orient(::Spherical, a, b, c; exact) = _rk_orient(booltype(exact), a, b, c)
@inline function _rk_orient(::True, a, b, c)
    # Equal vectors make the triple product exactly zero; skip predicate evaluation.
    (_usp_eq(a, b) || _usp_eq(a, c) || _usp_eq(b, c)) && return 0
    return UnitSpherical.exact_spherical_orient(_tup3(a), _tup3(b), _tup3(c))
end
@inline _rk_orient(::False, a, b, c) = cross(a, b) ⋅ c

# ## Exact-aware 3-vector arithmetic
#
# `_vec3` selects `Rational{BigInt}` for exact composite predicates or `Float64` otherwise.
# Tuple cross and dot products share both paths.
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

# A point lies on `[q0, q1]` if it is coplanar and within the minor-arc span. For `p = α q0 +
# β q1`, require α, β ≥ 0.
#
# With `n = q0×q1`, `sign(β) = sign((q0×p)·n)` and `sign(α) = sign((p×q1)·n)`. These
# determinant signs also apply to non-unit inputs.
function rk_point_on_segment(m::Spherical, p, q0, q1; exact)
    rk_orient(m, q0, q1, p; exact) == 0 || return false
    return _on_arc_span(booltype(exact), p, q0, q1)
end
# Filter span signs in `Float64`; use the rational authority when the bound cannot certify the
# result.
@inline function _on_arc_span(bt::True, p, q0, q1)
    r = _on_arc_span_filter(p, q0, q1)
    r === nothing || return r
    return _on_arc_span_authority(bt, p, q0, q1)
end
# Evaluate approximate span membership in `Float64` without filtering.
@inline _on_arc_span(bt::False, p, q0, q1) = _on_arc_span_authority(bt, p, q0, q1)

@inline function _on_arc_span_authority(bt, p, q0, q1)
    P = _vec3(bt, p); Q0 = _vec3(bt, q0); Q1 = _vec3(bt, q1)
    n = _cross3(Q0, Q1)
    if _iszero3(n)
        # Parallel endpoints define a degenerate arc. Accept only endpoint directions; a zero
        # normal would otherwise make every span test pass.
        return (_iszero3(_cross3(P, Q0)) && _dot3(P, Q0) > 0) ||
               (_iszero3(_cross3(P, Q1)) && _dot3(P, Q1) > 0)
    end
    return _dot3(_cross3(Q0, P), n) >= 0 && _dot3(_cross3(P, Q1), n) >= 0
end

# For degree-4 span determinants `(u × v) · n`, running-error analysis gives about
# `9μ·Σᵢ|wᵢ_terms|·|nᵢ_terms|`, where `μ = ½eps`.
#
# Use `16μ` to cover higher-order terms and rounding in the magnitude sum. Homogeneity makes
# the bound valid for non-unit inputs.
const _SPAN_ERR_C = 16 * (eps(Float64) / 2)

# Certify the span decision only when `|value|` exceeds its error bound. A certified negative
# sign returns `false`; two certified positive signs return `true`.
#
# Return `nothing` for unresolved signs, including degenerate normals. Exact endpoint matches
# are handled separately.
@inline function _on_arc_span_filter(p, q0, q1)
    # An endpoint belongs to its closed arc. Test identity before filtering its zero span
    # determinant.
    (_usp_eq(p, q0) || _usp_eq(p, q1)) && return true
    x0 = GI.x(q0); y0 = GI.y(q0); z0 = GI.z(q0)
    x1 = GI.x(q1); y1 = GI.y(q1); z1 = GI.z(q1)
    xp = GI.x(p);  yp = GI.y(p);  zp = GI.z(p)
    # n = q0 × q1, with per-component abs-magnitude sums Nᵢ
    n1 = y0*z1 - z0*y1;  N1 = abs(y0*z1) + abs(z0*y1)
    n2 = z0*x1 - x0*z1;  N2 = abs(z0*x1) + abs(x0*z1)
    n3 = x0*y1 - y0*x1;  N3 = abs(x0*y1) + abs(y0*x1)
    # s1 = (q0 × p) · n
    w1 = y0*zp - z0*yp;  W1 = abs(y0*zp) + abs(z0*yp)
    w2 = z0*xp - x0*zp;  W2 = abs(z0*xp) + abs(x0*zp)
    w3 = x0*yp - y0*xp;  W3 = abs(x0*yp) + abs(y0*xp)
    s1 = w1*n1 + w2*n2 + w3*n3
    e1 = _SPAN_ERR_C * (W1*N1 + W2*N2 + W3*N3)
    # s2 = (p × q1) · n
    v1 = yp*z1 - zp*y1;  V1 = abs(yp*z1) + abs(zp*y1)
    v2 = zp*x1 - xp*z1;  V2 = abs(zp*x1) + abs(xp*z1)
    v3 = xp*y1 - yp*x1;  V3 = abs(xp*y1) + abs(yp*x1)
    s2 = v1*n1 + v2*n2 + v3*n3
    e2 = _SPAN_ERR_C * (V1*N1 + V2*N2 + V3*N3)
    (s1 < -e1 || s2 < -e2) && return false     # one span factor certainly < 0
    (s1 > e1 && s2 > e2) && return true          # both span factors certainly > 0
    return nothing                               # Near a boundary; escalate.
end

# ## Ingest and interaction bounds

# Normalize off-unit input. Skip division when the squared norm is within 4 ULPs of 1 to keep
# ingestion idempotent. Repeated normalization can otherwise alternate between rounded values.
@inline function rk_normalize_usp(u)
    s = u[1] * u[1] + u[2] * u[2] + u[3] * u[3]
    abs(s - one(s)) <= 4 * eps(one(s)) && return UnitSphericalPoint(u)
    return UnitSphericalPoint(normalize(u))
end

# Convert lon/lat to unit xyz, or normalize 3D input, then canonicalize signed zeros.
# Ingestion and extent calculation share this conversion for bit-identical vertices.
@inline function _spherical_kernel_point(p)
    u = GI.is3d(p) ?
        UnitSphericalPoint(Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p))) :
        UnitSphereFromGeographic()((Float64(GI.x(p)), Float64(GI.y(p))))
    return _node_point(rk_normalize_usp(u))
end

# Ingest spherical points through `_spherical_kernel_point`, matching extent conversion.
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
# Great circles meet at `±d`, where `d = (a0×a1)×(b0×b1)`. A proper crossing requires the same
# candidate strictly inside both minor arcs. Endpoint incidences use exact arc membership.
#
# ### Float-fast path: the four-orient reduction
#
# With `na = a0×a1`, `nb = b0×b1`, and `d = na×nb`, BAC–CAB gives:
#
#     (a0×d)·na = (a0·nb)|na|²         (d×a1)·na = −(a1·nb)|na|²
#     (b0×d)·nb = −(b0·na)|nb|²        (d×b1)·nb =  (b1·na)|nb|²
#
# Writing `[u,v,w] = u·(v×w)`, for nonzero `d`:
#
#     _strictly_in_arc3(d , a,·) ⟺ [b0,b1,a0]>0 ∧ [b0,b1,a1]<0
#     _strictly_in_arc3(d , b,·) ⟺ [a0,a1,b0]<0 ∧ [a0,a1,b1]>0
#
# Negating `d` flips all four signs, yielding the S2 crossing pattern. All-zero orientations
# use the exact same-circle/degenerate classifier.
#
# Zero-length arcs cannot satisfy the strict pattern because two orientation signs coincide.
rk_classify_intersection(m::Spherical, a0, a1, b0, b1; exact) =
    _rk_classify_intersection(booltype(exact), m, a0, a1, b0, b1)

# Classify from four filtered exact orientation signs, using the reduction above.
function _rk_classify_intersection(bt::True, m, a0, a1, b0, b1)
    # Reuse `sABi = sign[a0,a1,bi]` and `sBAi = sign[b0,b1,ai]` for arc membership.
    sAB0 = rk_orient(m, a0, a1, b0; exact = bt)
    sAB1 = rk_orient(m, a0, a1, b1; exact = bt)
    sBA0 = rk_orient(m, b0, b1, a0; exact = bt)
    sBA1 = rk_orient(m, b0, b1, a1; exact = bt)
    # Arc membership: `p` on the arc iff coplanar (orient 0) and within the
    # minor-arc span. `_on_arc_span` (float-filtered) only fires on the coplanar
    # `== 0` cases, i.e. shared vertices / T-junctions.
    a0_on_b = sBA0 == 0 && _on_arc_span(bt, a0, b0, b1)
    a1_on_b = sBA1 == 0 && _on_arc_span(bt, a1, b0, b1)
    b0_on_a = sAB0 == 0 && _on_arc_span(bt, b0, a0, a1)
    b1_on_a = sAB1 == 0 && _on_arc_span(bt, b1, a0, a1)
    if sAB0 == 0 && sAB1 == 0 && sBA0 == 0 && sBA1 == 0
        # Use the exact authority for a shared great circle or degenerate arc.
        return _sph_classify(bt, a0, a1, b0, b1, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    end
    # d ≠ 0 (proven exactly by a nonzero orient). Endpoint incidence ⇒ touch.
    if a0_on_b || a1_on_b || b0_on_a || b1_on_a
        return SegSegClass(SS_TOUCH, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    end
    # Proper crossing ⟺ +d or −d strictly interior to both arcs, i.e. the
    # four-orient near-crossing pattern (equal to `_strictly_in_arc3(±d,…)`).
    proper = _proper_crossing_from_orients(sAB0, sAB1, sBA0, sBA1)
    return proper ? SegSegClass(SS_PROPER, false, false, false, false) :
                    SegSegClass(SS_DISJOINT, false, false, false, false)
end

# The approximate path uses direct `Float64` candidate directions and span tests.
function _rk_classify_intersection(bt::False, m, a0, a1, b0, b1)
    a0_on_b = rk_point_on_segment(m, a0, b0, b1; exact = bt)
    a1_on_b = rk_point_on_segment(m, a1, b0, b1; exact = bt)
    b0_on_a = rk_point_on_segment(m, b0, a0, a1; exact = bt)
    b1_on_a = rk_point_on_segment(m, b1, a0, a1; exact = bt)
    return _sph_classify(bt, a0, a1, b0, b1, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
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
# Choose the coordinate axis `r` least aligned with apex `n`. In the tangent frame, `u = r -
# (r·n̂)n̂` and `v = n × r`.
#
# Use `sign((p·r)(n·n) - (r·n)(p·n))` and `sign((n×r)·p)` for the quadrant. These signs are
# scale-invariant, so `n` need not be unit.
#
# Within a quadrant, `rk_orient(m, n, q, p)` supplies the CCW ordering, as in
# PolygonNodeTopology.

# Return the unit coordinate axis least aligned with `n3`, using its element type.
# Choose the first index on ties, as in `argmin`.
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
    su, sv = _sph_tangent_signs(n3, r3, P3)
    (su == 0 && sv == 0) &&
        throw(ArgumentError("cannot compute the quadrant of a zero-length direction"))
    if su >= 0
        return sv >= 0 ? 0 : 3
    else
        return sv >= 0 ? 1 : 2
    end
end

#=
Compute exact tangent-coordinate signs for the node quadrant. Cancellation near frame axes can
corrupt edge ordering and wedge tests.

Filter both signs and escalate unresolved cases. Since `r3` is a coordinate axis, `P·r`,
`n·r`, and `n×r` are exact; only `n·n` and `P·n` contribute rounding error.
=#
@inline _sph_tangent_signs(n3, r3, P3) = _sph_tangent_signs_exact(n3, r3, P3)

@inline function _sph_tangent_signs(n3::NTuple{3, Float64}, r3::NTuple{3, Float64},
        P3::NTuple{3, Float64})
    nn = _dot3(n3, n3); nr = _dot3(n3, r3); pn = _dot3(P3, n3); pr = _dot3(P3, r3)
    su = pr * nn - nr * pn
    w = _cross3(n3, r3)
    sv = _dot3(w, P3)
    #-- `nn` and the |Pᵢnᵢ| sum are the only inexact inputs; 8·eps covers the
    #-- three dot-product roundings, the two products and the subtraction
    e = 8 * eps(Float64)
    pn_mag = abs(P3[1] * n3[1]) + abs(P3[2] * n3[2]) + abs(P3[3] * n3[3])
    su_bound = e * (abs(pr) * nn + abs(nr) * pn_mag)
    sv_bound = e * (abs(w[1] * P3[1]) + abs(w[2] * P3[2]) + abs(w[3] * P3[3]))
    (abs(su) > su_bound && abs(sv) > sv_bound) && return (sign(su), sign(sv))
    return _sph_tangent_signs_exact(_rat3(n3), _rat3(r3), _rat3(P3))
end

@inline _rat3(t) = (Rational{BigInt}(t[1]), Rational{BigInt}(t[2]), Rational{BigInt}(t[3]))

@inline function _sph_tangent_signs_exact(n3, r3, P3)
    nn = _dot3(n3, n3); nr = _dot3(n3, r3); pn = _dot3(P3, n3); pr = _dot3(P3, r3)
    return (sign(pr * nn - nr * pn), sign(_dot3(_cross3(n3, r3), P3)))
end

function rk_quadrant(::Spherical, origin, p)
    n3 = _tup3(origin)
    return _sph_quadrant3(n3, _ref_axis(n3), _tup3(p))
end

# Order directions around explicit crossing apex `n3`: compare quadrants, then
# `sign((n×q)·p)`. Use determinant arithmetic because the constructed apex need not be
# `Float64`.
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

#=
Return the crossing direction `±(na×nb)` strictly inside both minor arcs. `bt` selects
coordinate arithmetic; exact orientation signs always select the antipodal candidate.

With `[u,v,w] = u·(v×w)`, `na = a0×a1`, and `nb = b0×b1`, BAC–CAB gives:

    d = na×nb = [b0,b1,a0]·a1 − [b0,b1,a1]·a0 = [a0,a1,b1]·b0 − [a0,a1,b0]·b1

For a proper crossing, `d` is a positive combination of each arc's endpoints iff `[b0,b1,a0] >
0`. The four signs have pattern `(+,−,−,+)` or its negation, so any one selects the candidate.
=#
function _sph_crossing_dir(bt, node::NodeKey)
    A0 = _vec3(bt, node.pt); A1 = _vec3(bt, node.a1)
    B0 = _vec3(bt, node.b0); B1 = _vec3(bt, node.b1)
    d = _cross3(_cross3(A0, A1), _cross3(B0, B1))
    return _crossing_dir_is_positive(node) ? d : _neg3(d)
end

# Compute the on-arc direction exactly, then scale and round it to a unit-sphere point for
# location queries. The rounded representative is not an exact node key.
function _crossing_locate_point(::Spherical, key::NodeKey)
    d = _sph_crossing_dir(True(), key)
    scale = max(abs(d[1]), abs(d[2]), abs(d[3]))
    x = Float64(d[1] / scale); y = Float64(d[2] / scale); z = Float64(d[3] / scale)
    s = sqrt(x * x + y * y + z * z)
    return UnitSphericalPoint(x / s, y / s, z / s)
end

# Select `+(na×nb)` from the first nonzero exact orientation. Proper crossings have four
# nonzero signs. An all-zero set implies `d == 0`, so the choice is immaterial.
function _crossing_dir_is_positive(node::NodeKey)
    s = _rk_orient(True(), node.b0, node.b1, node.pt); s != 0 && return s > 0
    s = _rk_orient(True(), node.b0, node.b1, node.a1); s != 0 && return s < 0
    s = _rk_orient(True(), node.pt, node.a1, node.b0); s != 0 && return s < 0
    s = _rk_orient(True(), node.pt, node.a1, node.b1); s != 0 && return s > 0
    return true
end

function rk_compare_edge_dir(m::Spherical, node::NodeKey, p, q; exact)
    node.is_crossing || return _compare_angle(m, node.pt, p, q; exact)
    # Order crossing-node edges around the exact crossing direction. Spherical tangent
    # directions do not permit planar endpoint substitution.
    bt = booltype(exact)
    return _sph_compare_around(bt, _sph_crossing_dir(bt, node), p, q)
end

# ## rk_nodes_coincide (exact slow path)
#
# Two nodes coincide iff their directions have zero cross product and positive dot product.
# Use stored vertex directions or exact on-arc crossing directions in rational arithmetic.
@inline _exact_node_dir(bt, k::NodeKey) =
    k.is_crossing ? _sph_crossing_dir(bt, k) : _vec3(bt, k.pt)

function rk_nodes_coincide(::Spherical, k1::NodeKey, k2::NodeKey; exact)
    k1 == k2 && return true
    bt = booltype(exact)
    d1 = _exact_node_dir(bt, k1); d2 = _exact_node_dir(bt, k2)
    return _iszero3(_cross3(d1, d2)) && _dot3(d1, d2) > 0
end

# ## Node ordering along an arc (design §2.5)
#
# Along minor arc `s0 → s1`, `da` precedes `db` iff `(da × db) · N > 0`, where `N = s0 × s1`.
#
# Certify the floating-point sign with a bound including errors in both directions and `N`;
# otherwise recompute with rational directions.

# Bound each cross-product component error by `2u·(|aⱼbₖ| + |aₖbⱼ|)`. Sum these component
# bounds to bound the Euclidean error norm.
@inline function _cross3_err(a, b)
    p1 = a[2]*b[3]; q1 = a[3]*b[2]
    p2 = a[3]*b[1]; q2 = a[1]*b[3]
    p3 = a[1]*b[2]; q3 = a[2]*b[1]
    e = eps(Float64) * ((abs(p1) + abs(q1)) + (abs(p2) + abs(q2)) + (abs(p3) + abs(q3)))
    return ((p1 - q1, p2 - q2, p3 - q3), e)
end

#=
Return the floating-point node direction and its relative error bound. Stored vertex
directions have zero conversion error.

For crossings, `d = ±(na × nb)` with `na = a0×a1` and `nb = b0×b1`. Short edges increase
relative normal error; near-parallel normals amplify it by `1/sin θ`.

Writing `Δ` for absolute errors:

Both fall out of one bound. Writing `Δ` for absolute errors,

    |Δd| ≤ |Δna||nb| + |na||Δnb| + (rounding of na×nb)

and `|d| = |na||nb| sin θ`, so the relative error of the direction is

    ε = (|Δna||nb| + |na||Δnb| + Δ(na×nb)) / |d|

The bound grows as edges shorten or approach tangency. A zero `|d|` gives `Inf` and forces
exact evaluation.
=#
@inline function _float_node_dir_err(k::NodeKey)
    k.is_crossing || return (_vec3(False(), k.pt), 0.0)
    A0 = _vec3(False(), k.pt); A1 = _vec3(False(), k.a1)
    B0 = _vec3(False(), k.b0); B1 = _vec3(False(), k.b1)
    na, e_na = _cross3_err(A0, A1)
    nb, e_nb = _cross3_err(B0, B1)
    d,  e_d  = _cross3_err(na, nb)
    dn = sqrt(_dot3(d, d))
    #-- the sign choice is exact on both paths (`_crossing_dir_is_positive`) and
    #-- does not touch the magnitudes
    dir = _crossing_dir_is_positive(k) ? d : _neg3(d)
    dn == 0 && return (dir, Inf)
    rel = (e_na * sqrt(_dot3(nb, nb)) + e_nb * sqrt(_dot3(na, na)) + e_d) / dn
    return (dir, rel)
end

#=
Bound the discriminant error using relative errors in both node directions and the arc normal,
plus cross-and-dot rounding.

Since `|disc| ≤ |da||db||N| = mag`, `rel ≥ 1` makes the bound at least `|disc|` and forces
escalation. This also covers near-tangent crossings without a separate threshold.
=#
function rk_compare_along_segment(m::Spherical, s0, s1, na::NodeKey, nb::NodeKey; exact)
    S0 = _vec3(False(), s0); S1 = _vec3(False(), s1)
    N, e_N = _cross3_err(S0, S1)
    da, rel_a = _float_node_dir_err(na)
    db, rel_b = _float_node_dir_err(nb)
    disc = _dot3(_cross3(da, db), N)
    nN = sqrt(_dot3(N, N))
    mag = sqrt(_dot3(da, da) * _dot3(db, db)) * nN
    #-- 16 ulp covers the cross-then-dot that forms `disc` (6 products, 5 adds)
    rel = rel_a + rel_b + (nN == 0 ? Inf : e_N / nN) + 16 * eps(Float64)
    tol = rel * mag
    #-- `tol` is NaN only if `rel` is `Inf` and `mag` is 0, i.e. fully degenerate;
    #-- either way the comparison is false and the exact path decides
    abs(disc) > tol && return disc > 0 ? -1 : 1
    #-- exact fallback (lazy): rational directions, exact for Float64 inputs.
    Se0 = _vec3(True(), s0); Se1 = _vec3(True(), s1); Ne = _cross3(Se0, Se1)
    ea = _exact_node_dir(True(), na); eb = _exact_node_dir(True(), nb)
    o = _dot3(_cross3(ea, eb), Ne)
    return o > 0 ? -1 : (o < 0 ? 1 : 0)
end

# ## Ring orientation

#=
Determine spherical winding from geodesic curvature using S2 `GetCurvature`
(`s2loop_measures.cc`). Gauss–Bonnet gives left-side area `2π − curvature`; nonnegative
curvature selects at most a hemisphere.

Only adjacent vertex pairs enter each turn, so non-adjacent antipodal vertices remain valid.
`_orient_ring`, `rk_point_in_ring`, and interaction bounds share this region choice.

Allow S2's curvature error bound of `11.25ε` per vertex. Exact hemispheres therefore count as
CCW in both windings.

`exact` is unused: turn signs always use exact orientation. Normalize vertices because
`robust_cross_product` requires unit input.
=#
function _ring_is_ccw(::Spherical, ring::Vector; exact)
    loop = _prune_loop_degeneracies([rk_normalize_usp(p) for p in ring])
    n = length(loop)
    n < 3 && return false   # bounds no area (JTS convention for flat rings)
    return _spherical_loop_curvature(loop) >= -(11.25 * eps(Float64) * n)
end

#=
Port S2 `PruneDegeneracies` (`s2loop_measures.cc`): remove repeated vertices (`AA → A`) and
retraced edges (`ABA → A`), including across closure.

Remaining vertices have distinct, non-retracing neighbors. A fully degenerate loop returns
fewer than three vertices.
=#
function _prune_loop_degeneracies(pts::Vector, same = ==)
    vertices = empty(pts)
    sizehint!(vertices, length(pts))
    for v in pts
        if !isempty(vertices)
            same(v, vertices[end]) && continue                       # AA → A
            if length(vertices) >= 2 && same(v, vertices[end - 1])   # ABA → A
                pop!(vertices)
                continue
            end
        end
        push!(vertices, v)
    end
    length(vertices) > 1 && same(vertices[1], vertices[end]) && pop!(vertices)
    m = length(vertices)
    m < 3 && return vertices
    # Remove retraced edge pairs across closure. A non-degenerate portion remains, so this
    # stops before consuming the loop.
    k = 0
    while same(vertices[k + 2], vertices[m - k]) || same(vertices[k + 1], vertices[m - k - 1])
        k += 1
    end
    return k == 0 ? vertices : vertices[(k + 1):(m - k)]
end

# Port S2 `TurnAngle` (`s2measures.cc`): return the turn at `b` along `a → b → c`, positive
# for CCW. Use robust edge normals for magnitude and exact orientation for sign.
function _sph_turn_angle(a, b, c)
    angle = _usp_angle(robust_cross_product(a, b), robust_cross_product(b, c))
    return _rk_orient(True(), a, b, c) > 0 ? angle : -angle
end

# S2 `Vector3.Angle`: atan2(|u×v|, u·v), stable near both parallel and
# antiparallel (acos of the dot is not).
_usp_angle(u, v) = atan(norm(cross(u, v)), u ⋅ v)

#=
Port S2 `GetCurvature` (`s2loop_measures.cc`): sum turns in canonical order with Kahan
compensation, then restore the stored direction's sign.

Positive curvature means left-side area `2π − curvature` is smaller than a hemisphere.
=#
function _spherical_loop_curvature(loop)
    n = length(loop)
    i, dir = _canonical_loop_order(loop)
    at(k) = loop[mod1(k, n)]
    total = _sph_turn_angle(at(i - dir), at(i), at(i + dir))
    compensation = 0.0
    for _ in 1:(n - 1)
        i += dir
        angle = _sph_turn_angle(at(i - dir), at(i), at(i + dir)) + compensation
        old_total = total
        total += angle
        compensation = (old_total - total) + angle
    end
    return dir * (total + compensation)
end

#=
Port S2 `GetCanonicalLoopOrder` (`s2loop_measures.cc`): minimize the vertex sequence
lexicographically over rotations in both directions.

Canonical summation and sign restoration make curvature invariant under rotation and negated
under reversal.
=#
function _canonical_loop_order(loop)
    n = length(loop)
    min_indices = [1]
    for i in 2:n
        if _tup3(loop[i]) <= _tup3(loop[min_indices[1]])
            _tup3(loop[i]) < _tup3(loop[min_indices[1]]) && empty!(min_indices)
            push!(min_indices, i)
        end
    end
    best = (min_indices[1], 1)
    for i in min_indices
        _loop_order_less((i, 1), best, loop) && (best = (i, 1))
        _loop_order_less((i, -1), best, loop) && (best = (i, -1))
    end
    return best
end

# Port of S2 `IsOrderLess`: whether traversal `o1` yields a lexicographically
# smaller vertex sequence than `o2` (both start at the same minimal vertex).
function _loop_order_less(o1, o2, loop)
    o1 == o2 && return false
    n = length(loop)
    (i1, d1) = o1
    (i2, d2) = o2
    for _ in 1:(n - 1)
        i1 += d1; i2 += d2
        p1 = _tup3(loop[mod1(i1, n)]); p2 = _tup3(loop[mod1(i2, n)])
        p1 < p2 && return true
        p1 > p2 && return false
    end
    return false
end

#=
## The same orientation test over EXACT vertex directions

Overlay vertices can have exact rational directions closer than one floating-point ULP.
Rounding can merge vertices or change edge directions and turn signs.

Compute curvature directly from exact directions. Cross products and turn signs remain
rational; only angle magnitudes use `BigFloat` square roots and `atan`.

Normalization is unnecessary: each turn is invariant under positive scaling of its vertex
directions.
=#
function _ring_is_ccw_dirs(dirs::Vector)
    loop = _prune_loop_degeneracies(dirs, _dirs_same_point)
    n = length(loop)
    n < 3 && return false   # bounds no area (JTS convention for flat rings)
    return _spherical_loop_curvature_exact(loop) >= -(11.25 * eps(Float64) * n)
end

# Whether two exact directions denote the same point of the sphere. `==` would
# be wrong: the directions are unnormalized, so the same point has infinitely
# many representations, all positive multiples of each other.
@inline _dirs_same_point(a, b) = _cross3(a, b) == (zero(a[1]), zero(a[1]), zero(a[1])) &&
    _dot3(a, b) > 0

# Sum exact-direction curvature in canonical order with a `BigFloat` accumulator,
# then restore the stored direction's sign.
function _spherical_loop_curvature_exact(loop)
    n = length(loop)
    i, dir = _canonical_loop_order(loop)
    at(k) = loop[mod1(k, n)]
    total = zero(BigFloat)
    for _ in 1:n
        total += _sph_turn_angle_exact(at(i - dir), at(i), at(i + dir))
        i += dir
    end
    return dir * total
end

# `_sph_turn_angle` over exact directions: exact sign, `BigFloat` magnitude.
function _sph_turn_angle_exact(a, b, c)
    angle = _usp_angle_exact(_cross3(a, b), _cross3(b, c))
    return _dot3(_cross3(a, b), c) > 0 ? angle : -angle
end

# `_usp_angle` with both arguments exact: `atan2(|u×v|, u·v)`, the cross and dot
# taken in rational arithmetic and only the final `atan` (and the `sqrt` feeding
# it) in `BigFloat`.
function _usp_angle_exact(u, v)
    w = _cross3(u, v)
    return atan(sqrt(BigFloat(_dot3(w, w))), BigFloat(_dot3(u, v)))
end

#=
For `Spherical(; oriented = true)`, stored winding defines the region as in S2 `InitOriented`
(`s2polygon.h`). Polygon interior lies left of every ring.

Shells denote their left region; holes denote their right cavity. Reversing a ring denotes the
complement, allowing regions larger than a hemisphere.
=#
_ring_interior_on_left(m::Spherical, pts::Vector, is_hole::Bool; exact) =
    m.oriented ? !is_hole : _ring_is_ccw(m, pts; exact)

# ## rk_point_in_ring (definitional-exterior crossing parity, winding-independent)

# Test whether both minor arcs contain the same antipodal circle intersection in their
# interiors. Use four exact orientations; straddling alone can select opposite antipodes.
@inline _proper_crossing_from_orients(sa, sb, sq, sm) =
    (sq > 0 && sm < 0 && sa < 0 && sb > 0) ||
    (sq < 0 && sm > 0 && sa > 0 && sb < 0)

function _arcs_cross_properly(bt::True, p0, p1, q0, q1)
    sa = _rk_orient(bt, p0, p1, q0)
    sb = _rk_orient(bt, p0, p1, q1)
    sq = _rk_orient(bt, q0, q1, p0)
    sm = _rk_orient(bt, q0, q1, p1)
    return _proper_crossing_from_orients(sa, sb, sq, sm)
end

# Preserve the approximate cross/dot formulation: its rounding differs from
# the four-orient reduction near degeneracies.
function _arcs_cross_properly(bt::False, p0, p1, q0, q1)
    P0 = _vec3(bt, p0); P1 = _vec3(bt, p1); Q0 = _vec3(bt, q0); Q1 = _vec3(bt, q1)
    na = _cross3(P0, P1); nb = _cross3(Q0, Q1)
    d = _cross3(na, nb)
    _iszero3(d) && return false
    (_strictly_in_arc3(d, P0, P1, na) && _strictly_in_arc3(d, Q0, Q1, nb)) && return true
    nd = _neg3(d)
    return _strictly_in_arc3(nd, P0, P1, na) && _strictly_in_arc3(nd, Q0, Q1, nb)
end

# Read 3D ring vertices unchanged to preserve exact boundary predicates. Convert 2D lon/lat
# vertices to unit xyz.
_ring_kernel_pts(ring) = _ring_kernel_pts(booltype(GI.is3d(GI.getpoint(ring, 1))), ring)
_ring_kernel_pts(::True, ring) = _node_points(ring)
_ring_kernel_pts(::False, ring) = _ring_usp(ring)

# Locate `p` by exact boundary membership, then crossing parity with kernel orientation and
# transversality predicates.
#
# Default mode uses winding-independent `spherical_ring_encloses`, anchored at the antipode of
# vertex mass. Self-intersections use even-odd semantics.
#
# If the anchor or test arc is degenerate, fall back to `spherical_ring_contains` with
# `_ring_interior_on_left`. Near-hemisphere winding uses the curvature tolerance.
#
# Oriented mode uses the stored winding and `is_hole` role. It selects the appropriate side of
# `spherical_ring_contains`; if every anchor is degenerate, reject the query.
rk_point_in_ring(m::Spherical, p, ring; exact, is_hole::Bool = false) =
    rk_point_in_ring(m, p, SphericalKernelRing(m, ring; exact, is_hole); exact)

"""
    SphericalKernelRing(m::Spherical, ring; exact, is_hole = false)

Cache converted vertices (`pts`), the deduplicated open parity walk (`ded`, `n`), the
interior-side bit, and the exterior parity anchor. `ded` aliases `pts` when no deduplication
is needed.

The interior-side bit matches edge topology and bounds. The anchor is `nothing` in oriented
mode or for degenerate vertex mass; the latter uses the wedge fallback.

Point-in-area locators convert each ring once. Remove consecutive duplicates before parity
traversal; fewer than three distinct vertices enclose no area.
"""
struct SphericalKernelRing
    pts::Vector{UnitSphericalPoint{Float64}}
    ded::Vector{UnitSphericalPoint{Float64}}
    n::Int
    interior_on_left::Bool
    anchor::Union{Nothing, UnitSphericalPoint{Float64}}
end

function SphericalKernelRing(m::Spherical, ring; exact, is_hole::Bool = false)
    pts = _ring_kernel_pts(ring)
    n = length(pts)
    n > 1 && pts[end] == pts[1] && (n -= 1)
    ded, n = _drop_repeated_ring_pts(pts, n)
    interior_on_left = n >= 3 && _ring_interior_on_left(m, ded, is_hole; exact)
    anchor = (!m.oriented && n >= 3) ? spherical_exterior_anchor(ded, n) : nothing
    return SphericalKernelRing(pts, ded, n, interior_on_left, anchor)
end

# Typed predicate functors keep `spherical_ring_contains` calls specialized without per-call
# closure allocation.
struct _RKOrient{M <: Spherical, E} <: Function
    m::M
    exact::E
end
(f::_RKOrient)(a, b, c) = rk_orient(f.m, a, b, c; exact = f.exact)

struct _RKProperCrossing{BT} <: Function
    bt::BT
end
(f::_RKProperCrossing)(q, mid, a, b) = _arcs_cross_properly(f.bt, q, mid, a, b) ? 1 : 0

# Reuse signs only when BOTH injected predicates use the exact kernel. Custom
# orientations (or approximate ones) need not have produced these exact signs.
@inline function UnitSpherical._proper_crossing_with_orients(
        f::_RKProperCrossing{True}, orient::_RKOrient{M, True},
        q, mid, a, b, sa, sb, sq, sm) where {M}
    return _proper_crossing_from_orients(sa, sb, sq, sm) ? 1 : 0
end

# Exact span test for the anchor walk's vertex-grazing resolution
# (`_anchor_crossing_parity`): whether `p`, already known to lie on the
# great circle of `(a, b)`, lies on the closed minor arc.
struct _RKOnTestArc{BT} <: Function
    bt::BT
end
(f::_RKOnTestArc)(p, a, b) = _on_arc_span(f.bt, p, a, b)

function rk_point_in_ring(m::Spherical, p, kr::SphericalKernelRing; exact)
    pts = kr.pts
    @inbounds for i in 1:length(pts)-1
        rk_point_on_segment(m, p, pts[i], pts[i+1]; exact) && return LOC_BOUNDARY
    end
    #-- Include the implicit closing edge in boundary tests.
    if length(pts) > 1 && pts[end] != pts[1] &&
            rk_point_on_segment(m, p, pts[end], pts[1]; exact)
        return LOC_BOUNDARY
    end
    kr.n < 3 && return LOC_EXTERIOR
    orient = _RKOrient(m, booltype(exact))
    on_arc = Returns(false)   # boundary classified exactly above
    proper_crossing = _RKProperCrossing(booltype(exact))
    if !m.oriented
        #-- Use exterior-anchor parity. An undecidable result uses the wedge fallback.
        enc = spherical_ring_encloses(kr.ded, kr.n, p; anchor = kr.anchor,
            orient, on_arc, proper_crossing,
            on_test_arc = _RKOnTestArc(booltype(exact)))
        enc === nothing || return enc ? LOC_INTERIOR : LOC_EXTERIOR
    end
    inside = spherical_ring_contains(kr.ded, kr.n, p; orient, on_arc, proper_crossing)
    inside === nothing && _throw_degenerate_point_in_ring(p)
    return inside == kr.interior_on_left ? LOC_INTERIOR : LOC_EXTERIOR
end

@noinline _throw_degenerate_point_in_ring(p) = throw(ArgumentError(
    "rk_point_in_ring: every anchor edge of the ring is degenerate with " *
    "respect to the query point $(_tup3(p)) — the ring is degenerate at " *
    "this point"))

# Compute edge and region bounds from kernel-converted points. Rings use linework bounds; pad
# boxes for rounding between conversion paths.
#
# Trust stored 3D extents as kernel-space bounds, including the full region for polygons.
# Cached extents avoid repeated region-bound calculation.
function rk_interaction_bounds(m::Spherical, geom)
    _reusable_stored_extent(m, geom) && return geom.extent
    return _pad_bounds(_sph_interaction_extent(m, GI.trait(geom), geom))
end

_sph_interaction_extent(m::Spherical, ::GI.AbstractPointTrait, geom) =
    GI.extent(_spherical_kernel_point(geom))
# A preparation-only pass: validate each represented edge while accumulating
# its bounds. An optional local point vector lets polygon preparation reuse the
# shell's ingest conversion without changing persistent coordinate ownership.
function _sph_validated_curve_extent(geom, points)
    n = GI.npoint(geom)
    prev = _spherical_kernel_point(GI.getpoint(geom, 1))
    points === nothing || push!(points, prev)
    ext = spherical_arc_extent(prev, prev)
    for i in 2:n
        cur = _spherical_kernel_point(GI.getpoint(geom, i))
        _exactly_antipodal(prev, cur) && _throw_antipodal_edge(prev, cur)
        points === nothing || push!(points, cur)
        ext = Extents.union(ext, spherical_arc_extent(prev, cur))
        prev = cur
    end
    return ext
end

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
    # Orient the shell's denoted region to the left before `_spherical_region_extent`.
    # Oriented shells already satisfy this convention; complement shells may bound nearly the
    # whole sphere.
    pts = _orient_ring(m, _ring_usp(GI.getexterior(geom)), false, false; exact = True())
    ext = _spherical_region_extent(pts)
    # Include stray holes outside the shell, matching JTS element envelopes and preserving
    # them during segment extraction.
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

# Collect converted vertices into a plain `Vector`. Comprehensions can inherit StaticArrays
# axes and return `SizedVector`, which point-list consumers reject.
function _ring_usp(ring)
    pts = Vector{UnitSphericalPoint{Float64}}()
    sizehint!(pts, GI.npoint(ring))
    for p in GI.getpoint(ring)
        push!(pts, _spherical_kernel_point(p))
    end
    return pts
end

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
