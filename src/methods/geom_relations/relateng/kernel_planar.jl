# # Planar RelateKernel

#=
Planar implementation of the RelateKernel contract declared in `kernel.jl`.
Orientation goes through `Predicates.orient` (AdaptivePredicates when
`exact = True()`), point-in-ring reuses the existing Hao–Sun ray-crossing
machinery (`_point_filled_curve_orientation`), and bounds are plain extents.
No coordinates are ever constructed.
=#

rk_orient(::Planar, a, b, c; exact) = Predicates.orient(a, b, c; exact)

function rk_point_on_segment(m::Planar, p, q0, q1; exact)
    rk_orient(m, q0, q1, p; exact) == 0 || return false
    return _collinear_between(p, q0, q1)
end

function rk_point_in_ring(m::Planar, p, ring; exact)
    o = _point_filled_curve_orientation(m, p, ring; in = point_in, on = point_on, out = point_out, exact)
    o == point_in && return LOC_INTERIOR
    o == point_on && return LOC_BOUNDARY
    return LOC_EXTERIOR
end

#=
Planar interaction bounds are the plain GI extent: stored extents are
returned as-is (the wrapper tree built by `_relate_cache_extents`, or a
user's `GO.tuples(x; calc_extent = true)` input). The computed fallback,
however, does not go through `GI.calc_extent`, which makes two separate
closure-`extrema` passes over the points — on the ring-heavy extent-cache
pass that dominates unprepared point-in-area queries a single min/max sweep
is ~13× faster with identical results (same trait dispatch, same 2D/3D
bounds, same union-of-stored-extents semantics for GeometryCollections).
=#
function rk_interaction_bounds(m::Planar, geom)
    ex = GI.extent(geom; fallback = false)
    ex === nothing || return ex
    return _planar_sweep_extent(m, GI.trait(geom), geom)
end

function _planar_sweep_extent(::Planar, ::GI.AbstractPointTrait, p)
    x, y = GI.x(p), GI.y(p)
    GI.is3d(p) || return Extents.Extent(X = (x, x), Y = (y, y))
    z = GI.z(p)
    return Extents.Extent(X = (x, x), Y = (y, y), Z = (z, z))
end

#-- concrete GeometryCollections union their members' extents (reading
#-- stored ones), exactly as GI's `calc_extent` does
_planar_sweep_extent(m::Planar, ::GI.GeometryCollectionTrait, geom) =
    reduce(Extents.union, (rk_interaction_bounds(m, g) for g in GI.getgeom(geom)))

function _planar_sweep_extent(::Planar, ::GI.AbstractGeometryTrait, geom)
    itr = GI.getpoint(geom)
    st = iterate(itr)
    #-- empty geometry: defer to GI's own (throwing) computation path
    st === nothing && return GI.extent(geom, fallback = true)
    p, s = st
    xlo = xhi = GI.x(p)
    ylo = yhi = GI.y(p)
    if GI.is3d(geom)
        zlo = zhi = GI.z(p)
        while (st = iterate(itr, s)) !== nothing
            p, s = st
            x, y, z = GI.x(p), GI.y(p), GI.z(p)
            xlo = min(xlo, x); xhi = max(xhi, x)
            ylo = min(ylo, y); yhi = max(yhi, y)
            zlo = min(zlo, z); zhi = max(zhi, z)
        end
        return Extents.Extent(X = (xlo, xhi), Y = (ylo, yhi), Z = (zlo, zhi))
    end
    while (st = iterate(itr, s)) !== nothing
        p, s = st
        x, y = GI.x(p), GI.y(p)
        xlo = min(xlo, x); xhi = max(xhi, x)
        ylo = min(ylo, y); yhi = max(yhi, y)
    end
    return Extents.Extent(X = (xlo, xhi), Y = (ylo, yhi))
end

# Exact coordinate equality of two points.
_equals2(p, q) = GI.x(p) == GI.x(q) && GI.y(p) == GI.y(q)

function rk_classify_intersection(m::Planar, a0, a1, b0, b1; exact)
    oa0 = rk_orient(m, b0, b1, a0; exact)
    oa1 = rk_orient(m, b0, b1, a1; exact)
    ob0 = rk_orient(m, a0, a1, b0; exact)
    ob1 = rk_orient(m, a0, a1, b1; exact)
    # fully collinear configuration (handles zero-length segments too)
    if oa0 == 0 && oa1 == 0 && ob0 == 0 && ob1 == 0
        a0_on_b = _collinear_between(a0, b0, b1)
        a1_on_b = _collinear_between(a1, b0, b1)
        b0_on_a = _collinear_between(b0, a0, a1)
        b1_on_a = _collinear_between(b1, a0, a1)
        n_inc = a0_on_b + a1_on_b + b0_on_a + b1_on_a
        n_inc == 0 && return SegSegClass(SS_DISJOINT, false, false, false, false)
        # single shared endpoint counts twice (one endpoint of each on the other)
        shared_endpoint_only = n_inc == 2 &&
            ((a0_on_b || a1_on_b) && (b0_on_a || b1_on_a)) &&
            (_equals2(a0, b0) || _equals2(a0, b1) || _equals2(a1, b0) || _equals2(a1, b1))
        kind = shared_endpoint_only ? SS_TOUCH : SS_COLLINEAR
        # zero-length degenerate: a point on a segment is a touch, not an overlap
        if _equals2(a0, a1) || _equals2(b0, b1)
            kind = SS_TOUCH
        end
        return SegSegClass(kind, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    end
    a0_on_b = oa0 == 0 && _collinear_between(a0, b0, b1)
    a1_on_b = oa1 == 0 && _collinear_between(a1, b0, b1)
    b0_on_a = ob0 == 0 && _collinear_between(b0, a0, a1)
    b1_on_a = ob1 == 0 && _collinear_between(b1, a0, a1)
    if a0_on_b || a1_on_b || b0_on_a || b1_on_a
        return SegSegClass(SS_TOUCH, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    end
    if (oa0 > 0) != (oa1 > 0) && oa0 != 0 && oa1 != 0 &&
       (ob0 > 0) != (ob1 > 0) && ob0 != 0 && ob1 != 0
        return SegSegClass(SS_PROPER, false, false, false, false)
    end
    return SegSegClass(SS_DISJOINT, false, false, false, false)
end

# Edge ordering around nodes: port of JTS PolygonNodeTopology
# (algorithm/PolygonNodeTopology.java), with the apex generalized to a
# symbolic NodeKey. Vertex-node apexes are a direct port; crossing-node
# apexes are handled exactly via the original segment endpoints (see
# rk_compare_edge_dir) — no apex coordinate is ever constructed.

# Port of PolygonNodeTopology.quadrant / Quadrant.quadrant: NE=0, NW=1,
# SW=2, SE=3, numbered CCW from the positive X-axis; axis directions belong
# to the `dx >= 0` / `dy >= 0` side. Pure coordinate comparisons, exact.
function rk_quadrant(::Planar, origin, p)
    ox, oy = GI.x(origin), GI.y(origin)
    px, py = GI.x(p), GI.y(p)
    (px == ox && py == oy) &&
        throw(ArgumentError("cannot compute the quadrant of a zero-length direction"))
    if px >= ox
        return py >= oy ? 0 : 3   # NE : SE
    else
        return py >= oy ? 1 : 2   # NW : SW
    end
end

# Port of PolygonNodeTopology.compareAngle with a vertex apex: angles
# increase CCW from the positive X-axis; different quadrants decide the
# comparison, same-quadrant ties are resolved by orientation (P > Q if P is
# CCW of Q).
#
# This and the helpers below (`_is_angle_greater`, `_is_between`,
# `_compare_between`, `rk_crossing_dirs_ccw`, `rk_is_crossing`,
# `rk_is_interior_segment`) depend only on `rk_quadrant` and `rk_orient`, both
# manifold-dispatched, so they are manifold-generic (`m::Manifold`). The
# spherical kernel supplies a tangent-plane `rk_quadrant`; the same-quadrant
# orient tiebreak `rk_orient(m, origin, q, p)` is already the tangent-plane CCW
# sign there. Only `rk_quadrant` and the crossing-apex `rk_compare_edge_dir`
# are manifold-specific.
function _compare_angle(m::Manifold, origin, p, q; exact)
    quadrant_p = rk_quadrant(m, origin, p)
    quadrant_q = rk_quadrant(m, origin, q)
    quadrant_p > quadrant_q && return 1
    quadrant_p < quadrant_q && return -1
    # vectors are in the same quadrant: check relative orientation
    o = rk_orient(m, origin, q, p; exact)
    o > 0 && return 1
    o < 0 && return -1
    return 0
end

# Port of PolygonNodeTopology.isAngleGreater.
function _is_angle_greater(m::Manifold, origin, p, q; exact)
    quadrant_p = rk_quadrant(m, origin, p)
    quadrant_q = rk_quadrant(m, origin, q)
    quadrant_p > quadrant_q && return true
    quadrant_p < quadrant_q && return false
    # vectors are in the same quadrant: P > Q if it is CCW of Q
    return rk_orient(m, origin, q, p; exact) > 0
end

# Port of PolygonNodeTopology.isBetween: whether edge p is inside the arc
# from e0 to e1 (the arc not including the origin direction). Edges assumed
# distinct (non-collinear).
function _is_between(m::Manifold, origin, p, e0, e1; exact)
    _is_angle_greater(m, origin, p, e0; exact) || return false
    return !_is_angle_greater(m, origin, p, e1; exact)
end

# Port of PolygonNodeTopology.compareBetween: 1 if p is inside the arc from
# e0 to e1 (the arc not crossing the positive X-axis), -1 if outside, 0 if
# collinear with either edge.
function _compare_between(m::Manifold, origin, p, e0, e1; exact)
    comp0 = _compare_angle(m, origin, p, e0; exact)
    comp0 == 0 && return 0
    comp1 = _compare_angle(m, origin, p, e1; exact)
    comp1 == 0 && return 0
    (comp0 > 0 && comp1 < 0) && return 1
    return -1
end

# The opposite endpoint of incident endpoint `p` on its defining segment of
# crossing node `k`, or `nothing` if `p` is not one of the four endpoints
# (a direction from a foreign segment pair, on a D3 coincidence-merged node).
# Coordinate-equality matching is unambiguous because the four endpoints of
# a proper crossing are pairwise distinct.
function _crossing_opposite(k::NodeKey, p)
    _equals2(p, k.pt) && return k.a1
    _equals2(p, k.a1) && return k.pt
    _equals2(p, k.b0) && return k.b1
    _equals2(p, k.b1) && return k.b0
    return nothing
end

function rk_compare_edge_dir(m::Planar, node::NodeKey, p, q; exact)
    node.is_crossing || return _compare_angle(m, node.pt, p, q; exact)
    #=
    Crossing apex (needed by RelateEdge/NodeSection edge ordering, where the
    node may be a proper crossing): the directions to compare are normally
    among the four endpoints of the defining segments. Because the symbolic
    apex lies *strictly* inside both segments (SS_PROPER), for any incident
    endpoint `x` with opposite endpoint `opp(x)` on the same segment:
      - the vector apex → x is a positive multiple of opp(x) → x, so
        quadrant(apex, x) == quadrant(opp(x), x), and
      - the directed line apex → x equals the directed line opp(x) → x, so
        sign(orient(apex, x, y)) == sign(orient(opp(x), x, y)).
    Substituting these into compareAngle reproduces the Java comparison
    (anchored at the positive X-axis of the apex) exactly, without ever
    constructing the apex coordinate.
    =#
    popp = _crossing_opposite(node, p)
    qopp = _crossing_opposite(node, q)
    if popp !== nothing && qopp !== nothing
        quadrant_p = rk_quadrant(m, popp, p)
        quadrant_q = rk_quadrant(m, qopp, q)
        quadrant_p > quadrant_q && return 1
        quadrant_p < quadrant_q && return -1
        # same quadrant: orient(apex, q, p) has the sign of orient(opp(q), q, p).
        # Zero only when p == q (distinct incident endpoints in the same quadrant
        # are never collinear through the apex of a proper crossing).
        o = rk_orient(m, qopp, q, p; exact)
        o > 0 && return 1
        o < 0 && return -1
        return 0
    end
    #=
    A direction point from a foreign segment pair: the node is a D3
    coincidence-merged node (TopologyComputer's self-noding merge pass),
    whose incident edges come from several segment pairs crossing at the
    same point. The endpoint substitution above does not apply, so compare
    around the exact rational apex instead (slow path; only reachable on
    the rare self-noding merge path).
    =#
    apex = _exact_crossing_point(node.pt, node.a1, node.b0, node.b1)
    return _compare_angle_exact(apex, p, q)
end

"""
CCW cyclic order of the four half-edge directions incident to the
proper crossing of (a0,a1) × (b0,b1), starting from a1. Since the
crossing is proper, b0/b1 are strictly on opposite sides of line(a0,a1):
if b1 is to the left, CCW order is (a1, b1, a0, b0), else (a1, b0, a0, b1).
"""
function rk_crossing_dirs_ccw(m::Manifold, a0, a1, b0, b1; exact)
    if rk_orient(m, a0, a1, b1; exact) > 0
        return (a1, b1, a0, b0)
    else
        return (a1, b0, a0, b1)
    end
end

# Port of PolygonNodeTopology.isCrossing, apex = vertex node coordinate.
# Crossing-node apexes are rejected: a proper crossing is a crossing by
# construction, and the only caller (TopologyComputer.updateAreaAreaCross)
# short-circuits proper intersections before asking.
function rk_is_crossing(m::Manifold, node::NodeKey, a0, a1, b0, b1; exact)
    node.is_crossing &&
        throw(ArgumentError("rk_is_crossing requires a vertex-node apex; proper crossings cross by construction"))
    nodept = node.pt
    a_lo, a_hi = a0, a1
    if _is_angle_greater(m, nodept, a_lo, a_hi; exact)
        a_lo, a_hi = a1, a0
    end
    # Find positions of b0 and b1. The edges cross if the positions are
    # different. If any edge is collinear they are reported as not crossing.
    comp_between0 = _compare_between(m, nodept, b0, a_lo, a_hi; exact)
    comp_between0 == 0 && return false
    comp_between1 = _compare_between(m, nodept, b1, a_lo, a_hi; exact)
    comp_between1 == 0 && return false
    return comp_between0 != comp_between1
end

# Port of PolygonNodeTopology.isInteriorSegment, apex = vertex node
# coordinate: whether segment node→b lies in the interior of the ring corner
# a0–node–a1 (ring interior on the right, i.e. a CW shell or CCW hole).
function rk_is_interior_segment(m::Manifold, node::NodeKey, a0, a1, b; exact)
    node.is_crossing &&
        throw(ArgumentError("rk_is_interior_segment requires a vertex-node apex"))
    nodept = node.pt
    a_lo, a_hi = a0, a1
    is_interior_between = true
    if _is_angle_greater(m, nodept, a_lo, a_hi; exact)
        a_lo, a_hi = a1, a0
        is_interior_between = false
    end
    is_between = _is_between(m, nodept, b, a_lo, a_hi; exact)
    return (is_between && is_interior_between) || (!is_between && !is_interior_between)
end

# Node coincidence, rational slow path (design D3). Float64 values are
# dyadic rationals, so Rational{BigInt} conversion and arithmetic are exact.

# Precondition: the segments cross *properly* (SS_PROPER), so their
# direction vectors are non-parallel and the denominator below is nonzero.
# crossing_node keys are only ever constructed for proper crossings.
"Exact intersection point of two properly crossing segments, as rationals."
function _exact_crossing_point(a0, a1, b0, b1)
    R = Rational{BigInt}
    ax0, ay0 = R(GI.x(a0)), R(GI.y(a0)); ax1, ay1 = R(GI.x(a1)), R(GI.y(a1))
    bx0, by0 = R(GI.x(b0)), R(GI.y(b0)); bx1, by1 = R(GI.x(b1)), R(GI.y(b1))
    dax, day = ax1 - ax0, ay1 - ay0
    dbx, dby = bx1 - bx0, by1 - by0
    denom = dax * dby - day * dbx          # nonzero for a proper crossing
    t = ((bx0 - ax0) * dby - (by0 - ay0) * dbx) // denom
    return (ax0 + t * dax, ay0 + t * day)
end

_exact_node_point(k::NodeKey) = k.is_crossing ?
    _exact_crossing_point(k.pt, k.a1, k.b0, k.b1) :
    (Rational{BigInt}(GI.x(k.pt)), Rational{BigInt}(GI.y(k.pt)))

function rk_nodes_coincide(::Planar, k1::NodeKey, k2::NodeKey; exact)
    k1 == k2 && return true
    # Slow path (design D3, follow-up F1): exact rational comparison.
    return _exact_node_point(k1) == _exact_node_point(k2)
end

#=
`_compare_angle` with an exact rational apex: the slow path of
`rk_compare_edge_dir` for direction points incident to a D3
coincidence-merged crossing node which are not endpoints of the node's
defining segments. All arithmetic is over `Rational{BigInt}` (Float64
inputs convert exactly), so the comparison is exact. Mirrors
`_compare_angle` (quadrant first, then `orient(origin, q, p)`).

Precondition: the direction points `p` and `q` must differ from the apex
`origin` — the quadrant of a zero vector is undefined. This is unreachable
in practice because section direction points are adjacent ring/line
vertices, which are distinct from the node.
=#
function _compare_angle_exact(origin, p, q)
    R = Rational{BigInt}
    ox, oy = origin
    px, py = R(GI.x(p)), R(GI.y(p))
    qx, qy = R(GI.x(q)), R(GI.y(q))
    _quad(x, y) = x >= ox ? (y >= oy ? 0 : 3) : (y >= oy ? 1 : 2)
    quadrant_p = _quad(px, py)
    quadrant_q = _quad(qx, qy)
    quadrant_p > quadrant_q && return 1
    quadrant_p < quadrant_q && return -1
    #-- same quadrant: orient(origin, q, p) as an exact rational cross product
    o = (qx - ox) * (py - oy) - (qy - oy) * (px - ox)
    o > 0 && return 1
    o < 0 && return -1
    return 0
end
