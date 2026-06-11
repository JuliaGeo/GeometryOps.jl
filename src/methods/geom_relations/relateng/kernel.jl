# # RelateKernel API contract

#=
The geometry layer of RelateNG (design doc D1/D2): every coordinate-level
question the topology layer may ask, answered with exact predicates and no
constructed coordinates. Each function takes the manifold as its first
argument; the planar implementations live in `kernel_planar.jl`. A future
`Spherical` kernel implements the same functions and must pass the same
conformance testset (Task 9).

All kernel functions are prefixed `rk_` (RelateKernel) and are internal —
nothing here is exported. Points are coordinate tuples (typically
`Tuple{Float64, Float64}`) obtained via `_tuple_point`; the `exact` flag is a
keyword taking `True()`/`False()` (GeometryOpsCore BoolsAsTypes), threaded
exactly like `Predicates.orient`.

The contract — what every manifold implementation must provide:

    rk_orient(m, a, b, c; exact)

Orientation of point `c` relative to the oriented segment `(a, b)`,
returned as a sign-valued number: `> 0` (counterclockwise, `c` to the
left), `< 0` (clockwise, `c` to the right), `== 0` (collinear, or
`a == b`). With `exact = True()` the sign must be correct even for
adversarial near-collinear inputs.

    rk_point_on_segment(m, p, q0, q1; exact)::Bool

Whether point `p` lies on the closed segment `[q0, q1]`, endpoints included.

    rk_point_in_ring(m, p, ring; exact)::Int8

Location of point `p` relative to the area enclosed by the closed `ring`
(a GeoInterface linestring/linearring, assumed closed regardless of a
repeated last point): one of `LOC_INTERIOR`, `LOC_BOUNDARY`, `LOC_EXTERIOR`.

    rk_interaction_bounds(m, geom)::Extents.Extent

The bounding region within which `geom` can interact with another geometry.
On the plane this is the ordinary extent; other manifolds may need to widen
it (e.g. great-circle edges bulge outside the coordinate box of their
endpoints).

    rk_bounds_disjoint(extA, extB)::Bool
    rk_bounds_covers(extA, extB)::Bool

Conservative interaction-bounds tests used for short-circuiting:
`rk_bounds_disjoint` must only return `true` when no interaction is possible;
`rk_bounds_covers` must only return `true` when `extA` covers `extB` in the
X/Y dimensions. These operate on the extents produced by
`rk_interaction_bounds` and are manifold-generic.

    rk_classify_intersection(m, a0, a1, b0, b1; exact)::SegSegClass

Combinatorial classification of the intersection of the closed segments
`(a0, a1)` and `(b0, b1)` (replaces JTS's `RobustLineIntersector`, design
D2). No intersection coordinate is ever constructed: a proper interior
crossing is reported purely symbolically as `SS_PROPER`, and all vertex
incidences are reported via the boolean `*_on_*` flags of the returned
`SegSegClass`, whose coordinates are exact input vertices. With
`exact = True()` the classification must be correct even for adversarial
near-collinear inputs.

    vertex_node(pt)::NodeKey
    crossing_node(a0, a1, b0, b1)::NodeKey

Manifold-generic constructors for symbolic node identities (design D2).
`vertex_node` keys a node by its exact coordinate; `crossing_node` keys a
proper-crossing node by the canonicalized defining segment pair, never by a
computed intersection coordinate. Keys constructed from the same vertex, or
from the same segment pair in any order/orientation, are `==` and hash
equal, so they can be used directly as `Dict` keys.

    rk_nodes_coincide(m, k1, k2; exact)::Bool

Whether two node keys denote the same point of the manifold. Same-kind keys
that are `==` trivially coincide; cross-kind coincidence (a vertex lying
exactly on a proper crossing, or two distinct crossings meeting at one
point) is decided exactly — on the plane via `Rational{BigInt}` arithmetic
(design D3). Only invoked on self-noding paths, so the slow path is
acceptable.
=#

# Symbolic segment-pair intersection classification (replaces RobustLineIntersector).
# JTS LineIntersector outcome mapping (for porting reference):
#   SS_DISJOINT  ↔ NO_INTERSECTION
#   SS_PROPER    ↔ POINT_INTERSECTION with isProper()
#   SS_TOUCH     ↔ POINT_INTERSECTION, not proper (incl. collinear abutment)
#   SS_COLLINEAR ↔ COLLINEAR_INTERSECTION
@enum SegSegKind::Int8 SS_DISJOINT SS_PROPER SS_TOUCH SS_COLLINEAR

"""
    SegSegClass

Combinatorial classification of the intersection of closed segments
(a0,a1) × (b0,b1). `kind` is `SS_PROPER` only for a crossing in both
segments' interiors (the node is *symbolic*: no coordinate exists for
it anywhere in the engine). All vertex incidences are reported via the
`*_on_*` flags, whose coordinates are exact input vertices.
"""
struct SegSegClass
    kind::SegSegKind
    a0_on_b::Bool
    a1_on_b::Bool
    b0_on_a::Bool
    b1_on_a::Bool
end

# Manifold-generic helpers

# Symbolic node identity (design D2). One concrete isbits key type for both
# node kinds so Dict{NodeKey{P}, ...} is type-stable. Equality and hashing
# are the default bit-pattern (egal) semantics for isbits structs; this is
# safe because the constructors normalize the only Float64 values whose
# numeric equality disagrees with bit equality: signed zeros (-0.0 → 0.0,
# via `x + zero(x)`, exact in IEEE arithmetic).
"""
    NodeKey{P}

Symbolic identity of a node (design D2). Vertex nodes key exactly by their
coordinate (`is_crossing == false`, all point fields equal to the vertex);
proper-crossing nodes key by their canonicalized defining segment pair
(`is_crossing == true`, fields `(pt, a1)` and `(b0, b1)` are the two
segments). No intersection coordinate is ever computed for the key.
Construct via [`vertex_node`](@ref) and [`crossing_node`](@ref).
"""
struct NodeKey{P}
    is_crossing::Bool
    pt::P          # vertex nodes: the coordinate. crossing nodes: canonical a0.
    a1::P
    b0::P
    b1::P
end

# Normalize signed zeros: -0.0 + 0.0 == +0.0 exactly, every other finite
# value is unchanged. Keeps bit-pattern key equality == coordinate equality.
@inline _pos_zero(x) = x + zero(x)
@inline _node_point(p) = (_pos_zero(GI.x(p)), _pos_zero(GI.y(p)))

"Node key of a vertex node: keyed exactly by its coordinate."
function vertex_node(pt)
    p = _node_point(pt)
    return NodeKey(false, p, p, p, p)
end

"""
    crossing_node(a0, a1, b0, b1)::NodeKey

Node key of the proper crossing of segments `(a0, a1)` and `(b0, b1)`.
Canonicalize: each segment ordered lexicographically by (x, y); segments
ordered by their first point — so any order/orientation of the same pair
produces an identical key.
"""
function crossing_node(a0, a1, b0, b1)
    a0, a1 = _seg_canon(_node_point(a0), _node_point(a1))
    b0, b1 = _seg_canon(_node_point(b0), _node_point(b1))
    if (GI.x(b0), GI.y(b0), GI.x(b1), GI.y(b1)) < (GI.x(a0), GI.y(a0), GI.x(a1), GI.y(a1))
        a0, a1, b0, b1 = b0, b1, a0, a1
    end
    return NodeKey(true, a0, a1, b0, b1)
end

# Order a segment's endpoints lexicographically by (x, y).
_seg_canon(p, q) = (GI.x(p), GI.y(p)) <= (GI.x(q), GI.y(q)) ? (p, q) : (q, p)

# Whether `p` lies within the coordinate bounding box of segment `(q0, q1)`.
# Valid as an on-segment test only when `p` is already known collinear with
# `(q0, q1)`; shared by manifolds whose segments are coordinate-monotone.
@inline function _collinear_between(p, q0, q1)
    (min(GI.x(q0), GI.x(q1)) <= GI.x(p) <= max(GI.x(q0), GI.x(q1))) &&
    (min(GI.y(q0), GI.y(q1)) <= GI.y(p) <= max(GI.y(q0), GI.y(q1)))
end
