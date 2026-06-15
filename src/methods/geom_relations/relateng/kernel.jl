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

    rk_quadrant(m, origin, p)::Int

Quadrant of the direction from `origin` to `p`, in the JTS `Quadrant`
convention: `0` = NE, `1` = NW, `2` = SW, `3` = SE, numbered CCW from the
positive X-axis, with axis directions belonging to the `dx >= 0` /
`dy >= 0` side. Throws `ArgumentError` for a zero-length direction.

    rk_compare_edge_dir(m, node::NodeKey, p, q; exact)::Int

Compare the angles of the edge directions `node → p` and `node → q` around
the (possibly symbolic) apex `node`: negative / zero / positive as the
direction toward `p` has angle less than / equal to / greater than the
direction toward `q`, angles increasing CCW from the positive X-axis at the
apex (port of JTS `PolygonNodeTopology.compareAngle` with a `NodeKey` apex).
For vertex nodes the apex coordinate is exact and the port is direct. For
crossing nodes, directions which are endpoints of the node's defining
segments are compared exactly from the original endpoints, never from a
constructed apex coordinate; foreign directions (incident edges of a D3
coincidence-merged node, from other segment pairs crossing at the same
point) are compared exactly around the rational apex (slow path).

    rk_crossing_dirs_ccw(m, a0, a1, b0, b1; exact)

CCW cyclic order of the four half-edge directions incident to the proper
crossing of `(a0, a1)` × `(b0, b1)`, as a 4-tuple of the original segment
endpoints starting from `a1`. Derived from a single orientation sign; no
crossing coordinate is constructed.

    rk_is_crossing(m, node, a0, a1, b0, b1; exact)::Bool

Whether the rings entering vertex node `node` along segments `a0–node–a1`
and `b0–node–b1` cross at the node (port of JTS
`PolygonNodeTopology.isCrossing`). If any segment is collinear with another,
returns `false`. Vertex-node apexes only: proper crossings are crossings by
construction and are short-circuited by the caller
(JTS `TopologyComputer.updateAreaAreaCross`).

    rk_is_interior_segment(m, node, a0, a1, b; exact)::Bool

Whether the segment `node → b` lies in the interior of the ring corner
`a0–node–a1`, the ring interior being on the right of the corner (i.e. a CW
shell or CCW hole); port of JTS `PolygonNodeTopology.isInteriorSegment`.
The test segment must not be collinear with the corner segments.
Vertex-node apexes only.
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

# A node point is the engine's canonical, signed-zero-normalized representation
# of a coordinate. Planar points stay 2-tuples (byte-identical to the original
# implementation, so `NodeKey` bytes are unchanged); 3D points (e.g. a spherical
# `UnitSphericalPoint`) keep all three components and their concrete type.
@inline _node_point(p) = _node_point(booltype(GI.is3d(p)), p)
@inline _node_point(::False, p) = (_pos_zero(GI.x(p)), _pos_zero(GI.y(p)))
@inline function _node_point(::True, p)
    return _rebuild_point(p, _pos_zero(GI.x(p)), _pos_zero(GI.y(p)), _pos_zero(GI.z(p)))
end
# Default rebuild: a plain 3-tuple. Concrete point types that must survive node
# construction add their own method (e.g. `UnitSphericalPoint` in
# `kernel_spherical.jl`).
@inline _rebuild_point(::Any, x, y, z) = (x, y, z)

# Collect a curve's coordinates as node points into a plain `Vector`. (A
# typed comprehension over `GI.getpoint` is not enough: for geometries backed
# by StaticArrays — e.g. the `extent_to_polygon` output — the iterator has
# static axes, so `collect` returns a `SizedVector`, which downstream code
# expecting `Vector` point lists rejects.)
function _node_points(geom)
    p1 = GI.getpoint(geom, 1)
    pts = Vector{typeof(_node_point(p1))}()
    sizehint!(pts, GI.npoint(geom))
    for p in GI.getpoint(geom)
        push!(pts, _node_point(p))
    end
    return pts
end

"Node key of a vertex node: keyed exactly by its coordinate."
function vertex_node(pt)
    p = _node_point(pt)
    return NodeKey(false, p, p, p, p)
end

"""
    crossing_node(a0, a1, b0, b1)::NodeKey

Node key of the proper crossing of segments `(a0, a1)` and `(b0, b1)`.
Canonicalize: each segment ordered lexicographically by (x, y); segments
ordered lexicographically by their endpoint tuples — so any
order/orientation of the same pair produces an identical key.

Only construct crossing keys for properly crossing segments (`SS_PROPER`
from `rk_classify_intersection`): the exact rational slow path in
`rk_nodes_coincide` divides by the segments' direction cross product, which
is nonzero precisely when the crossing is proper.
"""
function crossing_node(a0, a1, b0, b1)
    a0, a1 = _canonical_segment(_node_point(a0), _node_point(a1))
    b0, b1 = _canonical_segment(_node_point(b0), _node_point(b1))
    if (GI.x(b0), GI.y(b0), GI.x(b1), GI.y(b1)) < (GI.x(a0), GI.y(a0), GI.x(a1), GI.y(a1))
        a0, a1, b0, b1 = b0, b1, a0, a1
    end
    return NodeKey(true, a0, a1, b0, b1)
end

# Order a segment's endpoints lexicographically by (x, y).
_canonical_segment(p, q) = (GI.x(p), GI.y(p)) <= (GI.x(q), GI.y(q)) ? (p, q) : (q, p)

# Whether `p` lies within the coordinate bounding box of segment `(q0, q1)`.
# Valid as an on-segment test only when `p` is already known collinear with
# `(q0, q1)`; shared by manifolds whose segments are coordinate-monotone.
@inline function _collinear_between(p, q0, q1)
    (min(GI.x(q0), GI.x(q1)) <= GI.x(p) <= max(GI.x(q0), GI.x(q1))) &&
    (min(GI.y(q0), GI.y(q1)) <= GI.y(p) <= max(GI.y(q0), GI.y(q1)))
end
