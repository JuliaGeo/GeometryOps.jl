# # Result ring linking — maximal rings, minimal-ring split, edge rings
#
# Phase 2b of the OverlayNG port (design doc §3). Ports two JTS files kept
# together because they form one pipeline over the marked result-area edges:
#   - `MaximalEdgeRing.java`  → maximal-ring linking + the minimal-ring split
#     at self-touching nodes (`_MaxEdgeRing` + `_link_result_area_max_ring_at_node!`).
#   - `OverlayEdgeRing.java`  → one minimal result ring: its coordinates, its
#     shell/hole role (decided by exact kernel predicates over the arrangement's
#     nodes — `_ring_is_ccw_exact`), and hole containment (design §3
#     amendment 5 — never planar even-odd on emitted coordinates; see the KNOWN
#     GAP note there, containment is still decided on emitted coordinates).
#
# Ring linkage lives in the phase-2a `OverlayEdge` handles: `next_result_max`
# links maximal rings, `next_result` links minimal rings, and `max_edge_ring` /
# `edge_ring` are integer handles (`0` = null) into a builder's ring vectors.
# Each `NodedEdge`/`OverlayEdge` is a single straight segment between two nodes,
# so a ring's coordinates are exactly the sequence of its edges' node points —
# no intermediate `addCoordinates` bookkeeping is needed.
#
# Everything here is internal to GeometryOps — nothing is exported.

# A maximal edge ring: a cycle of result-area half-edges linked by
# `next_result_max`. Identified by `id` (stored in each member's `max_edge_ring`
# handle) so `_attach_max_edges!` can detect revisits (port of JTS
# `MaximalEdgeRing`).
mutable struct _MaxEdgeRing
    id::Int32
    start_edge::Int32
end

# One minimal result ring (port of JTS `OverlayEdgeRing`): the emitted output
# coordinates (`ring_pts`, closed), the arrangement node ids the ring visits
# (`node_ids`, open — the ring's EXACT identity, from which the shell/hole role
# is decided), its shell/hole role, a bounding box for extent pruning, its
# assigned shell / contained holes (handles, `0` = null), and a lazily-built
# indexed point-in-area locator over its own ring.
#
# `node_ids` and `ring_pts` are not in bijection: two distinct nodes can realize
# to the same output coordinate, which `_ring_add!` drops from `ring_pts` (see
# below). That is exactly why the role is taken from `node_ids`.
mutable struct _OverlayEdgeRing{P}
    id::Int32
    start_edge::Int32
    ring_pts::Vector{Tuple{Float64, Float64}}
    node_ids::Vector{Int32}
    is_hole::Bool
    bbox::NTuple{4, Float64}   # (xmin, xmax, ymin, ymax) of ring_pts
    shell::Int32
    holes::Vector{Int32}
    locator::Any               # Union{Nothing, IndexedPointInAreaLocator}, lazy
end

# The polygon-builder working context (JTS `PolygonBuilder`'s mutable state). Held
# together so the ring-linking and placement functions share the graph edge store,
# the arrangement (for `node_point`), the manifold/exact predicate context, and the
# growing ring collections. Parameterized on `M`/`P`/`E` so `m`/`exact` stay
# concrete and the `Planar`/`Spherical` methods dispatch.
mutable struct _PolyBuilderCtx{M <: Manifold, P, E}
    m::M
    edges::Vector{OverlayEdge{P}}
    arr::NodedArrangement{P}
    exact::E
    max_rings::Vector{_MaxEdgeRing}
    edge_rings::Vector{_OverlayEdgeRing{P}}
    shell_list::Vector{Int32}      # handles into edge_rings
    free_hole_list::Vector{Int32}
    #-- per-node memo of the exact kernel position used by the orientation
    #-- predicate (spherical only; see `_node_kernel_point`). Allocated on
    #-- first use so the planar path pays nothing.
    kernel_cache::Vector{P}
    kernel_ok::Vector{Bool}
end

_PolyBuilderCtx(m::M, edges::Vector{OverlayEdge{P}}, arr::NodedArrangement{P}, exact::E,
        max_rings, edge_rings, shell_list, free_hole_list) where {M <: Manifold, P, E} =
    _PolyBuilderCtx{M, P, E}(m, edges, arr, exact, max_rings, edge_rings, shell_list,
                             free_hole_list, P[], Bool[])

@inline _ctx_point_type(::_PolyBuilderCtx{M, P}) where {M, P} = P

# ## Maximal-ring linking at a node (port of `linkResultAreaMaxRingAtNode`)
#
# Design §3 amendment 4 (the known trap): this is called UNGATED for every
# in-result edge (JTS's `// TODO: skip already-linked` is deliberately
# unfulfilled — gating on already-linked loses degree-2 nodes whose lone
# out-edge was pre-linked as an in-edge). The per-node scan's own early return
# on an already-linked in-edge provides the necessary idempotency.
function _link_result_area_max_ring_at_node!(edges, node_edge::Integer)
    #-- precondition: node_edge is in the result area
    end_out = he_onext(edges, node_edge)
    curr_out = end_out
    #-- state machine: 1 = find an incoming result edge, 2 = link to an outgoing
    state = 1
    curr_result_in = Int32(0)
    while true
        #-- if the found in-edge is already linked, this node is done
        (curr_result_in != 0 && oe_is_result_max_linked(edges, curr_result_in)) && return nothing
        if state == 1
            curr_in = he_sym(edges, curr_out)
            if oe_in_result_area(edges, curr_in)
                curr_result_in = curr_in
                state = 2
            end
        else # state == 2
            if oe_in_result_area(edges, curr_out)
                oe_set_next_result_max!(edges, curr_result_in, curr_out)
                state = 1
            end
        end
        curr_out = he_onext(edges, curr_out)
        curr_out == end_out && break
    end
    state == 2 && throw(_OverlayTopologyError("no outgoing edge found"))
    return nothing
end

# Attach the edges of a maximal ring, tagging each with `mr.id` (port of the
# `MaximalEdgeRing` constructor / `attachEdges`).
function _attach_max_edges!(ctx, mr::_MaxEdgeRing)
    edges = ctx.edges
    edge = mr.start_edge
    while true
        edge == 0 && throw(_OverlayTopologyError("Ring edge is null"))
        edges[edge].max_edge_ring == mr.id &&
            throw(_OverlayTopologyError("Ring edge visited twice at max-ring build"))
        oe_next_result_max(edges, edge) == 0 &&
            throw(_OverlayTopologyError("Ring edge missing at max-ring build"))
        edges[edge].max_edge_ring = mr.id
        edge = oe_next_result_max(edges, edge)
        edge == mr.start_edge && break
    end
    return nothing
end

# ## Minimal-ring split (ports of `buildMinimalRings` / `linkMinimalRings`)
#
# Splits a self-touching maximal ring into OGC-valid minimal rings by relinking
# the max-ring edges in the OPPOSITE (CW) orientation via `next_result`. This is
# exactly the piece the spike prototypes faked (which produced invalid unions on
# many-island geometries).

# Build the minimal rings of a maximal ring, returning their handles.
function _build_minimal_rings!(ctx, mr::_MaxEdgeRing)
    _link_minimal_rings!(ctx, mr)
    min_rings = Int32[]
    edges = ctx.edges
    e = mr.start_edge
    while true
        edges[e].edge_ring == 0 && push!(min_rings, _new_edge_ring!(ctx, e))
        e = oe_next_result_max(edges, e)
        e == mr.start_edge && break
    end
    return min_rings
end

function _link_minimal_rings!(ctx, mr::_MaxEdgeRing)
    edges = ctx.edges
    e = mr.start_edge
    while true
        _link_min_ring_edges_at_node!(ctx, e, mr)
        e = oe_next_result_max(edges, e)
        e == mr.start_edge && break
    end
    return nothing
end

# Port of `linkMinRingEdgesAtNode`: relink this max ring's edges around one node
# into minimal rings (CW orientation, via `next_result`).
function _link_min_ring_edges_at_node!(ctx, node_edge::Integer, mr::_MaxEdgeRing)
    edges = ctx.edges
    end_out = Int32(node_edge)
    curr_max_ring_out = Int32(node_edge)
    curr_out = he_onext(edges, node_edge)
    while true
        _is_already_linked_min(edges, he_sym(edges, curr_out), mr) && return nothing
        if curr_max_ring_out == 0
            curr_max_ring_out = _select_max_out_edge(edges, curr_out, mr)
        else
            curr_max_ring_out = _link_max_in_edge!(edges, curr_out, curr_max_ring_out, mr)
        end
        curr_out = he_onext(edges, curr_out)
        curr_out == end_out && break
    end
    curr_max_ring_out != 0 &&
        throw(_OverlayTopologyError("Unmatched edge found during min-ring linking"))
    return nothing
end

@inline _is_already_linked_min(edges, edge::Integer, mr::_MaxEdgeRing) =
    edges[edge].max_edge_ring == mr.id && oe_is_result_linked(edges, edge)

@inline _select_max_out_edge(edges, curr_out::Integer, mr::_MaxEdgeRing) =
    edges[curr_out].max_edge_ring == mr.id ? Int32(curr_out) : Int32(0)

@inline function _link_max_in_edge!(edges, curr_out::Integer, curr_max_ring_out::Integer,
        mr::_MaxEdgeRing)
    curr_in = he_sym(edges, curr_out)
    edges[curr_in].max_edge_ring != mr.id && return Int32(curr_max_ring_out)
    oe_set_next_result!(edges, curr_in, curr_max_ring_out)
    return Int32(0)
end

# ## Minimal ring construction (port of the `OverlayEdgeRing` constructor)

#=
Port of JTS `CoordinateList.add(coord, allowRepeated = false)`, the accumulator
`OverlayEdge.addCoordinates` fills a ring's coordinates through: a point equal
to the current last one is not appended.

This is not cosmetic here. The arrangement's nodes are symbolic and exact, so two
*distinct* nodes can realize to the same Float64 coordinate at emission (design
§2.6 rounds once, there); the ring then carries a repeated point, which is not a
legal `LinearRing` vertex sequence. JTS never sees the situation because its
noder rounds while noding, so the two nodes are already one.
=#
@inline function _ring_add!(pts::Vector{Tuple{Float64, Float64}}, p::Tuple{Float64, Float64})
    (isempty(pts) || pts[end] != p) && push!(pts, p)
    return nothing
end

#=
Whether a minimal ring collapsed at emission: after repeated-point removal it has
fewer than three distinct vertices, so it bounds no area and is not a legal
`LinearRing` (JTS's `GeometryFactory.createLinearRing` rejects it outright).

Such a ring is the emitted image of an arrangement ring whose two sides are
distinct exactly but round to the same output coordinates — the exact analogue of
JTS's `DIM_COLLAPSE` edge label, arrived at one stage later because this engine
rounds at emission rather than while noding. It is excluded from the area result
for the same reason JTS excludes a collapse: it contributes zero area, and
keeping it produces an invalid result geometry. No tolerance is involved — the
test is exact equality of the emitted coordinates.
=#
@inline _ring_is_collapsed(r::_OverlayEdgeRing) = length(r.ring_pts) < 4

function _new_edge_ring!(ctx, start::Integer)
    P = _ctx_point_type(ctx)
    id = Int32(length(ctx.edge_rings) + 1)
    ring = _OverlayEdgeRing{P}(id, Int32(start), Tuple{Float64, Float64}[], Int32[],
                               false, (0.0, 0.0, 0.0, 0.0), Int32(0), Int32[], nothing)
    push!(ctx.edge_rings, ring)
    _compute_ring!(ctx, ring)
    return id
end

# Port of `computeRingPts` + `computeRing`: walk the minimal ring via
# `next_result`, collecting the arrangement node ids it visits and their emitted
# points; then derive the shell/hole role and the bounding box.
function _compute_ring!(ctx, ring::_OverlayEdgeRing)
    edges = ctx.edges
    pts = Tuple{Float64, Float64}[]
    ids = Int32[]
    origin = he_origin(edges, ring.start_edge)
    push!(ids, Int32(origin))
    _ring_add!(pts, node_point(ctx.arr, origin))
    edge = ring.start_edge
    while true
        edges[edge].edge_ring == ring.id &&
            throw(_OverlayTopologyError("Edge visited twice during ring-building"))
        dest = he_dest(edges, edge)
        edges[edge].edge_ring = ring.id
        ne = oe_next_result(edges, edge)
        ne == 0 && throw(_OverlayTopologyError("Found null edge in ring"))
        edge = ne
        #-- `ids` stays OPEN: the final dest is the start origin, already pushed
        edge == ring.start_edge && (_ring_add!(pts, node_point(ctx.arr, dest)); break)
        push!(ids, Int32(dest))
        _ring_add!(pts, node_point(ctx.arr, dest))
    end
    #-- the last dest is the start origin, so pts is already closed; be defensive
    pts[end] == pts[1] || push!(pts, pts[1])

    ring.ring_pts = pts
    ring.node_ids = ids
    ring.is_hole = _ring_is_ccw_exact(ctx, ids)
    ring.bbox = _ring_bbox(pts)
    return nothing
end

#=
## Shell-vs-hole from exact kernel predicates (design §0)

Result assembly is a decision, so it is decided on the arrangement, not on its
emitted image. `_ring_is_ccw_exact` answers "is this minimal ring wound CCW?"
from the ring's NODE IDS — the symbolic, exact identities — never from
`ring_pts`, which is `node_point` output: rounded once at emission (design
§2.6). The distinction is not academic. A ring whose exact width is a fraction
of an ULP rounds to a self-touching or inverted image, and `Orientation.isCCW`
over that image reports the opposite role; the ring is then filed as a second
shell or as a free hole no shell contains, which is where the
`OverlayTopologyError`s came from.

Both manifolds take the sign of the ring's exact signed area — the planar
shoelace, the spherical geodesic curvature. Neither reads an output coordinate.
=#

#=
Planar: the sign of the doubled signed area over the exact node coordinates,
as a certified Float64 filter with an exact `Rational{BigInt}` fallback — the
same filter/escalate shape every other predicate in the port uses.

The filter runs on `node_point`, which is the CORRECTLY ROUNDED exact node
(design §2.6: certified double-double, else the rounded rational), so each
coordinate carries at most ½ ulp of error. The bound below covers both that
perturbation (`cmax * len`, a coordinate displacement acting on the ring's
`L1` extent) and the accumulated summation error (`n * mag`). Clearing it
certifies the sign of the EXACT area, not merely of the emitted one; failing
it escalates.

The signed area is used rather than JTS's extreme-vertex "cap" test because
its sign is the ring's net winding whether or not the ring is simple — the cap
test presumes simplicity, which holds for arrangement minimal rings but is one
more premise to carry. Translating to the first vertex costs nothing and keeps
the filter well-conditioned on the corpus's ~1e6-magnitude coordinates.
=#
function _ring_is_ccw_exact(ctx::_PolyBuilderCtx{<:Planar}, ids::Vector{Int32})
    length(ids) < 3 && return false
    (area, certain) = _ring_area2_filter(ctx.arr, ids)
    certain && return area > 0
    return _ring_area2_exact(ctx.arr, ids) > 0
end

# Float64 doubled signed area of the emitted ring, with the certificate above.
function _ring_area2_filter(arr, ids::Vector{Int32})
    n = length(ids)
    ox, oy = node_point(arr, ids[1])
    acc = 0.0; mag = 0.0; len = 0.0; cmax = max(abs(ox), abs(oy))
    ax = 0.0; ay = 0.0                       # translated previous vertex
    @inbounds for i in 2:n
        (px, py) = node_point(arr, ids[i])
        cmax = max(cmax, abs(px), abs(py))
        bx = px - ox; by = py - oy
        t1 = ax * by; t2 = ay * bx
        acc += t1 - t2
        mag += abs(t1) + abs(t2)
        len += abs(bx) + abs(by)
        ax = bx; ay = by
    end
    #-- the closing term (last → first) has both translated factors zero
    u = 0.5 * eps(Float64)
    bound = 8 * u * (n * mag + cmax * len)
    return (acc, abs(acc) > bound)
end

# Exact doubled signed area: the plain shoelace over the arrangement's exact
# node coordinates (`_exact_node_point` — the input vertex for a vertex node,
# `_exact_crossing_point` for a crossing). Rational arithmetic on Float64-derived
# values is exact, so the returned sign is the ring's true winding.
function _ring_area2_exact(arr, ids::Vector{Int32})
    n = length(ids)
    acc = zero(Rational{BigInt})
    (ax, ay) = _exact_node_point(arr.nodes.keys[ids[1]])
    (x1, y1) = (ax, ay)
    @inbounds for i in 2:n
        (bx, by) = _exact_node_point(arr.nodes.keys[ids[i]])
        acc += ax * by - ay * bx
        ax = bx; ay = by
    end
    return acc + (ax * y1 - ay * x1)
end

#=
Spherical: the loop's geodesic curvature (`_ring_is_ccw(::Spherical, …)`, our
port of S2 `GetCurvature`) over the exact node directions. The turning-angle
formulation is apex-free — each term sees only three adjacent vertices, and its
sign comes from the exact `Sign` predicate — so it is the formulation this
codebase has already settled on for spherical orientation, and the planar
extreme-vertex cap is deliberately NOT ported to it.

What changes here is only its input: each vertex is the node's own direction —
the ingested unit vector for a vertex node, the exact crossing direction
`±(na × nb)` for a crossing node — instead of the emitted `(lon, lat)` converted
back through trigonometry. Distinct nodes therefore stay distinct, which
`_prune_loop_degeneracies` would otherwise collapse.
=#
_ring_is_ccw_exact(ctx::_PolyBuilderCtx{<:Spherical}, ids::Vector{Int32}) =
    _ring_is_ccw(ctx.m, [_node_kernel_point(ctx, id) for id in ids]; exact = ctx.exact)

# The exact sphere position of node `id` as a unit kernel point, memoized per
# builder (a node is shared by every ring through it, and the exact crossing
# direction is `Rational{BigInt}` work).
function _node_kernel_point(ctx::_PolyBuilderCtx{<:Spherical, P}, id::Integer) where {P}
    if isempty(ctx.kernel_ok)
        n = num_nodes(ctx.arr)
        ctx.kernel_cache = Vector{P}(undef, n)
        ctx.kernel_ok = fill(false, n)
    end
    i = Int(id)
    @inbounds ctx.kernel_ok[i] && return ctx.kernel_cache[i]
    k = ctx.arr.nodes.keys[i]
    p = if k.is_crossing
        d = _sph_crossing_dir(booltype(ctx.exact), k)
        rk_normalize_usp(UnitSphericalPoint(Float64(d[1]), Float64(d[2]), Float64(d[3])))
    else
        k.pt
    end
    @inbounds ctx.kernel_cache[i] = p
    @inbounds ctx.kernel_ok[i] = true
    return p
end

function _ring_bbox(pts::Vector{Tuple{Float64, Float64}})
    xmin = xmax = pts[1][1]
    ymin = ymax = pts[1][2]
    for p in pts
        xmin = min(xmin, p[1]); xmax = max(xmax, p[1])
        ymin = min(ymin, p[2]); ymax = max(ymax, p[2])
    end
    return (xmin, xmax, ymin, ymax)
end

# ## Hole containment (ports of `OverlayEdgeRing.locate` / `contains` / …)
#
# KNOWN GAP, tracked separately from the shell/hole role above: this half of
# result assembly still decides on `ring_pts`, i.e. on emitted coordinates. The
# ray crossing itself is robust (`rk_orient`, never naive even-odd), but its
# inputs are rounded, and so are the bbox/RTree prefilters that reject candidate
# shells. Making it exact needs a point-in-ring predicate over the arrangement's
# node coordinates on both manifolds (`Rational{BigInt}` crossing counts on the
# plane; a rational-direction crossing-parity scan on the sphere), plus outward
# padding of the two float prefilters so they cannot reject a genuinely contained
# hole, plus a way to index rational coordinates (or to demote the existing
# Float64 y-interval index to a candidate-edge filter). Measured over the robust
# corpus, all 405 free-hole placements are decided by a strictly-INTERIOR hole
# vertex — none falls through the all-BOUNDARY path — so nothing currently turns
# on it.

# Lazily builds an indexed point-in-area locator over this ring and locates `p`
# (design §3 amendment 5: robust ray crossing via `rk_orient`, over the ring's
# emitted coordinates — never naive even-odd).
function _ring_locate(ctx, ring::_OverlayEdgeRing, p)
    if ring.locator === nothing
        ring.locator = IndexedPointInAreaLocator(ctx.m, GI.Polygon([ring.ring_pts]);
                                                 exact = ctx.exact)
    end
    return locate(ring.locator, p)
end

# Whether `shell` contains `hole` (port of `contains` + `isPointInOrOut`). On the
# plane a bounding-box reject prefilters before the point tests; on the sphere the
# lon/lat box is unreliable near the poles/antimeridian, so containment is decided
# purely by the point tests (free holes are rare, this path is cold).
#
# Adaptation to the non-self-noding substrate: JTS uses `containsProperly` here,
# because a hole touching its own shell would have been connected into one maximal
# ring (JTS nodes A against itself). This substrate does NOT self-node a single
# input (design §2.2), so such a hole surfaces as a disconnected free hole whose
# bbox touches its shell's; the prefilter must therefore be non-strict
# (`_bbox_contains`), letting the point-in-area test — the real decision — run.
_ring_contains(ctx::_PolyBuilderCtx{<:Planar}, shell::_OverlayEdgeRing, hole::_OverlayEdgeRing) =
    _bbox_contains(shell.bbox, hole.bbox) && _is_point_in_or_out(ctx, shell, hole)
_ring_contains(ctx::_PolyBuilderCtx{<:Spherical}, shell::_OverlayEdgeRing, hole::_OverlayEdgeRing) =
    _is_point_in_or_out(ctx, shell, hole)

function _is_point_in_or_out(ctx, shell::_OverlayEdgeRing, hole::_OverlayEdgeRing)
    for p in hole.ring_pts
        loc = _ring_locate(ctx, shell, p)
        loc == LOC_INTERIOR && return true
        loc == LOC_EXTERIOR && return false
        #-- LOC_BOUNDARY: inconclusive, keep checking
    end
    return false
end

@inline _bbox_contains(outer::NTuple{4, Float64}, inner::NTuple{4, Float64}) =
    outer[1] <= inner[1] && outer[2] >= inner[2] && outer[3] <= inner[3] && outer[4] >= inner[4]

# Port of `findEdgeRingContaining`: the innermost (smallest-envelope) shell in
# `candidates` that contains `hole`, or `0`.
function _find_edge_ring_containing(ctx, hole::_OverlayEdgeRing, candidates)
    min_containing = Int32(0)
    for sh in candidates
        shell = ctx.edge_rings[sh]
        if _ring_contains(ctx, shell, hole)
            if min_containing == 0 ||
               _bbox_contains(ctx.edge_rings[min_containing].bbox, shell.bbox)
                min_containing = Int32(sh)
            end
        end
    end
    return min_containing
end
