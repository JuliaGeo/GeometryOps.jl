# # Result ring linking — maximal rings, minimal-ring split, edge rings
#
# Phase 2b of the OverlayNG port (design doc §3). Ports two JTS files kept
# together because they form one pipeline over the marked result-area edges:
#   - `MaximalEdgeRing.java`  → maximal-ring linking + the minimal-ring split
#     at self-touching nodes (`_MaxEdgeRing` + `_link_result_area_max_ring_at_node!`).
#   - `OverlayEdgeRing.java`  → one minimal result ring: its coordinates, its
#     shell/hole role (via `_ring_is_ccw`), and hole containment (design §3
#     amendment 5 — never planar even-odd on emitted coordinates).
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
# coordinates (`ring_pts`, closed), the manifold kernel points backing the
# orientation and PIP predicates (`kernel_pts`; identical to `ring_pts` on the
# plane), its shell/hole role, a bounding box for extent pruning, its assigned
# shell / contained holes (handles, `0` = null), and a lazily-built indexed
# point-in-area locator over its own ring.
mutable struct _OverlayEdgeRing{P}
    id::Int32
    start_edge::Int32
    ring_pts::Vector{Tuple{Float64, Float64}}
    kernel_pts::Vector{P}
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
end

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

function _new_edge_ring!(ctx, start::Integer)
    P = _ctx_point_type(ctx)
    id = Int32(length(ctx.edge_rings) + 1)
    ring = _OverlayEdgeRing{P}(id, Int32(start), Tuple{Float64, Float64}[], P[],
                               false, (0.0, 0.0, 0.0, 0.0), Int32(0), Int32[], nothing)
    push!(ctx.edge_rings, ring)
    _compute_ring!(ctx, ring)
    return id
end

# Port of `computeRingPts` + `computeRing`: walk the minimal ring via
# `next_result`, collecting node points; then derive the shell/hole role and the
# bounding box.
function _compute_ring!(ctx, ring::_OverlayEdgeRing)
    edges = ctx.edges
    pts = Tuple{Float64, Float64}[]
    push!(pts, node_point(ctx.arr, he_origin(edges, ring.start_edge)))
    edge = ring.start_edge
    while true
        edges[edge].edge_ring == ring.id &&
            throw(_OverlayTopologyError("Edge visited twice during ring-building"))
        push!(pts, node_point(ctx.arr, he_dest(edges, edge)))
        edges[edge].edge_ring = ring.id
        ne = oe_next_result(edges, edge)
        ne == 0 && throw(_OverlayTopologyError("Found null edge in ring"))
        edge = ne
        edge == ring.start_edge && break
    end
    #-- the last dest is the start origin, so pts is already closed; be defensive
    pts[end] == pts[1] || push!(pts, pts[1])

    ring.ring_pts = pts
    ring.kernel_pts = _ring_kernel_pts(ctx.m, pts)
    ring.is_hole = _ring_is_ccw(ctx.m, ring.kernel_pts; exact = ctx.exact)
    ring.bbox = _ring_bbox(pts)
    return nothing
end

# Emitted output coordinates ARE the planar kernel points; the sphere converts
# each realized (lon, lat) back to a unit vector for the kernel predicates.
_ring_kernel_pts(::Planar, pts::Vector{Tuple{Float64, Float64}}) = pts
_ring_kernel_pts(m::Spherical, pts::Vector{Tuple{Float64, Float64}}) =
    [_to_kernel_point(m, p) for p in pts]

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
