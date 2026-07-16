# # Overlay graph — edge merge, half-edge pairs, and the topology graph
#
# Phase 2a of the OverlayNG port (design doc §3). Consumes the phase-1
# `NodedArrangement` and the per-string `EdgeSourceInfo` table and produces the
# `OverlayGraph`: for every distinct noded edge a symmetric pair of `OverlayEdge`
# half-edges sharing one `OverlayLabel`, with each node's star ordered CCW.
#
# Three JTS files are ported here, kept together because they form one pipeline:
#   - `Edge.java`        → `MergeEdge` + `merge!` (depth-delta summing, collapse
#                          detection, hole-role merge, label creation).
#   - `EdgeMerger.java`  → `_merge_noded_edges` (coincident edges grouped by the
#                          unordered node-id pair — design §3 amendment 2 —
#                          never by coordinates).
#   - `OverlayEdge.java` / `OverlayGraph.java` → `OverlayEdge` + `OverlayGraph`.
#
# Everything here is internal to GeometryOps — nothing is exported.

# ==========================================================================
# Edge (JTS `Edge`): the merge intermediate + label creation
# ==========================================================================
#
# Renamed `MergeEdge` to avoid collision with `OverlayEdge` and RelateNG's
# `RelateEdge` in the single GeometryOps module. It carries no coordinates: its
# geometry is the node pair `(node_lo, node_hi)` (canonical direction = the base
# contributor's parent traversal order) plus the base contributor's parent
# `(string_idx, seg_idx)`, from which the half-edge direction points (the parent
# segment's original endpoints) are read at graph-build time. Topology info is the
# JTS `Edge` a/b `dim`/`depth_delta`/`is_hole` fields.

mutable struct MergeEdge
    node_lo     :: Int32   # canonical origin node id (base contributor's sub-edge start)
    node_hi     :: Int32   # canonical dest node id (base contributor's sub-edge end)
    string_idx  :: Int32   # base contributor's parent string (for direction points)
    seg_idx     :: Int32   # base contributor's parent segment

    a_dim         :: Int8
    a_depth_delta :: Int32  # summed across merges — widened from the Int8 source delta
    a_is_hole     :: Bool

    b_dim         :: Int8
    b_depth_delta :: Int32
    b_is_hole     :: Bool
end

# Build the base `MergeEdge` for a noded edge, carrying its single source's info
# on the matching input index (port of JTS `Edge(pts, info)` + `copyInfo`).
function _merge_edge(ne::NodedEdge, src::EdgeSourceInfo)
    if src.index == 0
        return MergeEdge(ne.node_lo, ne.node_hi, ne.string_idx, ne.seg_idx,
                         src.dim, Int32(src.depth_delta), src.is_hole,
                         DIM_NOT_PART, Int32(0), false)
    else
        return MergeEdge(ne.node_lo, ne.node_hi, ne.string_idx, ne.seg_idx,
                         DIM_NOT_PART, Int32(0), false,
                         src.dim, Int32(src.depth_delta), src.is_hole)
    end
end

# Whether this edge is part of a shell for input `gi` (a non-hole boundary).
_me_is_shell(e::MergeEdge, gi::Integer) = gi == 0 ? (e.a_dim == DIM_A && !e.a_is_hole) :
                                                    (e.b_dim == DIM_A && !e.b_is_hole)

# The merged hole role for input `gi`: a shell if any contributor is a shell
# (port of JTS `isHoleMerged`; `is_hole` is stored, hence the flip).
_me_is_hole_merged(gi::Integer, e1::MergeEdge, e2::MergeEdge) =
    !(_me_is_shell(e1, gi) || _me_is_shell(e2, gi))

# Merge a coincident contributor `inc` into `base` (port of JTS `Edge.merge`).
# Hole status is updated first (it reads the pre-update dims); dims take the max;
# depth deltas sum with a direction flip. The flip is exact from node ids: the
# unordered pair is unique, so `inc` runs the same direction as `base` iff their
# `node_lo` agree (design §3 amendment 2) — no coordinate comparison.
function _merge!(base::MergeEdge, inc::MergeEdge)
    base.a_is_hole = _me_is_hole_merged(0, base, inc)
    base.b_is_hole = _me_is_hole_merged(1, base, inc)
    inc.a_dim > base.a_dim && (base.a_dim = inc.a_dim)
    inc.b_dim > base.b_dim && (base.b_dim = inc.b_dim)
    flip = inc.node_lo == base.node_lo ? Int32(1) : Int32(-1)
    base.a_depth_delta += flip * inc.a_depth_delta
    base.b_depth_delta += flip * inc.b_depth_delta
    return base
end

# ## Label creation from a merged edge (port of JTS `Edge.createLabel` + statics)

@inline _del_sign(d::Integer) = d > 0 ? 1 : (d < 0 ? -1 : 0)

# Positive delta ⇒ Left = EXTERIOR, Right = INTERIOR (JTS `locationLeft/Right`).
@inline function _location_left(depth_delta::Integer)
    s = _del_sign(depth_delta)
    return s == 1 ? LOC_EXTERIOR : (s == -1 ? LOC_INTERIOR : LOC_NONE)
end
@inline function _location_right(depth_delta::Integer)
    s = _del_sign(depth_delta)
    return s == 1 ? LOC_INTERIOR : (s == -1 ? LOC_EXTERIOR : LOC_NONE)
end

# The effective edge-state dimension of a source (port of JTS `Edge.labelDim`):
# an area edge with zero summed depth delta is a collapse.
@inline function _label_dim(dim::Int8, depth_delta::Integer)
    dim == DIM_NOT_PART && return DIM_NOT_PART
    dim == DIM_L && return DIM_L
    #-- dim == DIM_A (area)
    return depth_delta == 0 ? DIM_COLLAPSE : DIM_A
end

# Initialize one source's slot of a label (port of JTS `Edge.initLabel`).
function _init_label!(l::OverlayLabel, gi::Integer, dim::Int8, depth_delta::Integer, is_hole::Bool)
    dl = _label_dim(dim, depth_delta)
    if dl == DIM_NOT_PART
        init_not_part!(l, gi)
    elseif dl == DIM_A            # boundary
        init_boundary!(l, gi, _location_left(depth_delta), _location_right(depth_delta), is_hole)
    elseif dl == DIM_COLLAPSE
        init_collapse!(l, gi, is_hole)
    elseif dl == DIM_L            # line
        init_line!(l, gi)
    end
    return l
end

# The shared `OverlayLabel` for a merged edge (port of JTS `Edge.createLabel`).
function _create_label(e::MergeEdge)
    l = OverlayLabel()
    _init_label!(l, 0, e.a_dim, e.a_depth_delta, e.a_is_hole)
    _init_label!(l, 1, e.b_dim, e.b_depth_delta, e.b_is_hole)
    return l
end

# ## EdgeMerger (port of JTS `EdgeMerger.merge`, keyed by unordered node pair)

@inline _edge_key(ne::NodedEdge) =
    ne.node_lo < ne.node_hi ? (ne.node_lo, ne.node_hi) : (ne.node_hi, ne.node_lo)

# Merge all noded edges of the arrangement. Coincident edges (same unordered node
# pair; the segment/minor-arc between two nodes is unique, antipodal edges
# excluded at ingest) collapse to one `MergeEdge`, the first seen setting the
# canonical direction and later ones merged into it.
function _merge_noded_edges(arr::NodedArrangement, sources::Vector{EdgeSourceInfo})
    edgemap = Dict{Tuple{Int32, Int32}, Int}()
    merged = MergeEdge[]
    for ne in arr.edges
        src = sources[ne.string_idx]
        key = _edge_key(ne)
        idx = get(edgemap, key, 0)
        if idx == 0
            push!(merged, _merge_edge(ne, src))
            edgemap[key] = length(merged)
        else
            _merge!(merged[idx], _merge_edge(ne, src))
        end
    end
    return merged
end

# ==========================================================================
# OverlayEdge (JTS `OverlayEdge`): the graph half-edge
# ==========================================================================
#
# A directed half-edge in the graph. Carries the winged-edge fields required by
# half_edge.jl (`origin`, `sym`, `o_next`, `dir_pt`), the shared `label`, the
# `is_forward` direction (which reinterprets the label's Left/Right), the result-
# marking flags, and the ring-linkage pointers phase 2b fills in. Ring pointers
# are integer half-edge indices (`0` = null); the `edge_ring`/`max_edge_ring`
# fields are integer handles into phase 2b's ring collections (`0` = null).

mutable struct OverlayEdge{P}
    origin :: Int32          # origin node id
    sym    :: Int32          # index of the symmetric half-edge
    o_next :: Int32          # next half-edge CCW around the origin (same origin)
    dir_pt :: P              # direction point: parent segment's far endpoint (original vertex)

    is_forward :: Bool       # direction relative to the shared label
    label      :: OverlayLabel

    in_result_area :: Bool
    in_result_line :: Bool
    visited        :: Bool

    next_result     :: Int32 # next half-edge in the result ring (0 = null)
    next_result_max :: Int32 # next half-edge in the result maximal ring (0 = null)
    edge_ring       :: Int32 # phase 2b OverlayEdgeRing handle (0 = null)
    max_edge_ring   :: Int32 # phase 2b MaximalEdgeRing handle (0 = null)
end

_overlay_edge(origin::Int32, dir_pt::P, is_forward::Bool, label::OverlayLabel) where {P} =
    OverlayEdge{P}(origin, Int32(0), Int32(0), dir_pt, is_forward, label,
                   false, false, false, Int32(0), Int32(0), Int32(0), Int32(0))

# ## OverlayEdge accessors (ports of the JTS `OverlayEdge` surface, index-based)

@inline oe_is_forward(edges, i::Integer) = @inbounds edges[i].is_forward
@inline oe_label(edges, i::Integer) = @inbounds edges[i].label

# Location of `position` for input `index`, resolved for this edge's orientation
# (port of JTS `OverlayEdge.getLocation`).
@inline oe_get_location(edges, i::Integer, index::Integer, position::Integer) =
    get_location((@inbounds edges[i].label), index, position, (@inbounds edges[i].is_forward))

@inline oe_get_location_boundary_or_line(edges, i::Integer, index::Integer, position::Integer) =
    get_location_boundary_or_line((@inbounds edges[i].label), index, position, (@inbounds edges[i].is_forward))

# Result-area marking (ports of the `markInResultArea*` / `unmark*` family).
@inline oe_in_result_area(edges, i::Integer) = @inbounds edges[i].in_result_area
@inline oe_in_result_area_both(edges, i::Integer) =
    (@inbounds edges[i].in_result_area) && (@inbounds edges[he_sym(edges, i)].in_result_area)
@inline function oe_mark_in_result_area!(edges, i::Integer)
    @inbounds edges[i].in_result_area = true; return nothing
end
@inline function oe_mark_in_result_area_both!(edges, i::Integer)
    @inbounds edges[i].in_result_area = true
    @inbounds edges[he_sym(edges, i)].in_result_area = true
    return nothing
end
@inline function oe_unmark_from_result_area_both!(edges, i::Integer)
    @inbounds edges[i].in_result_area = false
    @inbounds edges[he_sym(edges, i)].in_result_area = false
    return nothing
end

# Result-line marking (marks both members of the pair, per JTS).
@inline oe_in_result_line(edges, i::Integer) = @inbounds edges[i].in_result_line
@inline function oe_mark_in_result_line!(edges, i::Integer)
    @inbounds edges[i].in_result_line = true
    @inbounds edges[he_sym(edges, i)].in_result_line = true
    return nothing
end
@inline oe_in_result(edges, i::Integer) =
    (@inbounds edges[i].in_result_area) || (@inbounds edges[i].in_result_line)
@inline oe_in_result_either(edges, i::Integer) =
    oe_in_result(edges, i) || oe_in_result(edges, he_sym(edges, i))

# Visited flag (marks both members, per JTS `markVisitedBoth`).
@inline oe_is_visited(edges, i::Integer) = @inbounds edges[i].visited
@inline function oe_mark_visited_both!(edges, i::Integer)
    @inbounds edges[i].visited = true
    @inbounds edges[he_sym(edges, i)].visited = true
    return nothing
end

# Result-ring linkage getters/setters (`0` = null).
@inline oe_next_result(edges, i::Integer) = @inbounds edges[i].next_result
@inline function oe_set_next_result!(edges, i::Integer, e::Integer)
    @inbounds edges[i].next_result = Int32(e); return nothing
end
@inline oe_is_result_linked(edges, i::Integer) = (@inbounds edges[i].next_result) != 0
@inline oe_next_result_max(edges, i::Integer) = @inbounds edges[i].next_result_max
@inline function oe_set_next_result_max!(edges, i::Integer, e::Integer)
    @inbounds edges[i].next_result_max = Int32(e); return nothing
end
@inline oe_is_result_max_linked(edges, i::Integer) = (@inbounds edges[i].next_result_max) != 0

# ==========================================================================
# OverlayGraph (JTS `OverlayGraph`): build from arrangement + sources
# ==========================================================================

"""
    OverlayGraph{P}

The topology graph of an overlay operation (port of JTS `OverlayGraph`). Holds
the arrangement it was built from, the vector of half-edges (both orientations of
every merged edge), and one representative outgoing half-edge index per node id
(`node_edges[nid]`, `0` if the node has no incident edges). `P` is the manifold
kernel point type, so the graph is type-erased over input geometry types.
"""
struct OverlayGraph{P}
    arr        :: NodedArrangement{P}
    edges      :: Vector{OverlayEdge{P}}
    node_edges :: Vector{Int32}
end

"""
    OverlayGraph(m, arr::NodedArrangement, sources) -> OverlayGraph

Build the overlay graph from a noded arrangement and its `EdgeSourceInfo` table.
Coincident noded edges are merged (JTS `Edge.merge` semantics), each merged edge
becomes a symmetric `OverlayEdge` pair sharing one label, and every node's star
is ordered CCW about its symbolic apex via the exact kernel comparator.
"""
function OverlayGraph(m::Manifold, arr::NodedArrangement{P}, sources::Vector{EdgeSourceInfo};
        exact = True()) where {P}
    merged = _merge_noded_edges(arr, sources)
    nnodes = num_nodes(arr)
    edges = Vector{OverlayEdge{P}}()
    sizehint!(edges, 2 * length(merged))
    stars = [Int32[] for _ in 1:nnodes]
    for me in merged
        label = _create_label(me)
        ss = arr.segstrings[me.string_idx]
        #-- direction points are the parent segment's original endpoints, so the
        #-- forward half-edge (origin node_lo) heads toward the node_hi side and
        #-- the reverse (origin node_hi) toward the node_lo side (design §3.1).
        fwd_dir = ss.pts[me.seg_idx + 1]
        bwd_dir = ss.pts[me.seg_idx]
        push!(edges, _overlay_edge(me.node_lo, fwd_dir, true, label))
        i_fwd = Int32(length(edges))
        push!(edges, _overlay_edge(me.node_hi, bwd_dir, false, label))
        i_rev = Int32(length(edges))
        he_link!(edges, i_fwd, i_rev)
        push!(stars[me.node_lo], i_fwd)
        push!(stars[me.node_hi], i_rev)
    end
    node_edges = zeros(Int32, nnodes)
    _order_all_stars!(m, edges, arr.nodes.keys, stars, node_edges; exact)
    return OverlayGraph{P}(arr, edges, node_edges)
end

# Convenience: build the sources and the graph directly from an arrangement.
function OverlayGraph(m::Manifold, arr::NodedArrangement; exact = True())
    return OverlayGraph(m, arr, _edge_source_infos(m, arr; exact); exact)
end

# Order every node's star once (function barrier: the abstract `m` dispatches into
# the exact comparator here, off the type-stable build loop).
function _order_all_stars!(m::Manifold, edges, keys, stars, node_edges; exact)
    for nid in eachindex(stars)
        star = stars[nid]
        isempty(star) && continue
        he_order_star!(m, edges, keys, star; exact)
        @inbounds node_edges[nid] = star[1]
    end
    return nothing
end

# ## Graph queries (ports of the `OverlayGraph` accessor surface)

# All half-edges (both orientations), matching JTS `getEdges()`.
graph_edges(g::OverlayGraph) = g.edges

# A representative outgoing half-edge index for node `nid`, or `0` (JTS
# `getNodeEdge`).
graph_node_edge(g::OverlayGraph, nid::Integer) = @inbounds g.node_edges[nid]

# The representative node edges (one outgoing half-edge per non-empty node), for
# node-star iteration (JTS `getNodeEdges()`).
graph_node_edges(g::OverlayGraph) = Int32[e for e in g.node_edges if e != 0]

# The half-edge indices marked as being in the result area (JTS
# `getResultAreaEdges`).
graph_result_area_edges(g::OverlayGraph) =
    Int32[i for i in eachindex(g.edges) if g.edges[i].in_result_area]
