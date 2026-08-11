# NOTE: This functionality is experimental and may change at any time.

# # PolygonBuilder — result polygons from the marked result-area edges
#
# Phase 2b of the OverlayNG port (design doc §3). A port of
# `operation/overlayng/PolygonBuilder.java`: link the result-area edges into
# maximal rings, split those into OGC-valid minimal rings (the ring types and
# linking live in maximal_edge_ring.jl), classify shells vs holes, and assign
# each free hole to its containing shell (design §3 amendment 5: containment via
# the indexed point-in-area locators over kernel points, with an `RTree(STR())`
# over shell extents for pruning, mirroring JTS's HPRtree).
#
# Each `OverlayEdge` is a single segment, so a ring's coordinates are the node
# points of its edges — no `RingClipper`/`LineLimiter` (design §3: replaced by
# whole-ring extent pruning + one PIP per pruned ring).
#
# The second half of this file is the non-dissolving companion extraction — face
# enumeration and its opt-in dangle / cut-edge hygiene, ports of JTS's
# `Polygonizer` graph passes — specified in design doc §5.
#
# Everything here is internal to GeometryOps — nothing is exported.

# ## One extraction per graph (the shared precondition of both entry points)
#
# Ring extraction — the op pipeline below and the face enumeration further down —
# CONSUMES the graph: it writes the `next_result` / `next_result_max` /
# `edge_ring` / `max_edge_ring` linkage fields of the half-edges, which start
# null and are never reset. A second extraction over the same `OverlayGraph`
# therefore reads another extraction's linkage and silently returns a wrong
# answer rather than failing (measured: op-then-faces yields 3 of 4 face rings;
# faces-then-op yields 0 of 1 result polygons). Both entry points check first.

# Whether any ring-linkage field of the graph has been written, i.e. whether some
# extraction has already run over it. This is what is actually detected — not
# "which" extraction ran — and it is exact for both orders, because every
# extraction writes `next_result` for at least the edges it rings.
function _graph_is_ring_consumed(g::OverlayGraph)
    for e in g.edges
        (e.next_result != 0 || e.next_result_max != 0 ||
         e.edge_ring != 0 || e.max_edge_ring != 0) && return true
    end
    return false
end

# Whether any half-edge has been removed by the face-walk hygiene pass.
_graph_has_removed_edges(g::OverlayGraph) = any(e -> e.removed, g.edges)

function _assert_graph_extractable(g::OverlayGraph, who::AbstractString)
    _graph_is_ring_consumed(g) && throw(_OverlayTopologyError(
        "$who: this OverlayGraph has already been consumed by a ring extraction " *
        "(its half-edge ring-linkage fields are non-null). Ring extraction writes " *
        "linkage into the shared graph and never resets it, so only ONE extraction " *
        "(`_build_polygons` or `_build_faces`) is valid per graph — a second one " *
        "would silently return a wrong answer. Build a fresh `OverlayGraph` (and " *
        "re-run `_compute_labelling!`) for each extraction."))
    return nothing
end

# Build the result polygons from the graph's result-area edges (port of the
# `PolygonBuilder` constructor + `getPolygons`).
function _build_polygons(m::Manifold, g::OverlayGraph{P, T}, result_area_edges; exact) where {P, T}
    _assert_graph_extractable(g, "_build_polygons")
    #-- the op pipeline walks `onext`/`sym` directly and has no notion of a removed
    #-- edge, so a hygiene-filtered graph would silently ignore the removal
    _graph_has_removed_edges(g) && throw(_OverlayTopologyError(
        "_build_polygons: this OverlayGraph has had half-edges removed by the " *
        "face-walk hygiene pass (`remove_dangles` / `remove_cut_edges`). The op " *
        "pipeline does not honour removal, so its result would silently ignore it. " *
        "Hygiene is a face-enumeration facility — use `_build_faces` on that graph."))
    ctx = _PolyBuilderCtx(m, g.edges, g.arr, exact, _MaxEdgeRing[], _edge_ring_type(T)[],
                          Int32[], Int32[])
    _build_rings!(ctx, result_area_edges)
    return [_ring_to_polygon(ctx, sh) for sh in ctx.shell_list]
end

# Port of `buildRings`.
function _build_rings!(ctx, result_edges)
    #-- design §3 amendment 4: link UNGATED for every in-result edge
    for e in result_edges
        _link_result_area_max_ring_at_node!(ctx.edges, e)
    end
    max_rings = _build_maximal_rings!(ctx, result_edges)
    for mrid in max_rings
        min_rings = _build_minimal_rings!(ctx, ctx.max_rings[mrid])
        _assign_shells_and_holes!(ctx, min_rings)
    end
    _place_free_holes!(ctx)
    return nothing
end

# Port of `buildMaximalRings`: one `_MaxEdgeRing` per unprocessed in-result
# boundary edge. Returns the max-ring handles.
function _build_maximal_rings!(ctx, edges_iter)
    max_rings = Int32[]
    for e in edges_iter
        if oe_in_result_area(ctx.edges, e) && is_boundary_either(oe_label(ctx.edges, e)) &&
           ctx.edges[e].max_edge_ring == 0
            id = Int32(length(ctx.max_rings) + 1)
            mr = _MaxEdgeRing(id, Int32(e))
            push!(ctx.max_rings, mr)
            _attach_max_edges!(ctx, mr)
            push!(max_rings, id)
        end
    end
    return max_rings
end

# Port of `assignShellsAndHoles`: the minimal rings of one maximal ring are
# either a shell + its holes, or a set of (connected) holes whose shell is found
# later (free holes).
function _assign_shells_and_holes!(ctx, min_rings)
    #-- rings that collapsed at emission have no faithful, legal LinearRing image;
    #-- drop them before the shell/hole census, exactly as JTS drops its
    #-- DIM_COLLAPSE edges before ring building (`_ring_is_collapsed`)
    rings = filter(er -> !_ring_is_collapsed(ctx, ctx.edge_rings[er]), min_rings)
    isempty(rings) && return nothing
    shell = _find_single_shell(ctx, rings)
    if shell != 0
        for er in rings
            ctx.edge_rings[er].is_hole && _set_shell!(ctx, er, shell)
        end
        push!(ctx.shell_list, shell)
    else
        #-- No shell among the survivors. Either this maximal ring never had one
        #-- (JTS's free-hole case: a connected group of holes whose shell is a
        #-- different maximal ring), or it had one and the collapse test just
        #-- dropped it. Those two look identical here but are owed different
        #-- treatment when placement fails, so the second is recorded.
        if length(rings) < length(min_rings) &&
           any(er -> !ctx.edge_rings[er].is_hole, min_rings)
            append!(ctx.orphan_hole_list, rings)
        end
        append!(ctx.free_hole_list, rings)
    end
    return nothing
end

# Port of `findSingleShell`: the single non-hole ring, or `0` (all holes).
function _find_single_shell(ctx, min_rings)
    shell = Int32(0)
    shell_count = 0
    for er in min_rings
        if !ctx.edge_rings[er].is_hole
            shell = Int32(er)
            shell_count += 1
        end
    end
    shell_count <= 1 || throw(_OverlayTopologyError("found two shells in EdgeRing list"))
    return shell
end

# Port of `OverlayEdgeRing.setShell` (+ `addHole`).
function _set_shell!(ctx, hole_er::Integer, shell::Integer)
    ctx.edge_rings[hole_er].shell = Int32(shell)
    shell != 0 && push!(ctx.edge_rings[shell].holes, Int32(hole_er))
    return nothing
end

# ## Free-hole placement (port of `placeFreeHoles`)

#=
A free hole that no shell contains is normally a broken arrangement, and the
throw below is a real invariant check worth keeping.

There is exactly one way it can happen legitimately: the hole's own shell was a
sub-grid sliver and `_ring_is_collapsed` dropped it. The hole then has no
containing shell because the thing that contained it is gone, and it is not a
defect — it is the second half of the same drop.

Such holes are re-offered to placement first rather than dropped with their
shell, because the two are not equivalent: a sub-grid shell says nothing about
its holes' own dimensions (a compact hole can sit inside a long thin shell and
have a mean width the shell does not share), and if some SURVIVING shell really
does contain the hole then that is where it belongs. Only when placement genuinely
fails is the hole concluded to be the interior of the dropped sliver and dropped
with it. Measured on the JTS robust corpus this re-offer never succeeds — the
result's shells are disjoint faces, so nothing else contains it — but the check
costs one pass over a list that is almost always empty and it is the difference
between a justified drop and an assumed one.
=#
@inline _drop_unplaceable(ctx, hole_er) =
    hole_er in ctx.orphan_hole_list ? true :
        throw(_OverlayTopologyError("unable to assign free hole to a shell"))

# Planar: prune candidate shells with an `RTree(STR())` over shell extents
# (design §3 amendment 5, the HPRtree analogue).
function _place_free_holes!(ctx::_PolyBuilderCtx{<:Planar})
    isempty(ctx.free_hole_list) && return nothing
    shells = ctx.shell_list
    if isempty(shells)
        all(er -> er in ctx.orphan_hole_list, ctx.free_hole_list) && return nothing
        throw(_OverlayTopologyError("unable to assign free hole to a shell"))
    end
    exts = [_ext_of(ctx.edge_rings[s].bbox) for s in shells]
    index = RTree(STR(), collect(shells); extents = exts)
    for hole_er in ctx.free_hole_list
        ctx.edge_rings[hole_er].shell == 0 || continue
        hole_ext = _ext_of(ctx.edge_rings[hole_er].bbox)
        cand = Int32[]
        SpatialTreeInterface.depth_first_search(Base.Fix1(Extents.intersects, hole_ext), index) do i
            push!(cand, index.data[i])
        end
        shell = _find_edge_ring_containing(ctx, ctx.edge_rings[hole_er], cand)
        shell == 0 && (_drop_unplaceable(ctx, hole_er); continue)
        _set_shell!(ctx, hole_er, shell)
    end
    return nothing
end

# Spherical: the lon/lat extent prune is unreliable near the poles/antimeridian,
# so test every shell (free holes are rare — this path is cold).
function _place_free_holes!(ctx::_PolyBuilderCtx{<:Spherical})
    isempty(ctx.free_hole_list) && return nothing
    for hole_er in ctx.free_hole_list
        ctx.edge_rings[hole_er].shell == 0 || continue
        shell = _find_edge_ring_containing(ctx, ctx.edge_rings[hole_er], ctx.shell_list)
        shell == 0 && (_drop_unplaceable(ctx, hole_er); continue)
        _set_shell!(ctx, hole_er, shell)
    end
    return nothing
end

# The shell-extent key of the `RTree(STR())` prune above, and nothing else: it is
# built here, consumed by `Extents.intersects` two lines up, and never reaches a
# result geometry or a caller. Planar-only, hence `NTuple{4}` only — the
# spherical `_place_free_holes!` tests every shell instead of pruning, so an xyz
# box never needs an `Extent` spelling.
@inline _ext_of(bbox::NTuple{4, Float64}) =
    Extents.Extent(X = (bbox[1], bbox[2]), Y = (bbox[3], bbox[4]))

# ## Polygon assembly (port of `OverlayEdgeRing.toPolygon`)

# Emit the polygon of one shell ring and its assigned holes. Ring windings are
# left as the graph produced them (JTS does the same); `GO.area` on either
# manifold is orientation-independent, and validity does not depend on winding.
function _ring_to_polygon(ctx::_PolyBuilderCtx{M, P, E, T}, shell_handle::Integer) where {M, P, E, T}
    sh = ctx.edge_rings[shell_handle]
    rings = Vector{Vector{T}}()
    push!(rings, sh.ring_pts)
    for h in sh.holes
        push!(rings, ctx.edge_rings[h].ring_pts)
    end
    return GI.Polygon(rings)
end

# ## Face enumeration — all boundary cycles of the arrangement (Polygonizer-style)
#
# The op pipeline above extracts only the rings the op's result predicate
# selects, after dissolving interior boundaries (`unmarkDuplicateEdges`). Some
# consumers — antimeridian splitting, polygon-cut-by-line, polygonize — instead
# need the arrangement's FACES: every boundary cycle of the noded linework, each
# tracing the face on its RIGHT via the half-edge face traversal (successor =
# onext ∘ sym), with per-input face locations read off the shared labels. A
# dangling edge (a line dead-end) is traversed twice by its face's cycle, out
# and back — callers that need dangle-free rings opt into the hygiene pass
# below, callers like the antimeridian pole seam rely on exactly this doubling.
#
# CYCLES ARE NOT FACES. This enumerates one ring per boundary *cycle*; a face
# with cavities contributes several (its outer cycle plus one per cavity), and
# the unbounded outer face contributes one cycle per connected component. Cycle
# count is `E − V + 2C`, face count only `E − V + C + 1`. The grouping back into
# faces is re-derived geometrically by `_place_free_holes!` (which carries the
# KNOWN GAP note in maximal_edge_ring.jl — containment is decided on rounded
# emitted coordinates). Nothing here may assume a cycle↔face bijection.
#
# Winding: a bounded cycle comes out CW (`is_hole == false`) with its face on the
# right; an outer cycle comes out CCW (`is_hole == true`).
#
# Reuses the `_OverlayEdgeRing` pipeline unchanged: the face link fills the
# same `next_result` field the op pipeline links, so `_compute_ring!` (ring
# points, kernel points, shell/hole orientation, bbox) and the hole-placement
# machinery run identically. Like the op pipeline, this consumes the graph's
# ring-linkage fields — run one extraction per `OverlayGraph`
# (`_assert_graph_extractable` above enforces it).
#
# This path is areal-only: it never marks `in_result_area` / `in_result_line`, so
# `_build_lines` / `_build_points` are unreachable from it.

# The face-cycle successor of a half-edge: `onext ∘ sym`, skipping edges removed
# by the hygiene pass. This is the single primitive the whole face layer turns
# on — `_build_faces` writes it into `next_result`, and `_remove_cut_edges!`
# walks it to label cycles. The skip is JTS `PolygonizeGraph.computeNextCWEdges`
# passing over `marked` edges, expressed on the intact `o_next` star instead of
# by rewiring it. Termination: `sym(i)` is live whenever `i` is (removal is
# symmetric), so the star walk stops there at worst.
@inline function _face_successor(edges, i::Integer)
    e = he_onext(edges, he_sym(edges, i))
    while oe_is_removed(edges, e)
        e = he_onext(edges, e)
    end
    return e
end

# Link every live half-edge to its face-cycle successor and build one
# `_OverlayEdgeRing` per face cycle. Requires a labelled graph
# (`_compute_labelling!`). Returns the builder context; ring handles are
# `1:length(ctx.edge_rings)`.
#
# `remove_dangles` / `remove_cut_edges` are the opt-in hygiene pass documented
# below. Both default to `false`: the default enumeration is the raw cycle
# structure of the arrangement, dangle doubling included, which
# `antimeridian_split` structurally depends on (its pole pair IS the doubling at
# the meridian arc's degree-1 pole endpoint).
function _build_faces(m::Manifold, g::OverlayGraph{P, T}; exact,
        remove_dangles::Bool = false, remove_cut_edges::Bool = false) where {P, T}
    _assert_graph_extractable(g, "_build_faces")
    remove_dangles && _remove_dangles!(g)
    remove_cut_edges && _remove_cut_edges!(g)
    edges = g.edges
    for i in eachindex(edges)
        oe_is_removed(edges, i) && continue
        oe_set_next_result!(edges, i, _face_successor(edges, i))
    end
    ctx = _PolyBuilderCtx(m, edges, g.arr, exact, _MaxEdgeRing[],
                          _edge_ring_type(T)[], Int32[], Int32[])
    for i in eachindex(edges)
        (oe_is_removed(edges, i) || edges[i].edge_ring != 0) && continue
        _new_edge_ring!(ctx, i)
    end
    return ctx
end

# ## Face-walk hygiene — dangle and cut-edge removal (OPT-IN, default OFF)
#
# Ports of JTS `PolygonizeGraph.deleteDangles` / `deleteCutEdges`. Both are
# GRAPH-LEVEL passes: they mark half-edges removed *before* the face walk runs,
# so the walk produces clean rings directly. Nothing post-filters an emitted
# point list — a point list cannot express the difference between a dangle spike
# and a legitimate revisit of a node, and it cannot express a cut edge at all.
#
# The two JTS cases are genuinely distinct:
#
#   * DANGLES — an edge incident on a degree-1 node. The face walk runs out to
#     the dead end and back, so the spike lands in the ring. Removal is
#     iterative: deleting a dangle can expose a new degree-1 node behind it, so
#     whole dangle trees peel off (JTS uses an explicit node stack; so do we).
#     Purely structural — live degree only, no ring labels needed.
#
#   * CUT EDGES / bridges — an edge traversed twice *in the same cycle, in
#     opposite directions*, i.e. whose two directed halves carry the same face
#     cycle. JTS detects exactly this (`de.getLabel() == sym.getLabel()` after
#     `findLabeledEdgeRings`); we compare face-cycle ids computed with
#     `_face_successor`. A dangle is the degree-1 special case of a bridge, so
#     this pass subsumes dangle removal — including whole chains, since every
#     edge of a bridge chain lies on one cycle and one pass catches them all.
#
# Both are structural/exact: an edge is or is not incident on a degree-1 node, is
# or is not traversed twice by one cycle. No distance, tolerance or precision
# model is involved anywhere (design §0).
#
# Removal is non-destructive to the arrangement AND to the ring linkage: the
# `o_next` stars are never rewired and no `next_result*` / `*edge_ring` field is
# touched, so a hygiene-filtered graph is still an unconsumed graph. It is not
# reversible, though, and the op pipeline cannot honour it — `_build_polygons`
# rejects a filtered graph rather than quietly ignoring the removal.
#
# Both functions return the removed edges as one representative half-edge index
# per removed edge (the removed linework, JTS's returned `LineString`s).

# Port of `deleteDangles`: iteratively remove every edge incident on a degree-1
# node until no live degree-1 node remains.
function _remove_dangles!(g::OverlayGraph)
    _assert_graph_extractable(g, "_remove_dangles!")
    edges = g.edges
    removed = Int32[]
    stack = Int32[]
    for nid in eachindex(g.node_edges)
        e0 = @inbounds g.node_edges[nid]
        e0 == 0 && continue
        oe_live_degree(edges, e0) == 1 && push!(stack, Int32(nid))
    end
    while !isempty(stack)
        nid = pop!(stack)
        e0 = @inbounds g.node_edges[nid]
        e0 == 0 && continue
        #-- JTS `deleteAllEdges(node)`: the node was stacked at live degree 1 and
        #-- degree only ever decreases, so this removes that one edge (or none, if
        #-- the node was stacked twice and is already bare).
        e = e0
        while true
            if !oe_is_removed(edges, e)
                oe_remove_both!(edges, e)
                push!(removed, Int32(e))
                to = he_dest(edges, e)
                te = @inbounds g.node_edges[to]
                #-- the far node may have become a dangle in turn
                te != 0 && oe_live_degree(edges, te) == 1 && push!(stack, Int32(to))
            end
            e = he_onext(edges, e)
            e == e0 && break
        end
    end
    return removed
end

# Port of `deleteCutEdges`: label the live face cycles, then remove every edge
# whose two directed halves carry the same cycle id.
function _remove_cut_edges!(g::OverlayGraph)
    _assert_graph_extractable(g, "_remove_cut_edges!")
    edges = g.edges
    #-- scratch cycle ids, NOT the graph's `edge_ring` handles: this pass must
    #-- leave the graph unconsumed for the face walk that follows it.
    cycle = zeros(Int32, length(edges))
    ncycles = Int32(0)
    for i in eachindex(edges)
        (oe_is_removed(edges, i) || cycle[i] != 0) && continue
        ncycles += Int32(1)
        e = Int32(i)
        while true
            @inbounds cycle[e] = ncycles
            e = _face_successor(edges, e)
            e == i && break
        end
    end
    removed = Int32[]
    for i in eachindex(edges)
        oe_is_removed(edges, i) && continue
        s = he_sym(edges, i)
        i < s || continue                     # visit each edge pair once
        if (@inbounds cycle[i]) == (@inbounds cycle[s])
            oe_remove_both!(edges, i)
            push!(removed, Int32(i))
        end
    end
    return removed
end

# The location of ring `er`'s face — the face on the ring's RIGHT — for input
# `gi`. Prefers a boundary edge of `gi` (side locations are authoritative);
# falls back to the first edge's boundary-or-line location.
function _face_ring_location(ctx, er::Integer, gi::Integer)
    edges = ctx.edges
    start = ctx.edge_rings[er].start_edge
    e = start
    while true
        is_boundary(oe_label(edges, e), gi) &&
            return oe_get_location(edges, e, gi, POS_RIGHT)
        e = oe_next_result(edges, e)
        e == start && break
    end
    return oe_get_location_boundary_or_line(edges, start, gi, POS_RIGHT)
end

# Select the face cycles `keep(loc_a, loc_b)` accepts into the builder context's
# shell / free-hole lists and assign each free hole to its shell. `keep` sees the
# raw per-input face locations. Kept clockwise cycles are face shells; kept
# counter-clockwise cycles are cavities, placed into their shells by the same
# containment machinery the op pipeline uses.
#
# This is THE face-selection path — every face consumer runs through it, so that
# the collapse filter, the shell/hole split and the hole placement are stated
# once. Collapsed rings are dropped here for the same reason
# `_assign_shells_and_holes!` drops them from the op result: no Float64 image of
# them is a legal, faithful `LinearRing` — either they have fewer than three
# distinct emitted vertices, or their exact width is finer than the output
# format resolves where they sit (`_ring_is_collapsed`).
function _select_faces!(ctx, keep::F) where {F}
    for er in 1:length(ctx.edge_rings)
        ring = ctx.edge_rings[er]
        _ring_is_collapsed(ctx, ring) && continue
        keep(_face_ring_location(ctx, er, 0), _face_ring_location(ctx, er, 1)) || continue
        ring.is_hole ? push!(ctx.free_hole_list, Int32(er)) :
                       push!(ctx.shell_list, Int32(er))
    end
    _place_free_holes!(ctx)
    return ctx
end

# Build the polygons of the faces `keep` selects. `remove_dangles` /
# `remove_cut_edges` forward to the opt-in hygiene pass in `_build_faces`.
function _build_face_polygons(m::Manifold, g::OverlayGraph, keep::F; exact,
        remove_dangles::Bool = false, remove_cut_edges::Bool = false) where {F}
    ctx = _build_faces(m, g; exact, remove_dangles, remove_cut_edges)
    _select_faces!(ctx, keep)
    return [_ring_to_polygon(ctx, sh) for sh in ctx.shell_list]
end
