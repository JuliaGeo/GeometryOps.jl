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
# Everything here is internal to GeometryOps — nothing is exported.

# Build the result polygons from the graph's result-area edges (port of the
# `PolygonBuilder` constructor + `getPolygons`).
function _build_polygons(m::Manifold, g::OverlayGraph{P}, result_area_edges; exact) where {P}
    ctx = _PolyBuilderCtx(m, g.edges, g.arr, exact, _MaxEdgeRing[], _OverlayEdgeRing{P}[],
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
    shell = _find_single_shell(ctx, min_rings)
    if shell != 0
        for er in min_rings
            ctx.edge_rings[er].is_hole && _set_shell!(ctx, er, shell)
        end
        push!(ctx.shell_list, shell)
    else
        append!(ctx.free_hole_list, min_rings)
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

# Planar: prune candidate shells with an `RTree(STR())` over shell extents
# (design §3 amendment 5, the HPRtree analogue).
function _place_free_holes!(ctx::_PolyBuilderCtx{<:Planar})
    isempty(ctx.free_hole_list) && return nothing
    shells = ctx.shell_list
    if isempty(shells)
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
        shell == 0 && throw(_OverlayTopologyError("unable to assign free hole to a shell"))
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
        shell == 0 && throw(_OverlayTopologyError("unable to assign free hole to a shell"))
        _set_shell!(ctx, hole_er, shell)
    end
    return nothing
end

@inline _ext_of(bbox::NTuple{4, Float64}) =
    Extents.Extent(X = (bbox[1], bbox[2]), Y = (bbox[3], bbox[4]))

# ## Polygon assembly (port of `OverlayEdgeRing.toPolygon`)

# Emit the polygon of one shell ring and its assigned holes. Ring windings are
# left as the graph produced them (JTS does the same); `GO.area` on either
# manifold is orientation-independent, and validity does not depend on winding.
function _ring_to_polygon(ctx, shell_handle::Integer)
    sh = ctx.edge_rings[shell_handle]
    rings = Vector{Vector{Tuple{Float64, Float64}}}()
    push!(rings, sh.ring_pts)
    for h in sh.holes
        push!(rings, ctx.edge_rings[h].ring_pts)
    end
    return GI.Polygon(rings)
end

# ## Face enumeration — all minimal rings of the arrangement (Polygonizer-style)
#
# The op pipeline above extracts only the rings the op's result predicate
# selects, after dissolving interior boundaries (`unmarkDuplicateEdges`). Some
# consumers — antimeridian splitting, polygon-cut-by-line, polygonize — instead
# need the arrangement's FACES: every minimal ring of the noded linework, each
# tracing the face on its RIGHT via the half-edge face traversal (successor =
# onext ∘ sym), with per-input face locations read off the shared labels. A
# dangling edge (a line dead-end) is traversed twice by its face's ring, out
# and back — callers that need dangle-free rings filter, callers like the
# antimeridian pole seam rely on exactly this doubling.
#
# Reuses the `_OverlayEdgeRing` pipeline unchanged: the face link fills the
# same `next_result` field the op pipeline links, so `_compute_ring!` (ring
# points, kernel points, shell/hole orientation, bbox) and the hole-placement
# machinery run identically. Like the op pipeline, this consumes the graph's
# ring-linkage fields — run one extraction per `OverlayGraph`.

# Link every half-edge to its face-ring successor and build one `_OverlayEdgeRing`
# per face cycle. Requires a labelled graph (`_compute_labelling!`). Returns the
# builder context; ring handles are `1:length(ctx.edge_rings)`.
function _build_faces(m::Manifold, g::OverlayGraph{P}; exact) where {P}
    edges = g.edges
    for i in eachindex(edges)
        oe_set_next_result!(edges, i, he_onext(edges, he_sym(edges, i)))
    end
    ctx = _PolyBuilderCtx(m, edges, g.arr, exact, _MaxEdgeRing[],
                          _OverlayEdgeRing{P}[], Int32[], Int32[])
    for i in eachindex(edges)
        edges[i].edge_ring == 0 && _new_edge_ring!(ctx, i)
    end
    return ctx
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

# Build the polygons of the faces `keep(loc_a, loc_b)` selects (`keep` sees the
# raw per-input face locations): kept clockwise rings are face shells, kept
# counter-clockwise rings are cavities of kept faces, assigned to their shells
# by the same containment machinery the op pipeline uses.
function _build_face_polygons(m::Manifold, g::OverlayGraph, keep::F; exact) where {F}
    ctx = _build_faces(m, g; exact)
    for er in 1:length(ctx.edge_rings)
        keep(_face_ring_location(ctx, er, 0), _face_ring_location(ctx, er, 1)) || continue
        ring = ctx.edge_rings[er]
        ring.is_hole ? push!(ctx.free_hole_list, Int32(er)) :
                       push!(ctx.shell_list, Int32(er))
    end
    _place_free_holes!(ctx)
    return [_ring_to_polygon(ctx, sh) for sh in ctx.shell_list]
end
