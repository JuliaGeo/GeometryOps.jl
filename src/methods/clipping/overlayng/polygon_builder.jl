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
