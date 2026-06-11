# # RelateNG point location
#
# Point-location machinery for RelateNG. This file holds the ports of three
# small, tightly coupled JTS classes, in this order (JTS file boundaries
# preserved as clearly marked sections):
#
# 1. `LinearBoundary`        (JTS LinearBoundary.java)
# 2. `AdjacentEdgeLocator`   (JTS AdjacentEdgeLocator.java) — Task 11
# 3. `RelatePointLocator`    (JTS RelatePointLocator.java)  — Task 12

#==========================================================================
# LinearBoundary (port of JTS LinearBoundary.java)
==========================================================================#

"""
    LinearBoundary(lines, rule::BoundaryNodeRule)

Determines the boundary points of a linear geometry, using a
[`BoundaryNodeRule`](@ref). `lines` is an iterable of linestrings
(any GeoInterface linestring-like geometries); the endpoint degree of
every line endpoint is counted and the rule decides which degrees are
boundary points.

Coordinate keys are normalized via `_node_point` (kernel.jl): exact
`(Float64, Float64)` tuples with signed zeros normalized (`-0.0 → +0.0`),
so lookups here agree with the `NodeKey` vertex-node identity from the
kernel (Task 7) under Dict bit-pattern hashing.

Faithful to Java: only *empty* lines are skipped. Closed lines are NOT
special-cased — a closed line contributes degree 2 to its closure vertex
(both endpoints coincide), which is never a boundary under the Mod-2 or
monovalent rules but would be under e.g. the endpoint rule.
"""
struct LinearBoundary{BR <: BoundaryNodeRule, P}
    vertex_degree::Dict{P, Int}
    has_boundary::Bool
    rule::BR
end

function LinearBoundary(lines, rule::BoundaryNodeRule)
    # assert: dim(geom) == 1
    vertex_degree = _compute_boundary_points(lines)
    has_boundary = _check_boundary(vertex_degree, rule)
    return LinearBoundary(vertex_degree, has_boundary, rule)
end

function _check_boundary(vertex_degree::Dict, rule::BoundaryNodeRule)
    for degree in values(vertex_degree)
        if is_in_boundary(rule, degree)
            return true
        end
    end
    return false
end

has_boundary(lb::LinearBoundary) = lb.has_boundary

function is_boundary(lb::LinearBoundary, pt)
    key = _node_point(pt)
    haskey(lb.vertex_degree, key) || return false
    degree = lb.vertex_degree[key]
    return is_in_boundary(lb.rule, degree)
end

function _compute_boundary_points(lines)
    vertex_degree = Dict{Tuple{Float64, Float64}, Int}()
    for line in lines
        n = GI.npoint(line)
        n == 0 && continue
        _add_endpoint!(_node_point(GI.getpoint(line, 1)), vertex_degree)
        _add_endpoint!(_node_point(GI.getpoint(line, n)), vertex_degree)
    end
    return vertex_degree
end

function _add_endpoint!(p, degree::Dict)
    dim = get(degree, p, 0)
    dim += 1
    degree[p] = dim
    return nothing
end

#==========================================================================
# AdjacentEdgeLocator (port of JTS AdjacentEdgeLocator.java)
==========================================================================#

"""
    AdjacentEdgeLocator(m::Manifold, geom; exact)

Determines the location for a point which is known to lie on at least one
edge of a set of polygons. This provides the union-semantics for determining
point location in a GeometryCollection, which may have polygons with
adjacent edges which are effectively in the interior of the geometry.
Note that it is also possible to have adjacent edges which lie on the
boundary of the geometry (e.g. a polygon contained within another polygon
with adjacent edges).

The manifold `m` and the `exact` flag are stored in the struct (rather than
threaded through every call) for consistency with how `RelateGeometry`
holds them (Task 13); [`locate`](@ref locate(::AdjacentEdgeLocator, ::Any))
uses the stored values for all kernel queries.

The Java constructor signature is `AdjacentEdgeLocator(Geometry geom)`; the
manifold/exact parameters are the only additions.
"""
struct AdjacentEdgeLocator{M <: Manifold, E, P}
    m::M
    exact::E
    ring_list::Vector{Vector{P}}
end

function AdjacentEdgeLocator(m::Manifold, geom; exact)
    ring_list = Vector{Tuple{Float64, Float64}}[]
    _ael_init!(m, ring_list, geom; exact)
    return AdjacentEdgeLocator(m, exact, ring_list)
end

"""
    locate(ael::AdjacentEdgeLocator, p)

Location (`LOC_INTERIOR` or `LOC_BOUNDARY`) of point `p`, which must lie on
at least one polygon edge of the locator's geometry, under union semantics.
"""
function locate(ael::AdjacentEdgeLocator{M, E, P}, p) where {M, E, P}
    pt = _node_point(p)
    # Stand-in for `NodeSections(p)` (Task 15): a plain section list; node
    # assembly is specialized in `_create_node_edges` below.
    sections = NodeSection{P, Nothing}[]
    for ring in ael.ring_list
        _add_sections!(ael, pt, ring, sections)
    end
    # Java: `RelateNode node = sections.createNode();
    #        return node.hasExteriorEdge(true) ? BOUNDARY : INTERIOR;`
    node_edges = _create_node_edges(ael.m, vertex_node(pt), sections; exact = ael.exact)
    return _node_has_exterior_edge(node_edges) ? LOC_BOUNDARY : LOC_INTERIOR
end

# Port of AdjacentEdgeLocator.addSections.
function _add_sections!(ael::AdjacentEdgeLocator, p, ring, sections)
    for i in 1:(length(ring) - 1)
        p0 = ring[i]
        pnext = ring[i + 1]

        if _equals2(p, pnext)
            #-- segment final point is assigned to next segment
            continue
        elseif _equals2(p, p0)
            iprev = i > 1 ? i - 1 : length(ring) - 1
            pprev = ring[iprev]
            push!(sections, _create_section(ael, p, pprev, pnext))
        elseif rk_point_on_segment(ael.m, p, p0, pnext; exact = ael.exact)
            push!(sections, _create_section(ael, p, p0, pnext))
        end
    end
    return nothing
end

# Port of AdjacentEdgeLocator.createSection. (The Java prints a debug warning
# for zero-length section segments; here they are simply constructed — they
# only arise from invalid rings with repeated points.)
function _create_section(::AdjacentEdgeLocator, p, prev, next)
    return NodeSection(true, DIM_A, Int32(1), Int32(0), nothing, false, prev, vertex_node(p), next)
end

#=
Node-wheel construction, standing in for the Java pipeline
`NodeSections.createNode()` → `RelateNode` until the full node-topology
machinery lands (Tasks 15–17). This is a faithful slice of
`NodeSections.java` / `PolygonNodeConverter.java` / `RelateNode.java` /
`RelateEdge.java`, specialized to what AdjacentEdgeLocator can produce:
every section is an area (`DIM_A`) corner of geometry A at the node, in
canonical orientation (CW shells / CCW holes, i.e. polygon interior on the
right of travel `v0 → node → v1`), with `id = 1` and `ring_id = 0`.

The specialization is exact:

- `NodeSections.prepareSections` sorts by `NodeSection.compareTo`; with
  `is_a`/`dim`/`id`/`ring_id` all equal that reduces to comparing the edge
  vertices (`_ns_compare_vertices`).
- Since every section reports the same polygon (`id = 1`), the sections are
  routed through `PolygonNodeConverter.convert` whenever there are at least
  two of them. The converter sorts by the angle of the entering edge
  (`EdgeAngleComparator`), drops exact duplicates (`extractUnique`), and —
  because every AEL section is a shell section (`ring_id = 0`) — rewrites
  each section to itself (`convertShellAndHoles` finds no holes). So the
  conversion reduces to: stable angle sort + dedupe. (For a single section
  Java skips the converter; deduping a singleton is the identity, so the
  same code path is used here.)
- `RelateNode.addEdges` only ever sees `DIM_A` sections of geometry A, so
  the edge labels collapse to one (left, right, line) location triple
  (`_AelEdge`), and `RelateEdge.merge` reduces to its area-area branch
  (`mergeDimEdgeLoc` is a no-op between two area edges).

Note the construction is order-dependent (`updateEdgesInArea` only marks
edges already present in the wheel), which is why the Java processing order
is reproduced exactly. All angle comparisons go through
`rk_compare_edge_dir` around the symbolic node.
=#

# Slice of RelateEdge restricted to area edges of a single geometry:
# `dim` is always `DIM_A` and `is_a` is always `true`, so only the
# direction point and the location triple remain.
mutable struct _AelEdge{P}
    dir_pt::P
    loc_left::Int8
    loc_right::Int8
    loc_line::Int8
end

# Port of NodeSections.createNode, specialized as described above. Returns
# the wheel of edges around `node` in CCW order.
function _create_node_edges(m, node::NodeKey{P}, sections::Vector{<:NodeSection}; exact) where {P}
    edges = _AelEdge{P}[]
    isempty(sections) && return edges
    #-- NodeSections.prepareSections
    sort!(sections; lt = (a, b) -> _ns_compare_vertices(a, b) < 0)
    #-- PolygonNodeConverter.convert: EdgeAngleComparator sort (stable, so
    #-- equal-angle sections keep their prepareSections order) + extractUnique
    sort!(sections; alg = MergeSort,
        lt = (a, b) -> rk_compare_edge_dir(m, node, get_vertex(a, 0), get_vertex(b, 0); exact) < 0)
    unique_sections = _extract_unique(sections)
    for ns in unique_sections
        _add_edges!(m, node, edges, ns; exact)
    end
    return edges
end

# Specialization of NodeSection.compareTo for AdjacentEdgeLocator sections:
# `is_a`/`dim`/`id`/`ring_id` are identical across sections, so only the
# edge-vertex comparison remains (no `nothing` vertices occur for area
# sections).
function _ns_compare_vertices(a::NodeSection, b::NodeSection)
    comp_v0 = _compare_pt(get_vertex(a, 0), get_vertex(b, 0))
    comp_v0 != 0 && return comp_v0
    return _compare_pt(get_vertex(a, 1), get_vertex(b, 1))
end

# JTS Coordinate.compareTo: lexicographic on (x, y).
function _compare_pt(p, q)
    GI.x(p) < GI.x(q) && return -1
    GI.x(p) > GI.x(q) && return 1
    GI.y(p) < GI.y(q) && return -1
    GI.y(p) > GI.y(q) && return 1
    return 0
end

# Port of PolygonNodeConverter.extractUnique: drop consecutive duplicate
# sections (assumes the list is sorted so duplicates are adjacent).
function _extract_unique(sections::Vector{S}) where {S <: NodeSection}
    unique_sections = S[sections[1]]
    for ns in sections
        if _ns_compare_vertices(last(unique_sections), ns) != 0
            push!(unique_sections, ns)
        end
    end
    return unique_sections
end

# Port of RelateNode.addEdges(NodeSection), area case (the only case here).
function _add_edges!(m, node::NodeKey, edges::Vector{<:_AelEdge}, ns::NodeSection; exact)
    #-- assumes node edges have CW orientation (as per JTS norm)
    #-- entering edge - interior on L
    e0 = _add_area_edge!(m, node, edges, get_vertex(ns, 0), false; exact)
    #-- exiting edge - interior on R
    e1 = _add_area_edge!(m, node, edges, get_vertex(ns, 1), true; exact)
    # Zero-length edges are skipped (Java returns null from addEdge; they
    # only arise from invalid rings with repeated points).
    (e0 === nothing || e1 === nothing) && return nothing

    index0 = findfirst(e -> e === e0, edges)
    index1 = findfirst(e -> e === e1, edges)
    _update_edges_in_area!(edges, index0, index1)
    _update_if_area_prev!(edges, index0)
    _update_if_area_next!(edges, index1)
    return nothing
end

# Port of RelateNode.updateEdgesInArea: mark every edge strictly between
# the entering and exiting edge (in CCW order) as area interior.
function _update_edges_in_area!(edges, index_from, index_to)
    index = _next_index(edges, index_from)
    while index != index_to
        _set_area_interior!(edges[index])
        index = _next_index(edges, index)
    end
    return nothing
end

# Port of RelateNode.updateIfAreaPrev.
function _update_if_area_prev!(edges, index)
    index_prev = _prev_index(edges, index)
    if edges[index_prev].loc_left == LOC_INTERIOR
        _set_area_interior!(edges[index])
    end
    return nothing
end

# Port of RelateNode.updateIfAreaNext.
function _update_if_area_next!(edges, index)
    index_next = _next_index(edges, index)
    if edges[index_next].loc_right == LOC_INTERIOR
        _set_area_interior!(edges[index])
    end
    return nothing
end

#=
Port of RelateNode.addEdge restricted to area edges (addAreaEdge): adds or
merges an edge to the node wheel, keeping the wheel sorted by CCW angle
with the positive X-axis. `is_forward` is the direction of the edge.
Returns the created or merged edge for this point, or `nothing` for a
zero-length (malformed) edge.
=#
function _add_area_edge!(m, node::NodeKey, edges::Vector{<:_AelEdge}, dir_pt, is_forward::Bool; exact)
    #-- check for well-formed edge - skip zero-len input
    _equals2(node.pt, dir_pt) && return nothing

    insert_index = 0
    for (i, e) in pairs(edges)
        comp = rk_compare_edge_dir(m, node, e.dir_pt, dir_pt; exact)
        if comp == 0
            _merge_area!(e, is_forward)
            return e
        end
        if comp == 1
            #-- found further edge, so insert a new edge at this position
            insert_index = i
            break
        end
    end
    #-- add a new edge (RelateEdge.create / setLocationsArea)
    e = _AelEdge(dir_pt,
        is_forward ? LOC_EXTERIOR : LOC_INTERIOR,   # left
        is_forward ? LOC_INTERIOR : LOC_EXTERIOR,   # right
        LOC_BOUNDARY)                               # line
    if insert_index == 0
        #-- add edge at end of list
        push!(edges, e)
    else
        #-- add edge before higher edge found
        insert!(edges, insert_index, e)
    end
    return e
end

# Port of RelateEdge.merge, area-into-area case (the only case here:
# mergeDimEdgeLoc is a no-op between two area edges, and the on-location
# stays BOUNDARY).
function _merge_area!(e::_AelEdge, is_forward::Bool)
    loc_left = is_forward ? LOC_EXTERIOR : LOC_INTERIOR
    loc_right = is_forward ? LOC_INTERIOR : LOC_EXTERIOR
    #-- mergeSideLocation: INTERIOR takes precedence over EXTERIOR
    if e.loc_left != LOC_INTERIOR
        e.loc_left = loc_left
    end
    if e.loc_right != LOC_INTERIOR
        e.loc_right = loc_right
    end
    return nothing
end

# Port of RelateEdge.setAreaInterior (single-geometry form).
function _set_area_interior!(e::_AelEdge)
    e.loc_left = LOC_INTERIOR
    e.loc_right = LOC_INTERIOR
    e.loc_line = LOC_INTERIOR
    return nothing
end

# Ports of RelateNode.prevIndex / nextIndex (1-based).
_prev_index(edges, i) = i > 1 ? i - 1 : length(edges)
_next_index(edges, i) = i >= length(edges) ? 1 : i + 1

# Port of RelateNode.hasExteriorEdge(true). (An empty wheel — which the
# `locate` precondition rules out — yields `false`, as in Java.)
_node_has_exterior_edge(edges) =
    any(e -> e.loc_left == LOC_EXTERIOR || e.loc_right == LOC_EXTERIOR, edges)

# Port of AdjacentEdgeLocator.init + addRings: collect the polygon rings of
# the (possibly collection) geometry as canonically oriented coordinate
# vectors. (Java leaves `ringList` null for an empty geometry; here it just
# stays empty.)
function _ael_init!(m, ring_list, geom; exact)
    _add_rings!(m, GI.trait(geom), geom, ring_list; exact)
    return nothing
end

_add_rings!(m, geom, ring_list; exact) =
    _add_rings!(m, GI.trait(geom), geom, ring_list; exact)

function _add_rings!(m, ::GI.PolygonTrait, poly, ring_list; exact)
    shell = GI.getexterior(poly)
    _add_ring!(m, shell, true, ring_list; exact)
    for hole in GI.gethole(poly)
        _add_ring!(m, hole, false, ring_list; exact)
    end
    return nothing
end

#-- recurse through collections (Java `instanceof GeometryCollection` covers
#-- MultiPolygon etc.; multi-point/line elements fall through to the no-op)
function _add_rings!(m, ::GI.AbstractGeometryCollectionTrait, geom, ring_list; exact)
    for g in GI.getgeom(geom)
        _add_rings!(m, g, ring_list; exact)
    end
    return nothing
end

_add_rings!(m, ::GI.AbstractTrait, geom, ring_list; exact) = nothing

# Port of AdjacentEdgeLocator.addRing.
function _add_ring!(m, ring, require_cw::Bool, ring_list; exact)
    #TODO: remove repeated points?
    pts = Tuple{Float64, Float64}[_node_point(pt) for pt in GI.getpoint(ring)]
    pts = _orient_ring(m, pts, require_cw; exact)
    push!(ring_list, pts)
    return nothing
end

# Port of RelateGeometry.orient (static; needed here first, reused by the
# RelateGeometry port in Task 13): coordinate vector of `pts` with the
# requested orientation, reversing a copy only if needed.
function _orient_ring(m, pts::Vector, orient_cw::Bool; exact)
    is_flipped = orient_cw == _ring_is_ccw(m, pts; exact)
    return is_flipped ? reverse(pts) : pts
end

#=
Port of JTS `Orientation.isCCW(CoordinateSequence)` with the orientation
index routed through the kernel (`rk_orient`): whether the closed `ring`
(repeated end point required) is counterclockwise. Returns `false` for flat
or degenerate rings. The algorithm finds the highest point reached by an
upward segment, then the subsequent downward segment, and decides from the
"cap" they form — using only one exact orientation test, so it is robust
(unlike a floating signed-area sum).
=#
function _ring_is_ccw(m, ring::Vector; exact)
    # number of points without closing endpoint
    npts = length(ring) - 1
    # return default value if ring is flat
    npts < 3 && return false
    pt(i) = ring[i + 1]   # 0-based access, mirroring the Java indexing

    # Find first highest point after a lower point, if one exists
    # (e.g. a rising segment). If one does not exist, i_up_hi remains 0
    # and the ring must be flat.
    up_hi_pt = pt(0)
    prev_y = GI.y(up_hi_pt)
    up_low_pt = up_hi_pt   # only read when i_up_hi != 0, i.e. after assignment
    i_up_hi = 0
    for i in 1:npts
        py = GI.y(pt(i))
        # if segment is upwards and endpoint is higher, record it
        if py > prev_y && py >= GI.y(up_hi_pt)
            up_hi_pt = pt(i)
            i_up_hi = i
            up_low_pt = pt(i - 1)
        end
        prev_y = py
    end
    # check if ring is flat and return default value if so
    i_up_hi == 0 && return false

    # Find the next lower point after the high point (e.g. a falling
    # segment). This must exist since the ring is not flat.
    i_down_low = i_up_hi
    while true
        i_down_low = (i_down_low + 1) % npts
        (i_down_low != i_up_hi && GI.y(pt(i_down_low)) == GI.y(up_hi_pt)) || break
    end

    down_low_pt = pt(i_down_low)
    i_down_hi = i_down_low > 0 ? i_down_low - 1 : npts - 1
    down_hi_pt = pt(i_down_hi)

    if _equals2(up_hi_pt, down_hi_pt)
        # the high point is on a "pointed cap": its orientation decides.
        # Degenerate A-B-A caps (coincident segments / < 3 distinct points)
        # have orientation 0 and return false.
        (_equals2(up_low_pt, up_hi_pt) || _equals2(down_low_pt, up_hi_pt) ||
            _equals2(up_low_pt, down_low_pt)) && return false
        return rk_orient(m, up_low_pt, up_hi_pt, down_low_pt; exact) > 0
    else
        # flat cap - direction of flat top determines orientation
        del_x = GI.x(down_hi_pt) - GI.x(up_hi_pt)
        return del_x < 0
    end
end
