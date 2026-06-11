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
TODO(Task 17): replace this private slice with the real machinery once
`NodeSections`/`RelateNode` land — `locate` becomes "build `NodeSections`,
push sections, `create_node`, `has_exterior_edge(node, true)`" and every
`_Ael*`/`_create_node_edges`-family helper below is deleted, with the
ported AdjacentEdgeLocatorTest cases as the regression gate.

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
# sections). `_compare_pt` is the Coordinate.compareTo port from
# node_sections.jl.
function _ns_compare_vertices(a::NodeSection, b::NodeSection)
    comp_v0 = _compare_pt(get_vertex(a, 0), get_vertex(b, 0))
    comp_v0 != 0 && return comp_v0
    return _compare_pt(get_vertex(a, 1), get_vertex(b, 1))
end

# (`_extract_unique` — the PolygonNodeConverter.extractUnique port — lived
# here until Task 16; it now comes from polygon_node_converter.jl. Its
# `compare_to` ordering reduces to `_ns_compare_vertices` for AEL sections,
# whose other compared fields are all equal.)

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

# Port of AdjacentEdgeLocator.addRing. (`_orient_ring` — the port of
# RelateGeometry.orient — lived here until Task 13; it now resides with the
# rest of the RelateGeometry port in relate_geometry.jl, together with its
# helper `_ring_is_ccw`.)
function _add_ring!(m, ring, require_cw::Bool, ring_list; exact)
    #TODO: remove repeated points?
    pts = Tuple{Float64, Float64}[_node_point(pt) for pt in GI.getpoint(ring)]
    pts = _orient_ring(m, pts, require_cw; exact)
    push!(ring_list, pts)
    return nothing
end

#==========================================================================
# RelatePointLocator (port of JTS RelatePointLocator.java)
==========================================================================#

"""
    RelatePointLocator(m::Manifold, geom; exact, is_prepared = false,
                       boundary_rule = Mod2Boundary())

Locates a point on a geometry, including mixed-type collections.
The dimension of the containing geometry element is also determined.
GeometryCollections are handled with union semantics;
i.e. the location of a point is that location of that point
on the union of the elements of the collection.

Union semantics for GeometryCollections has the following behaviours:

1. For a mixed-dimension (heterogeneous) collection a point may lie on two
   geometry elements with different dimensions. In this case the location on
   the largest-dimension element is reported.
2. For a collection with overlapping or adjacent polygons, points on polygon
   element boundaries may lie in the effective interior of the collection
   geometry.

Supports specifying the [`BoundaryNodeRule`](@ref) to use for line endpoints
(`RelateGeometry` passes its rule down here; the default matches Java's
`BoundaryNodeRule.OGC_SFS_BOUNDARY_RULE`, i.e. Mod-2).

The Java constructor signature is `RelatePointLocator(geom, isPrepared,
bnRule)`; the manifold/`exact` parameters are the only additions (consistent
with [`AdjacentEdgeLocator`](@ref)). In JTS, prepared mode swaps the
per-polygon `SimplePointInAreaLocator` for a cached
`IndexedPointInAreaLocator`; here the `is_prepared` flag is stored but both
modes currently use the direct ring loop — prepared-mode spatial indexing is
a perf follow-up (Task 22).
"""
mutable struct RelatePointLocator{M <: Manifold, E, G, BR <: BoundaryNodeRule}
    const m::M
    const exact::E
    const geom::G
    const is_prepared::Bool
    const boundary_rule::BR
    # element collections extracted from the (possibly nested-GC) input.
    # Java leaves these null when no element of that kind exists; here they
    # are simply empty. Heterogeneous GI element types force `Any` element
    # eltypes.
    const points::Set{Tuple{Float64, Float64}}
    const lines::Vector{Any}
    const polygons::Vector{Any}
    const line_boundary::LinearBoundary{BR, Tuple{Float64, Float64}}
    const is_empty::Bool
    # lazily built on the first multi-boundary point (Java: adjEdgeLocator)
    adj_edge_locator::Union{Nothing, AdjacentEdgeLocator{M, E, Tuple{Float64, Float64}}}
end

function RelatePointLocator(m::Manifold, geom; exact,
        is_prepared::Bool = false, boundary_rule::BoundaryNodeRule = Mod2Boundary())
    #-- init(geom)
    points = Set{Tuple{Float64, Float64}}()
    lines = Any[]
    polygons = Any[]
    _extract_elements!(points, lines, polygons, geom)
    # Java caches `isEmpty = geom.isEmpty()` (recursive emptiness); since
    # `extractElements` skips empty elements, the input is recursively empty
    # iff nothing was extracted.
    is_empty = isempty(points) && isempty(lines) && isempty(polygons)
    # Java builds `lineBoundary` only when lines exist; an empty
    # LinearBoundary behaves identically (no boundary, no boundary points),
    # so it is built unconditionally here.
    line_boundary = LinearBoundary(lines, boundary_rule)
    return RelatePointLocator(m, exact, geom, is_prepared, boundary_rule,
        points, lines, polygons, line_boundary, is_empty, nothing)
end

has_boundary(loc::RelatePointLocator) = has_boundary(loc.line_boundary)

# Port of RelatePointLocator.extractElements + addPoint/addLine/addPolygonal:
# trait-dispatched traversal of the (possibly nested) collection structure.
_extract_elements!(points, lines, polygons, geom) =
    _extract_elements!(points, lines, polygons, GI.trait(geom), geom)

function _extract_elements!(points, lines, polygons, trait::GI.AbstractTrait, geom)
    GI.isempty(geom) && return nothing
    if trait isa GI.PointTrait
        #-- addPoint: normalized coordinate tuples, as in LinearBoundary
        push!(points, _node_point(geom))
    elseif trait isa GI.AbstractCurveTrait
        #-- addLine (Java LinearRing extends LineString, hence AbstractCurve)
        push!(lines, geom)
    elseif trait isa Union{GI.PolygonTrait, GI.MultiPolygonTrait}
        #-- addPolygonal: whole polygonal geometry kept as one element
        push!(polygons, geom)
    elseif trait isa GI.AbstractGeometryCollectionTrait
        #-- covers GeometryCollection, MultiPoint, MultiLineString
        for g in GI.getgeom(geom)
            _extract_elements!(points, lines, polygons, g)
        end
    end
    return nothing
end

"""
    locate(loc::RelatePointLocator, p)

The location (`LOC_*` code) of point `p` relative to the locator's geometry,
under GC union semantics.
"""
locate(loc::RelatePointLocator, p) = dimloc_location(locate_with_dim(loc, p))

"""
    locate_line_end_with_dim(loc::RelatePointLocator, p)

Locates a line endpoint, as a `DL_*` dimension-location code.
In a mixed-dim GC, the line end point may also lie in an area.
In this case the area location is reported.
Otherwise, the dimloc is either `DL_LINE_BOUNDARY` or `DL_LINE_INTERIOR`,
depending on the endpoint valence and the [`BoundaryNodeRule`](@ref) in place.
"""
function locate_line_end_with_dim(loc::RelatePointLocator, p)
    #-- if a GC with areas, check for point on area
    if !isempty(loc.polygons)
        loc_poly = locate_on_polygons(loc, p, false, nothing)
        loc_poly != LOC_EXTERIOR && return dimloc_area(loc_poly)
    end
    #-- not in area, so return line end location
    return is_boundary(loc.line_boundary, p) ? DL_LINE_BOUNDARY : DL_LINE_INTERIOR
end

"""
    locate_node(loc::RelatePointLocator, p, parent_polygonal)

The location (`LOC_*` code) of a point `p` which is known to be a node of
the geometry (i.e. a vertex or on an edge). `parent_polygonal` is the
polygonal element the point is a node of (or `nothing`).
"""
locate_node(loc::RelatePointLocator, p, parent_polygonal) =
    dimloc_location(locate_node_with_dim(loc, p, parent_polygonal))

"""
    locate_node_with_dim(loc::RelatePointLocator, p, parent_polygonal)

The dimension-location (`DL_*` code) of a point `p` which is known to be a
node of the geometry.
"""
locate_node_with_dim(loc::RelatePointLocator, p, parent_polygonal) =
    locate_with_dim(loc, p, true, parent_polygonal)

"""
    locate_with_dim(loc::RelatePointLocator, p)

Computes the topological location (`DL_*` dimension-location code) of a
single point in a geometry, including the dimension of the geometry element
the point is located in (if not in the exterior). It handles both
single-element and multi-element geometries. The algorithm for multi-part
geometries takes into account the SFS Boundary Determination Rule.
"""
locate_with_dim(loc::RelatePointLocator, p) = locate_with_dim(loc, p, false, nothing)

# Private 4-argument form (Java `locateWithDim(p, isNode, parentPolygonal)`):
# `is_node` indicates the coordinate is a node (on an edge) of the geometry.
function locate_with_dim(loc::RelatePointLocator, p, is_node::Bool, parent_polygonal)
    loc.is_empty && return DL_EXTERIOR

    #=
    In a polygonal geometry a node must be on the boundary.
    (This is not the case for a mixed collection, since
    the node may be in the interior of a polygon.)
    =#
    if is_node && GI.trait(loc.geom) isa Union{GI.PolygonTrait, GI.MultiPolygonTrait}
        return DL_AREA_BOUNDARY
    end

    dim_loc = compute_dim_location(loc, p, is_node, parent_polygonal)
    return dim_loc
end

# Port of RelatePointLocator.computeDimLocation.
function compute_dim_location(loc::RelatePointLocator, p, is_node::Bool, parent_polygonal)
    #-- check dimensions in order of precedence
    if !isempty(loc.polygons)
        loc_poly = locate_on_polygons(loc, p, is_node, parent_polygonal)
        loc_poly != LOC_EXTERIOR && return dimloc_area(loc_poly)
    end
    if !isempty(loc.lines)
        loc_line = locate_on_lines(loc, p, is_node)
        loc_line != LOC_EXTERIOR && return dimloc_line(loc_line)
    end
    if !isempty(loc.points)
        loc_pt = locate_on_points(loc, p)
        loc_pt != LOC_EXTERIOR && return dimloc_point(loc_pt)
    end
    return DL_EXTERIOR
end

# Port of RelatePointLocator.locateOnPoints.
function locate_on_points(loc::RelatePointLocator, p)
    return _node_point(p) in loc.points ? LOC_INTERIOR : LOC_EXTERIOR
end

# Port of RelatePointLocator.locateOnLines.
function locate_on_lines(loc::RelatePointLocator, p, is_node::Bool)
    if is_boundary(loc.line_boundary, p)
        return LOC_BOUNDARY
    end
    #-- must be on line, in interior
    is_node && return LOC_INTERIOR

    #TODO: index the lines
    for line in loc.lines
        #-- have to check every line, since any/all may contain point
        l = locate_on_line(loc, p, is_node, line)
        l != LOC_EXTERIOR && return l
        #TODO: minor optimization - some BoundaryNodeRules can short-circuit
    end
    return LOC_EXTERIOR
end

# Port of RelatePointLocator.locateOnLine. (Java first short-circuits on the
# cached line envelope; GI geometries do not cache extents, so the check is
# skipped — perf follow-up alongside prepared-mode indexing, Task 22.)
# `is_node` is unused, as in Java (kept for signature parity).
function locate_on_line(loc::RelatePointLocator, p, is_node::Bool, line)
    #-- Java: PointLocation.isOnLine over the coordinate sequence
    n = GI.npoint(line)
    q0 = _tuple_point(GI.getpoint(line, 1))
    for i in 2:n
        q1 = _tuple_point(GI.getpoint(line, i))
        if rk_point_on_segment(loc.m, p, q0, q1; exact = loc.exact)
            return LOC_INTERIOR
        end
        q0 = q1
    end
    return LOC_EXTERIOR
end

# Port of RelatePointLocator.locateOnPolygons.
function locate_on_polygons(loc::RelatePointLocator, p, is_node::Bool, parent_polygonal)
    num_bdy = 0
    #TODO: use a spatial index on the polygons
    for polygonal in loc.polygons
        l = locate_on_polygonal(loc, p, is_node, parent_polygonal, polygonal)
        if l == LOC_INTERIOR
            return LOC_INTERIOR
        end
        if l == LOC_BOUNDARY
            num_bdy += 1
        end
    end
    if num_bdy == 1
        return LOC_BOUNDARY
    #-- check for point lying on adjacent boundaries
    elseif num_bdy > 1
        if loc.adj_edge_locator === nothing
            loc.adj_edge_locator = AdjacentEdgeLocator(loc.m, loc.geom; exact = loc.exact)
        end
        return locate(loc.adj_edge_locator, p)
    end
    return LOC_EXTERIOR
end

# Port of RelatePointLocator.locateOnPolygonal (+ getLocator): Java
# dispatches to a per-polygonal PointOnGeometryLocator (indexed when
# prepared, simple otherwise); here the SimplePointInAreaLocator ring loop
# is used for both modes (prepared-mode indexing is Task 22 territory).
function locate_on_polygonal(loc::RelatePointLocator, p, is_node::Bool, parent_polygonal, polygonal)
    if is_node && parent_polygonal === polygonal
        return LOC_BOUNDARY
    end
    return _locate_point_in_polygonal(loc.m, p, GI.trait(polygonal), polygonal; exact = loc.exact)
end

#=
Port of the SimplePointInAreaLocator logic used by `locateOnPolygonal`
(SimplePointInAreaLocator.locate → locateInGeometry → locatePointInPolygon),
with point-in-ring routed through the kernel: shell first, then standard
even-odd composition over the holes. (The Java envelope short-circuit is
skipped, as in `locate_on_line`.)
=#
function _locate_point_in_polygonal(m, p, ::GI.PolygonTrait, poly; exact)
    GI.isempty(poly) && return LOC_EXTERIOR
    shell_loc = rk_point_in_ring(m, p, GI.getexterior(poly); exact)
    shell_loc != LOC_INTERIOR && return shell_loc
    #-- now test if the point lies in or on the holes
    for hole in GI.gethole(poly)
        hole_loc = rk_point_in_ring(m, p, hole; exact)
        hole_loc == LOC_BOUNDARY && return LOC_BOUNDARY
        hole_loc == LOC_INTERIOR && return LOC_EXTERIOR
        #-- if in EXTERIOR of this hole keep checking the other ones
    end
    return LOC_INTERIOR
end

function _locate_point_in_polygonal(m, p, ::GI.MultiPolygonTrait, mp; exact)
    for poly in GI.getgeom(mp)
        l = _locate_point_in_polygonal(m, p, GI.trait(poly), poly; exact)
        l != LOC_EXTERIOR && return l
    end
    return LOC_EXTERIOR
end
