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
## LinearBoundary (port of JTS LinearBoundary.java)
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

function LinearBoundary(m::Manifold, lines, rule::BoundaryNodeRule)
    # assert: dim(geom) == 1
    vertex_degree = _compute_boundary_points(m, lines)
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

function _compute_boundary_points(m::Manifold, lines)
    vertex_degree = Dict{_kernel_point_type(m), Int}()
    for line in lines
        n = GI.npoint(line)
        n == 0 && continue
        _add_endpoint!(_to_kernel_point(m, GI.getpoint(line, 1)), vertex_degree)
        _add_endpoint!(_to_kernel_point(m, GI.getpoint(line, n)), vertex_degree)
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
## AdjacentEdgeLocator (port of JTS AdjacentEdgeLocator.java)
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
    ring_list = Vector{_kernel_point_type(m)}[]
    _ael_init!(m, ring_list, geom; exact)
    return AdjacentEdgeLocator(m, exact, ring_list)
end

"""
    locate(ael::AdjacentEdgeLocator, p)

Location (`LOC_INTERIOR` or `LOC_BOUNDARY`) of point `p`, which must lie on
at least one polygon edge of the locator's geometry, under union semantics.
"""
function locate(ael::AdjacentEdgeLocator, p)
    pt = _node_point(p)
    sections = NodeSections(vertex_node(pt))
    for ring in ael.ring_list
        _add_sections!(ael, pt, ring, sections)
    end
    node = create_node(ael.m, sections; exact = ael.exact)
    #node.finish(false, false);
    return has_exterior_edge(node, true) ? LOC_BOUNDARY : LOC_INTERIOR
end

# Port of AdjacentEdgeLocator.addSections.
function _add_sections!(ael::AdjacentEdgeLocator, p, ring, sections::NodeSections)
    for i in 1:(length(ring) - 1)
        p0 = ring[i]
        pnext = ring[i + 1]

        if _equals2(p, pnext)
            #-- segment final point is assigned to next segment
            continue
        elseif _equals2(p, p0)
            iprev = i > 1 ? i - 1 : length(ring) - 1
            pprev = ring[iprev]
            add_node_section!(sections, _create_section(ael, p, pprev, pnext))
        elseif rk_point_on_segment(ael.m, p, p0, pnext; exact = ael.exact)
            add_node_section!(sections, _create_section(ael, p, p0, pnext))
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
    pts = _to_kernel_points(m, ring)
    pts = _orient_ring(m, pts, require_cw; exact)
    push!(ring_list, pts)
    return nothing
end

#==========================================================================
## RelatePointLocator (port of JTS RelatePointLocator.java)
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
with [`AdjacentEdgeLocator`](@ref)). As in JTS, prepared mode swaps the
per-polygon `SimplePointInAreaLocator` ring loop for a cached
[`IndexedPointInAreaLocator`](@ref) (indexed_point_in_area.jl), created
lazily on the first use per polygonal element (Task 22). Unprepared mode
deviates from Java (which keys indexing on `isPrepared` alone): the first
query on a polygonal element uses the direct ring loop, but repeat queries
build and reuse the indexed locator — one O(n) scan beats an O(n) index
build, while the many area-vertex locations of a multi-element relate
amortize the index (see `locate_on_polygonal`).
"""
mutable struct RelatePointLocator{M <: Manifold, E, G, BR <: BoundaryNodeRule, P}
    const m::M
    const exact::E
    const geom::G
    const is_prepared::Bool
    const boundary_rule::BR
    # element collections extracted from the (possibly nested-GC) input.
    # Java leaves these null when no element of that kind exists; here they
    # are simply empty. Heterogeneous GI element types force `Any` element
    # eltypes. `P` is the manifold's kernel point type (Phase 3).
    const points::Set{P}
    const lines::Vector{Any}
    const polygons::Vector{Any}
    const line_boundary::LinearBoundary{BR, P}
    const is_empty::Bool
    # per-polygonal-element indexed locators, created lazily by
    # `_get_poly_locator` (Java: polyLocator, filled by getLocator).
    # Prepared mode fills an entry on its first query; unprepared mode on
    # its second (see `locate_on_polygonal`).
    const poly_locator::Vector{Union{Nothing, IndexedPointInAreaLocator{M, E}}}
    # unprepared mode: queries seen per polygonal element, driving the lazy
    # index heuristic above
    const poly_query_count::Vector{Int32}
    # lazily built on the first multi-boundary point (Java: adjEdgeLocator)
    adj_edge_locator::Union{Nothing, AdjacentEdgeLocator{M, E, P}}
end

function RelatePointLocator(m::Manifold, geom; exact,
        is_prepared::Bool = false, boundary_rule::BoundaryNodeRule = Mod2Boundary())
    #-- init(geom)
    P = _kernel_point_type(m)
    points = Set{P}()
    lines = Any[]
    polygons = Any[]
    _extract_elements!(m, points, lines, polygons, geom)
    # Java caches `isEmpty = geom.isEmpty()` (recursive emptiness); since
    # `extractElements` skips empty elements, the input is recursively empty
    # iff nothing was extracted.
    is_empty = isempty(points) && isempty(lines) && isempty(polygons)
    # Java builds `lineBoundary` only when lines exist; an empty
    # LinearBoundary behaves identically (no boundary, no boundary points),
    # so it is built unconditionally here.
    line_boundary = LinearBoundary(m, lines, boundary_rule)
    # Java allocates `polyLocator` for both modes (Simple/Indexed); here both
    # modes may cache indexed locator objects (unprepared lazily, on repeat
    # queries), so it is allocated unconditionally.
    poly_locator = Vector{Union{Nothing, IndexedPointInAreaLocator{typeof(m), typeof(exact)}}}(
        nothing, length(polygons))
    poly_query_count = zeros(Int32, length(polygons))
    #-- P cannot be inferred from the `nothing` adj_edge_locator, so spell out
    #-- every type parameter
    return RelatePointLocator{typeof(m), typeof(exact), typeof(geom),
            typeof(boundary_rule), P}(
        m, exact, geom, is_prepared, boundary_rule,
        points, lines, polygons, line_boundary, is_empty, poly_locator,
        poly_query_count, nothing)
end

has_boundary(loc::RelatePointLocator) = has_boundary(loc.line_boundary)

# Port of RelatePointLocator.extractElements + addPoint/addLine/addPolygonal:
# trait-dispatched traversal of the (possibly nested) collection structure.
_extract_elements!(m, points, lines, polygons, geom) =
    _extract_elements!(m, points, lines, polygons, GI.trait(geom), geom)

function _extract_elements!(m, points, lines, polygons, ::GI.PointTrait, geom)
    GI.isempty(geom) && return nothing
    #-- addPoint: normalized kernel points, as in LinearBoundary
    push!(points, _to_kernel_point(m, geom))
    return nothing
end
function _extract_elements!(m, points, lines, polygons, ::GI.AbstractCurveTrait, geom)
    GI.isempty(geom) && return nothing
    #-- addLine (Java LinearRing extends LineString, hence AbstractCurve)
    push!(lines, geom)
    return nothing
end
function _extract_elements!(m, points, lines, polygons,
        ::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, geom)
    GI.isempty(geom) && return nothing
    #-- addPolygonal: whole polygonal geometry kept as one element
    push!(polygons, geom)
    return nothing
end
function _extract_elements!(m, points, lines, polygons,
        ::GI.AbstractGeometryCollectionTrait, geom)
    GI.isempty(geom) && return nothing
    #-- covers GeometryCollection, MultiPoint, MultiLineString
    for g in GI.getgeom(geom)
        _extract_elements!(m, points, lines, polygons, g)
    end
    return nothing
end
_extract_elements!(m, points, lines, polygons, ::GI.AbstractTrait, geom) = nothing

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

# Port of RelatePointLocator.locateOnLine, including Java's short-circuit on
# the cached line envelope (the lines come from the RelateGeometry wrapper
# tree, which carries a stored extent on every linework element).
# `is_node` is unused, as in Java (kept for signature parity).
function locate_on_line(loc::RelatePointLocator, p, is_node::Bool, line)
    #-- Java: lineEnv.intersects(p) short-circuit (p is already a kernel point)
    pt_ext = _kernel_point_box(p)
    if !Extents.intersects(rk_interaction_bounds(loc.m, line), pt_ext)
        return LOC_EXTERIOR
    end
    #-- Java: PointLocation.isOnLine over the coordinate sequence
    n = GI.npoint(line)
    q0 = _to_kernel_point(loc.m, GI.getpoint(line, 1))
    for i in 2:n
        q1 = _to_kernel_point(loc.m, GI.getpoint(line, i))
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
    for i in eachindex(loc.polygons)
        l = locate_on_polygonal(loc, p, is_node, parent_polygonal, i)
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

# Queries a polygonal element absorbs via the direct ring loop before its
# IndexedPointInAreaLocator is built. Both costs scale with the element's
# segment count, so one threshold fits all sizes: an unsorted index build
# costs ~10-13 ring scans (measured on Natural Earth coastlines), making
# the worst-case regret of switching at 8 about one build. Real relates are
# bimodal — a handful of queries (barely-touching neighbors, where indexing
# never pays) or hundreds (one area-vertex location per polygon element of
# the other geometry), so the threshold rarely sits near the break-even.
const _LAZY_INDEX_QUERY_THRESHOLD = Int32(8)

# Port of RelatePointLocator.locateOnPolygonal: Java dispatches to a
# per-polygonal PointOnGeometryLocator — a cached IndexedPointInAreaLocator
# when prepared, a SimplePointInAreaLocator otherwise. Prepared mode does
# the same here (Task 22). Unprepared mode deviates from Java: the first
# query on an element uses the direct SimplePointInAreaLocator ring loop
# (one O(n) scan beats an O(n) index build + query), but repeat queries —
# e.g. one area-vertex location per polygon element of the other geometry
# in a multipolygon/multipolygon relate — build and amortize the index.
function locate_on_polygonal(loc::RelatePointLocator, p, is_node::Bool, parent_polygonal, index::Int)
    polygonal = loc.polygons[index]
    if is_node && parent_polygonal === polygonal
        return LOC_BOUNDARY
    end
    #-- the RayCrossingCounter horizontal-ray sweep is coordinate-plane
    #-- logic (as is all of JTS), so a future non-planar kernel falls
    #-- through to its own rk_point_in_ring even when prepared
    if loc.m isa Planar
        use_index = loc.is_prepared
        if !use_index
            count = (loc.poly_query_count[index] += Int32(1))
            use_index = count > _LAZY_INDEX_QUERY_THRESHOLD
        end
        if use_index
            return locate(_get_poly_locator(loc, index), p)
        end
    end
    return _locate_point_in_polygonal(loc.m, p, GI.trait(polygonal), polygonal; exact = loc.exact)
end

# Port of RelatePointLocator.getLocator (indexed arm): lazily create and
# cache the indexed locator for polygonal element `index`. Prepared mode
# pays for the midpoint-sorted layout (build once, query forever); the
# unprepared lazy index skips the sort, which dominates the build cost
# (see `SortedPackedIntervalRTree`).
function _get_poly_locator(loc::RelatePointLocator, index::Int)
    locator = loc.poly_locator[index]
    if locator === nothing
        locator = IndexedPointInAreaLocator(loc.m, loc.polygons[index];
            exact = loc.exact, sort_leaves = loc.is_prepared)
        loc.poly_locator[index] = locator
    end
    return locator
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
