# # RelateNG input geometry facade
#
# Ports of two tightly coupled JTS classes (JTS file boundaries preserved as
# marked sections):
#
# 1. `RelateGeometry`       (JTS RelateGeometry.java)       — Task 13
# 2. `RelateSegmentString`  (JTS RelateSegmentString.java)  — Task 13
#
# `RelateGeometry` wraps one input geometry of a relate operation and caches
# the metadata the topology computer needs: dimension analysis, emptiness,
# extent, plus lazy access to unique points and a `RelatePointLocator`.
# `RelateSegmentString` models one linear edge (line or ring) extracted from
# it.

#==========================================================================
## RelateGeometry (port of JTS RelateGeometry.java)
==========================================================================#

# `GEOM_A`/`GEOM_B` (JTS `RelateGeometry.GEOM_A/GEOM_B`) are defined in
# relate_predicates.jl (Task 4) and are not redefined here.

# Port of RelateGeometry.name(isA).
geom_name(is_a::Bool) = is_a ? "A" : "B"

"""
    RelateGeometry(m::Manifold, geom; exact, is_prepared = false,
                   boundary_rule = Mod2Boundary())

The input-geometry facade of RelateNG: wraps one of the two operand
geometries and caches its metadata — recursive emptiness, extent,
dimension analysis (`has_points`/`has_lines`/`has_areas`), zero-length-line
degeneracy — plus lazily created unique points and a
[`RelatePointLocator`](@ref).

The Java constructor signature is `RelateGeometry(Geometry input, boolean
isPrepared, BoundaryNodeRule bnRule)`; the manifold/`exact` parameters are
the only additions (consistent with [`RelatePointLocator`](@ref)). Where the
Java caches `geomEnv = input.getEnvelopeInternal()`, here `extent` is the
union of the interaction bounds (`rk_interaction_bounds`) of the non-empty
elements, or `nothing` if the geometry is empty.
"""
mutable struct RelateGeometry{M <: Manifold, E, G, BR <: BoundaryNodeRule, X, P}
    const m::M
    const exact::E
    const geom::G
    const is_prepared::Bool
    const boundary_rule::BR
    const extent::X
    const dim::Int8
    const has_points::Bool
    const has_lines::Bool
    const has_areas::Bool
    const is_line_zero_len::Bool
    const is_geom_empty::Bool
    # id counter for extracted elements (Java `elementId`); never reset, so
    # repeated extraction (prepared mode) keeps producing distinct ids.
    element_id::Int32
    # lazy caches. `P` is the manifold's kernel point type (Phase 3): the
    # coordinate type of every node point and segment-string vertex.
    unique_points::Union{Nothing, Set{P}}
    locator::Union{Nothing, RelatePointLocator{M, E, G, BR, P}}
end

function RelateGeometry(m::Manifold, geom; exact,
        is_prepared::Bool = false, boundary_rule::BoundaryNodeRule = Mod2Boundary())
    #-- cache geometry metadata
    is_geom_empty = _relate_is_empty(geom)
    #-- Hold an extent-cached wrapper tree instead of the raw input: one
    #-- coordinate pass here makes every downstream extent consult (the
    #-- engine's envelope checks, extraction filters, line-end walks, point
    #-- locator short-circuits) O(1) — the GI equivalent of the envelope
    #-- cache JTS carries on every Geometry. Coordinates are never copied:
    #-- the wrappers share the original linework objects.
    geom = is_geom_empty ? geom : _relate_cache_extents(m, geom)
    extent = _relate_extent(m, geom)
    dim = _geom_dimension(geom)
    dim, has_points, has_lines, has_areas = _analyze_dimensions(geom, dim, is_geom_empty)
    is_line_zero_len = _is_zero_length_line(geom, dim)
    #-- P (the kernel point type) cannot be inferred from the `nothing` lazy
    #-- caches, so spell out every type parameter
    P = _kernel_point_type(m)
    return RelateGeometry{typeof(m), typeof(exact), typeof(geom), typeof(boundary_rule),
            typeof(extent), P}(
        m, exact, geom, is_prepared, boundary_rule, extent,
        dim, has_points, has_lines, has_areas, is_line_zero_len, is_geom_empty,
        Int32(0), nothing, nothing)
end

# Recursive emptiness (Java `Geometry.isEmpty()`): a collection is empty iff
# every element is empty. `GI.isempty` is not recursive for collections.
_relate_is_empty(geom) = _relate_is_empty(GI.trait(geom), geom)
function _relate_is_empty(::GI.AbstractGeometryCollectionTrait, geom)
    for g in GI.getgeom(geom)
        _relate_is_empty(g) || return false
    end
    return true
end
_relate_is_empty(::GI.AbstractTrait, geom) = GI.isempty(geom)

# Equivalent of Java `Geometry.getEnvelopeInternal()` as interaction bounds:
# the union of `rk_interaction_bounds` over the non-empty atomic elements
# (empty elements contribute nothing in Java too), or `nothing` if the
# geometry is (recursively) empty.
_relate_extent(m::Manifold, geom) = _relate_extent(m, GI.trait(geom), geom)
function _relate_extent(m::Manifold, ::GI.AbstractGeometryCollectionTrait, geom)
    ext = nothing
    for g in GI.getgeom(geom)
        e = _relate_extent(m, g)
        e === nothing && continue
        ext = ext === nothing ? e : Extents.union(ext, e)
    end
    return ext
end
function _relate_extent(m::Manifold, ::GI.AbstractTrait, geom)
    GI.isempty(geom) && return nothing
    return rk_interaction_bounds(m, geom)
end

#==========================================================================
## Extent caching (the stand-in for Java's per-Geometry envelope cache)

Rebuild the input as a GeoInterface wrapper tree with the interaction
bounds embedded at every level, in one coordinate pass. Wrappers share the
original linework objects (a `GI.LinearRing(ring; extent)` around a
same-trait geometry stores `ring`'s coordinate backing, copying nothing),
so this costs O(#elements) small allocations plus the one extent scan the
constructor performed anyway. Levels whose stored extent is usable as
interaction bounds (`_reusable_stored_extent`) are reused as-is, so
re-wrapping an already-cached tree does no coordinate work.
==========================================================================#

_has_stored_extent(geom) =
    geom isa GI.Wrappers.WrapperGeometry && hasproperty(geom, :extent) &&
    geom.extent isa Extents.Extent

# A stored extent is reusable as interaction bounds only if it lives in the
# space the kernel prunes in: any stored extent on `Planar`, but only a 3D
# `(X, Y, Z)` extent on `Spherical` — user inputs typically carry lon/lat
# boxes, which must not be compared against unit-sphere boxes. Our own cache
# pass always stores `(X, Y, Z)`; a user storing one is trusted to mean it.
_reusable_stored_extent(::Manifold, geom) = _has_stored_extent(geom)
_reusable_stored_extent(::Spherical, geom) =
    _has_stored_extent(geom) && geom.extent isa Extents.Extent{(:X, :Y, :Z)}

_relate_cache_extents(m::Manifold, geom) = _relate_cache_extents(m, GI.trait(geom), geom)

#-- point elements: their extent is themselves, nothing to cache
_relate_cache_extents(::Manifold, ::Union{GI.AbstractPointTrait, GI.AbstractMultiPointTrait}, geom) = geom

#-- linework leaves: lines and rings (the only level where coordinates are
#-- read, and hence where the manifold's edge validation runs)
function _relate_cache_extents(m::Manifold, trait::GI.AbstractCurveTrait, line)
    (GI.isempty(line) || _reusable_stored_extent(m, line)) && return line
    _validate_relate_edges(m, line)
    return GI.geointerface_geomtype(trait)(line;
        extent = rk_interaction_bounds(m, line), crs = GI.crs(line))
end

function _relate_cache_extents(m::Manifold, trait::GI.AbstractPolygonTrait, poly)
    GI.isempty(poly) && return poly
    if _reusable_stored_extent(m, poly) && all(r -> GI.isempty(r) || _reusable_stored_extent(m, r), GI.getring(poly))
        return poly
    end
    rings = [_relate_cache_extents(m, GI.trait(r), r) for r in GI.getring(poly)]
    ext = _union_stored_extents(m, rings)
    ext === nothing && return poly
    return GI.geointerface_geomtype(trait)(rings; extent = ext, crs = GI.crs(poly))
end

#-- collections (covers Multi* types too): recurse, union the child extents
function _relate_cache_extents(m::Manifold, trait::GI.AbstractGeometryCollectionTrait, geom)
    children = [_relate_cache_extents(m, GI.trait(g), g) for g in GI.getgeom(geom)]
    ext = _union_stored_extents(m, children)
    ext === nothing && return geom
    return GI.geointerface_geomtype(trait)(children; extent = ext, crs = GI.crs(geom))
end

#-- any other trait: leave untouched
_relate_cache_extents(::Manifold, ::GI.AbstractTrait, geom) = geom

# Union of the children's extents, reading stored ones and computing only
# for non-empty children that have none (e.g. point members of a GC);
# `nothing` when no child contributes one. Computed extents go through
# `rk_interaction_bounds` so they stay in the manifold's coordinate space.
function _union_stored_extents(m::Manifold, children)
    ext = nothing
    for c in children
        ce = if _reusable_stored_extent(m, c)
            c.extent
        elseif GI.isempty(c)
            nothing
        else
            rk_interaction_bounds(m, c)
        end
        ce === nothing && continue
        ext = ext === nothing ? ce : Extents.union(ext, ce)
    end
    return ext
end

# Equivalent of Java `Geometry.getDimension()`: the inherent dimension of the
# geometry type, including empty elements (an empty polygon still has
# dimension 2); collections report the maximum over their elements
# (`DIM_FALSE` when there are none).
_geom_dimension(geom) = _geom_dimension(GI.trait(geom), geom)
_geom_dimension(::Union{GI.AbstractPointTrait, GI.AbstractMultiPointTrait}, geom) = DIM_P
_geom_dimension(::Union{GI.AbstractCurveTrait, GI.AbstractMultiCurveTrait}, geom) = DIM_L
_geom_dimension(::Union{GI.AbstractPolygonTrait, GI.AbstractMultiPolygonTrait}, geom) = DIM_A
function _geom_dimension(::GI.AbstractGeometryCollectionTrait, geom)
    dim = DIM_FALSE
    for g in GI.getgeom(geom)
        d = _geom_dimension(g)
        d > dim && (dim = d)
    end
    return dim
end
_geom_dimension(::GI.AbstractTrait, geom) = DIM_FALSE

# Port of RelateGeometry.analyzeDimensions, returning
# `(dim, has_points, has_lines, has_areas)`. The Java `instanceof
# Point/LineString/Polygon` checks are widened to the corresponding GI
# abstract traits (e.g. `GI.AbstractPolygonTrait` also covers triangles
# etc.), per the Task 12 review.
function _analyze_dimensions(geom, dim0::Int8, is_geom_empty::Bool)
    is_geom_empty && return (dim0, false, false, false)
    return _analyze_dimensions(GI.trait(geom), geom, dim0)
end
_analyze_dimensions(::Union{GI.AbstractPointTrait, GI.AbstractMultiPointTrait}, geom, dim0) =
    (DIM_P, true, false, false)
_analyze_dimensions(::Union{GI.AbstractCurveTrait, GI.AbstractMultiCurveTrait}, geom, dim0) =
    (DIM_L, false, true, false)
_analyze_dimensions(::Union{GI.AbstractPolygonTrait, GI.AbstractMultiPolygonTrait}, geom, dim0) =
    (DIM_A, false, false, true)
#-- analyze a (possibly mixed type) collection
_analyze_dimensions(::GI.AbstractTrait, geom, dim0) =
    _analyze_collection_dimensions(geom, dim0, false, false, false)

# The recursive element walk of analyzeDimensions (Java uses a
# GeometryCollectionIterator; only atomic elements match the checks).
function _analyze_collection_dimensions(geom, dim, has_points, has_lines, has_areas)
    for g in GI.getgeom(geom)
        dim, has_points, has_lines, has_areas = _analyze_element_dimensions(
            GI.trait(g), g, dim, has_points, has_lines, has_areas)
    end
    return (dim, has_points, has_lines, has_areas)
end
_analyze_element_dimensions(::GI.AbstractGeometryCollectionTrait, g, dim, hp, hl, ha) =
    _analyze_collection_dimensions(g, dim, hp, hl, ha)
function _analyze_element_dimensions(::GI.AbstractPointTrait, g, dim, hp, hl, ha)
    GI.isempty(g) && return (dim, hp, hl, ha)
    return (max(dim, DIM_P), true, hl, ha)
end
function _analyze_element_dimensions(::GI.AbstractCurveTrait, g, dim, hp, hl, ha)
    GI.isempty(g) && return (dim, hp, hl, ha)
    return (max(dim, DIM_L), hp, true, ha)
end
function _analyze_element_dimensions(::GI.AbstractPolygonTrait, g, dim, hp, hl, ha)
    GI.isempty(g) && return (dim, hp, hl, ha)
    return (max(dim, DIM_A), hp, hl, true)
end
_analyze_element_dimensions(::GI.AbstractTrait, g, dim, hp, hl, ha) = (dim, hp, hl, ha)

# Port of RelateGeometry.isZeroLengthLine.
function _is_zero_length_line(geom, dim::Int8)
    #-- avoid expensive zero-length calculation if not linear
    dim == DIM_L || return false
    return _is_zero_length(geom)
end

# Port of RelateGeometry.isZeroLength(Geometry): tests if all linear elements
# are zero-length. For efficiency the test avoids computing actual length.
_is_zero_length(geom) = _is_zero_length(GI.trait(geom), geom)
function _is_zero_length(::GI.AbstractGeometryCollectionTrait, geom)
    for g in GI.getgeom(geom)
        _is_zero_length(g) || return false
    end
    return true
end
_is_zero_length(::GI.AbstractCurveTrait, geom) = _is_zero_length_linestring(geom)
_is_zero_length(::GI.AbstractTrait, geom) = true

# Port of RelateGeometry.isZeroLength(LineString): exact coordinate equality
# of every point with the first one.
function _is_zero_length_linestring(line)
    n = GI.npoint(line)
    if n >= 2
        p0 = _tuple_point(GI.getpoint(line, 1))
        for i in 1:n
            p = _tuple_point(GI.getpoint(line, i))
            #-- most non-zero-len lines will trigger this right away
            _equals2(p0, p) || return false
        end
    end
    return true
end

# Java getGeometry / isPrepared / getEnvelope / getDimension.
get_geometry(rg::RelateGeometry) = rg.geom
is_prepared(rg::RelateGeometry) = rg.is_prepared
get_extent(rg::RelateGeometry) = rg.extent
get_dimension(rg::RelateGeometry) = rg.dim

# Port of RelateGeometry.hasDimension(dim).
function has_dimension(rg::RelateGeometry, dim::Integer)
    dim == DIM_P && return rg.has_points
    dim == DIM_L && return rg.has_lines
    dim == DIM_A && return rg.has_areas
    return false
end

has_area_and_line(rg::RelateGeometry) = rg.has_areas && rg.has_lines

"""
    get_dimension_real(rg::RelateGeometry)

Gets the actual non-empty dimension of the geometry.
Zero-length LineStrings are treated as Points.
"""
function get_dimension_real(rg::RelateGeometry)
    rg.is_geom_empty && return DIM_FALSE
    get_dimension(rg) == DIM_L && rg.is_line_zero_len && return DIM_P
    rg.has_areas && return DIM_A
    rg.has_lines && return DIM_L
    return DIM_P
end

has_edges(rg::RelateGeometry) = rg.has_lines || rg.has_areas

# Port of RelateGeometry.getLocator (lazy).
function _get_locator(rg::RelateGeometry)
    loc = rg.locator
    loc === nothing || return loc
    loc = RelatePointLocator(rg.m, rg.geom; exact = rg.exact,
        is_prepared = rg.is_prepared, boundary_rule = rg.boundary_rule)
    rg.locator = loc
    return loc
end

# Port of RelateGeometry.isNodeInArea.
function is_node_in_area(rg::RelateGeometry, node_pt, parent_polygonal)
    loc = locate_node_with_dim(_get_locator(rg), node_pt, parent_polygonal)
    return loc == DL_AREA_INTERIOR
end

locate_line_end_with_dim(rg::RelateGeometry, p) =
    locate_line_end_with_dim(_get_locator(rg), p)

"""
    locate_area_vertex(rg::RelateGeometry, pt)

Locates a vertex of a polygon. A vertex of a Polygon or MultiPolygon is on
the boundary; but a vertex of an overlapped polygon in a GeometryCollection
may be in the interior.
"""
function locate_area_vertex(rg::RelateGeometry, pt)
    #=
    Can pass a `nothing` polygon, because the point is an exact vertex,
    which will be detected as being on the boundary of its polygon
    =#
    return locate_node(rg, pt, nothing)
end

locate_node(rg::RelateGeometry, pt, parent_polygonal) =
    locate_node(_get_locator(rg), pt, parent_polygonal)

locate_with_dim(rg::RelateGeometry, pt) = locate_with_dim(_get_locator(rg), pt)

"""
    is_self_noding_required(rg::RelateGeometry)

Indicates whether the geometry requires self-noding for correct evaluation
of specific spatial predicates. Self-noding is required for geometries which
may self-cross — i.e. lines, and overlapping elements in
GeometryCollections. Self-noding is not required for polygonal geometries,
since they can only touch at vertices.
"""
function is_self_noding_required(rg::RelateGeometry)
    trait = GI.trait(rg.geom)
    if trait isa Union{GI.AbstractPointTrait, GI.AbstractMultiPointTrait,
            GI.AbstractPolygonTrait, GI.AbstractMultiPolygonTrait}
        return false
    end
    #-- a GC with a single polygon does not need noding
    rg.has_areas && _num_geometries(rg.geom) == 1 && return false
    #-- GCs with only points do not need noding
    !rg.has_areas && !rg.has_lines && return false
    return true
end

# Java Geometry.getNumGeometries: element count for collections, 1 otherwise.
_num_geometries(geom) =
    GI.trait(geom) isa GI.AbstractGeometryCollectionTrait ? GI.ngeom(geom) : 1

"""
    is_polygonal(rg::RelateGeometry)

Tests whether the geometry has polygonal topology. This is not the case if
it is a GeometryCollection containing more than one polygon (since they may
overlap or be adjacent). The significance is that polygonal topology allows
more assumptions about the location of boundary vertices.
"""
function is_polygonal(rg::RelateGeometry)
    #TODO: also true for a GC containing one polygonal element (and possibly some lower-dimension elements)
    return GI.trait(rg.geom) isa Union{GI.AbstractPolygonTrait, GI.AbstractMultiPolygonTrait}
end

is_geom_empty(rg::RelateGeometry) = rg.is_geom_empty

has_boundary(rg::RelateGeometry) = has_boundary(_get_locator(rg))

function get_unique_points(rg::RelateGeometry)
    #-- will be re-used in prepared mode
    up = rg.unique_points
    up === nothing || return up
    up = _create_unique_points(rg.m, rg.geom)
    rg.unique_points = up
    return up
end

# Port of RelateGeometry.createUniquePoints. Only called on P geometries.
# (Java uses ComponentCoordinateExtracter, which records the first coordinate
# of each point/line component; for point geometries that is every point.)
function _create_unique_points(m::Manifold, geom)
    set = Set{_kernel_point_type(m)}()
    _add_component_coordinates!(set, m, geom)
    return set
end

_add_component_coordinates!(set, m, geom) =
    _add_component_coordinates!(set, m, GI.trait(geom), geom)
function _add_component_coordinates!(set, m, ::GI.AbstractGeometryCollectionTrait, geom)
    for g in GI.getgeom(geom)
        _add_component_coordinates!(set, m, g)
    end
    return nothing
end
function _add_component_coordinates!(set, m, ::GI.AbstractPointTrait, geom)
    GI.isempty(geom) && return nothing
    push!(set, _to_kernel_point(m, geom))
    return nothing
end
function _add_component_coordinates!(set, m, ::GI.AbstractCurveTrait, geom)
    GI.isempty(geom) && return nothing
    push!(set, _to_kernel_point(m, GI.getpoint(geom, 1)))
    return nothing
end
_add_component_coordinates!(set, m, ::GI.AbstractTrait, geom) = nothing

# Port of RelateGeometry.getEffectivePoints: the point elements which are not
# covered by an element of higher dimension. (This JTS version has no
# MAX_EFFECTIVE_POINTS cap; all points are checked.) Returns the point
# geometries themselves, as in Java.
function get_effective_points(rg::RelateGeometry)
    pt_list_all = Any[]
    _extract_point_elements!(pt_list_all, rg.geom)

    get_dimension_real(rg) <= DIM_P && return pt_list_all

    #-- only return Points not covered by another element
    pt_list = Any[]
    for p in pt_list_all
        GI.isempty(p) && continue
        loc_dim = locate_with_dim(rg, _to_kernel_point(rg.m, p))
        if dimloc_dimension(loc_dim) == DIM_P
            push!(pt_list, p)
        end
    end
    return pt_list
end

# Equivalent of Java PointExtracter.getPoints: every Point element, including
# those nested in collections.
_extract_point_elements!(list, geom) =
    _extract_point_elements!(list, GI.trait(geom), geom)
function _extract_point_elements!(list, ::GI.AbstractPointTrait, geom)
    push!(list, geom)
    return nothing
end
function _extract_point_elements!(list, ::GI.AbstractGeometryCollectionTrait, geom)
    for g in GI.getgeom(geom)
        _extract_point_elements!(list, g)
    end
    return nothing
end
_extract_point_elements!(list, ::GI.AbstractTrait, geom) = nothing

"""
    extract_segment_strings(rg::RelateGeometry, is_a::Bool, ext_filter)

Extract [`RelateSegmentString`](@ref)s from the geometry which intersect a
given extent (one per line, one per polygon ring). If `ext_filter` is
`nothing` all edges are extracted.

!!! warning
    `nothing` here means *no filter* (Java's prepared-mode `null`), while
    `get_extent(rg)` returns `nothing` for an *empty* geometry, where JTS's
    null Envelope intersects nothing. Never forward an empty geometry's
    extent as the filter — callers (the engine's `computeAtEdges` port)
    must early-return on empty inputs before extraction.
"""
function extract_segment_strings(rg::RelateGeometry, is_a::Bool, ext_filter)
    #-- `RelateSegmentString{P}` is concrete (its geometry references are
    #-- opaque), so the vector is concretely typed and the per-segment loops
    #-- downstream (`_segment_extent_table`, the NestedLoop enumerator) stay
    #-- statically dispatched, for any input geometry type
    seg_strings = Vector{RelateSegmentString{_kernel_point_type(rg.m)}}()
    _extract_segment_strings!(rg, is_a, ext_filter, rg.geom, seg_strings)
    return seg_strings
end

function _extract_segment_strings!(rg::RelateGeometry, is_a::Bool, ext_filter, geom, seg_strings)
    trait = GI.trait(geom)
    #-- record if parent is MultiPolygon
    parent_polygonal = trait isa GI.AbstractMultiPolygonTrait ? geom : nothing

    # Java iterates getGeometryN over getNumGeometries: for an atomic
    # geometry that yields the geometry itself.
    elements = trait isa GI.AbstractGeometryCollectionTrait ? GI.getgeom(geom) : (geom,)
    for g in elements
        # Java `instanceof GeometryCollection` covers the Multi* types too.
        if GI.trait(g) isa GI.AbstractGeometryCollectionTrait
            _extract_segment_strings!(rg, is_a, ext_filter, g, seg_strings)
        else
            #-- an atomic input geometry's extent is already cached on `rg`
            #-- (Java's getEnvelopeInternal cache); don't rescan it below
            elem_ext = g === rg.geom ? get_extent(rg) : missing
            _extract_segment_strings_from_atomic!(rg, is_a, g, parent_polygonal,
                ext_filter, seg_strings, elem_ext)
        end
    end
    return nothing
end

function _extract_segment_strings_from_atomic!(rg::RelateGeometry, is_a::Bool, geom,
        parent_polygonal, ext_filter, seg_strings, elem_ext = missing)
    GI.isempty(geom) && return nothing
    if ext_filter !== nothing
        if elem_ext === missing
            elem_ext = rk_interaction_bounds(rg.m, geom)
        end
        !Extents.intersects(ext_filter, elem_ext) && return nothing
    end

    rg.element_id += Int32(1)
    trait = GI.trait(geom)
    if trait isa GI.AbstractCurveTrait
        pts = _to_kernel_points(rg.m, geom)
        ss = _rss_create_line(pts, is_a, rg.element_id, rg)
        push!(seg_strings, ss)
    elseif trait isa GI.AbstractPolygonTrait
        parent_poly = parent_polygonal !== nothing ? parent_polygonal : geom
        #-- the exterior ring's extent is the element extent (for an invalid
        #-- polygon with a hole outside its shell it is a superset, which can
        #-- only under-prune — extracted non-interacting edges are harmless)
        _extract_ring_to_segment_string!(rg, is_a, GI.getexterior(geom), 0, ext_filter,
            parent_poly, seg_strings, elem_ext)
        for (i, hole) in enumerate(GI.gethole(geom))
            _extract_ring_to_segment_string!(rg, is_a, hole, i, ext_filter, parent_poly, seg_strings)
        end
    end
    return nothing
end

function _extract_ring_to_segment_string!(rg::RelateGeometry, is_a::Bool, ring, ring_id::Integer,
        ext_filter, parent_poly, seg_strings, ring_ext = missing)
    GI.isempty(ring) && return nothing
    if ext_filter !== nothing
        if ring_ext === missing
            ring_ext = rk_interaction_bounds(rg.m, ring)
        end
        !Extents.intersects(ext_filter, ring_ext) && return nothing
    end

    #-- orient the points if required
    require_cw = ring_id == 0
    pts = _to_kernel_points(rg.m, ring)
    pts = _orient_ring(rg.m, pts, require_cw, ring_id != 0; exact = rg.exact)
    ss = _rss_create_ring(pts, is_a, rg.element_id, ring_id, parent_poly, rg)
    push!(seg_strings, ss)
    return nothing
end

# Port of RelateGeometry.orient (static; moved here from point_locator.jl in
# Task 13 — `AdjacentEdgeLocator._add_ring!` also uses it): coordinate vector
# of `pts` with the ring's denoted region on the requested side (`orient_cw =
# true` ⇒ on the right), reversing a copy only if needed. Which side the
# region lies on in the stored order comes from `_ring_interior_on_left`.
function _orient_ring(m, pts::Vector, orient_cw::Bool, is_hole::Bool; exact)
    is_flipped = orient_cw == _ring_interior_on_left(m, pts, is_hole; exact)
    return is_flipped ? reverse(pts) : pts
end

#=
The one per-ring bit every consumer shares: whether the region the ring
DENOTES — the shell region for a shell or bare ring, the cavity for a hole
— lies on the left of the stored vertex order. On the plane, and on an
unoriented spherical manifold, the denoted region is the enclosed one and
the bit is the ring's winding (`_ring_is_ccw`); `Spherical(; oriented =
true)` overrides this (kernel_spherical.jl) with the declared role — there
the stored winding is authoritative.
=#
_ring_interior_on_left(m, pts::Vector, is_hole::Bool; exact) =
    _ring_is_ccw(m, pts; exact)

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

#==========================================================================
## RelateSegmentString (port of JTS RelateSegmentString.java)
==========================================================================#

"""
    RelateSegmentString{P}

Models a linear edge of a [`RelateGeometry`](@ref): the coordinate vector of
one line or one polygon ring, tagged with which input geometry it came from
(`is_a`), its dimension, the element/ring ids assigned during extraction,
and (for rings) the parent polygonal geometry.

In JTS this extends `BasicSegmentString`; here the coordinates are stored
directly in `pts`. Segment indices are 1-based: segment `i` runs from
`pts[i]` to `pts[i + 1]` (the Java equivalents are 0-based).

The geometry references are deliberately opaque (abstract field types, as
they are in Java): `parent_polygonal` is only ever compared by identity and
carried into [`NodeSection`](@ref)s, so parameterizing on the input geometry
type would only re-specialize the whole edge machinery per geometry-type
pair (a pure compile-time cost). Only `pts` — the per-segment hot path —
stays concretely typed, on the manifold's kernel point type `P`.
"""
struct RelateSegmentString{P}
    is_a::Bool
    dim::Int8
    id::Int32
    ring_id::Int32
    input_geom::RelateGeometry
    parent_polygonal::Any
    pts::Vector{P}
end

# Port of RelateSegmentString.createLine.
_rss_create_line(pts::Vector, is_a::Bool, element_id::Integer, parent::RelateGeometry) =
    _rss_create(pts, is_a, DIM_L, element_id, -1, nothing, parent)

# Port of RelateSegmentString.createRing.
_rss_create_ring(pts::Vector, is_a::Bool, element_id::Integer, ring_id::Integer,
        poly, parent::RelateGeometry) =
    _rss_create(pts, is_a, DIM_A, element_id, ring_id, poly, parent)

# Port of RelateSegmentString.createSegmentString.
function _rss_create(pts::Vector, is_a::Bool, dim::Int8, element_id::Integer,
        ring_id::Integer, poly, parent::RelateGeometry)
    pts = _remove_repeated_points(pts)
    return RelateSegmentString(is_a, dim, Int32(element_id), Int32(ring_id), parent, poly, pts)
end

# Port of CoordinateArrays.removeRepeatedPoints (via hasRepeatedPoints):
# drops consecutive coordinates that are exactly equal, returning the input
# vector unchanged (not copied) when there is nothing to remove.
function _remove_repeated_points(pts::Vector)
    any(i -> _equals2(pts[i], pts[i + 1]), 1:(length(pts) - 1)) || return pts
    out = [pts[1]]
    for i in 2:length(pts)
        _equals2(pts[i], last(out)) || push!(out, pts[i])
    end
    return out
end

get_geometry(ss::RelateSegmentString) = ss.input_geom
get_polygonal(ss::RelateSegmentString) = ss.parent_polygonal

# Port of BasicSegmentString.isClosed.
is_closed(ss::RelateSegmentString) = _equals2(ss.pts[1], ss.pts[end])

"""
    create_node_section(ss::RelateSegmentString, seg_index::Integer, node::NodeKey)

The [`NodeSection`](@ref) of this segment string at the node `node`, known
to lie on segment `seg_index`.

Port of RelateSegmentString.createNodeSection, with the symbolic twist
(design D2): the Java method takes the intersection `Coordinate`; here the
node is identified by its [`NodeKey`](@ref). For a vertex node (`SS_TOUCH`
or other vertex incidences) the key carries the exact coordinate and the
incident vertices are found as in Java. A proper-crossing node (`SS_PROPER`)
lies strictly inside the segment, so the incident vertices are the segment
endpoints and the node is never at a vertex.
"""
function create_node_section(ss::RelateSegmentString, seg_index::Integer, node::NodeKey)
    if node.is_crossing
        #-- a proper crossing is interior to its segment
        is_node_at_vertex = false
        prev = ss.pts[seg_index]
        next = ss.pts[seg_index + 1]
    else
        pt = node.pt
        is_node_at_vertex =
            _equals2(pt, ss.pts[seg_index]) || _equals2(pt, ss.pts[seg_index + 1])
        prev = prev_vertex(ss, seg_index, pt)
        next = next_vertex(ss, seg_index, pt)
    end
    return NodeSection(ss.is_a, ss.dim, ss.id, ss.ring_id, ss.parent_polygonal,
        is_node_at_vertex, prev, node, next)
end

# Port of RelateSegmentString.prevVertex: the vertex before the node lying
# on segment `seg_index`, or `nothing` if none exists.
function prev_vertex(ss::RelateSegmentString, seg_index::Integer, pt)
    seg_start = ss.pts[seg_index]
    _equals2(seg_start, pt) || return seg_start
    #-- pt is at segment start, so get previous vertex
    seg_index > 1 && return ss.pts[seg_index - 1]
    is_closed(ss) && return _prev_in_ring(ss, seg_index)
    return nothing
end

# Port of RelateSegmentString.nextVertex: the vertex after the node lying
# on segment `seg_index`, or `nothing` if none exists.
function next_vertex(ss::RelateSegmentString, seg_index::Integer, pt)
    seg_end = ss.pts[seg_index + 1]
    _equals2(seg_end, pt) || return seg_end
    #-- pt is at seg end, so get next vertex
    seg_index < length(ss.pts) - 1 && return ss.pts[seg_index + 2]
    is_closed(ss) && return _next_in_ring(ss, seg_index + 1)
    #-- segstring is not closed, so there is no next segment
    return nothing
end

# Ports of BasicSegmentString.prevInRing / nextInRing (1-based; the closing
# point duplicates the first, so the wraparound skips it).
function _prev_in_ring(ss::RelateSegmentString, index::Integer)
    prev_index = index - 1
    prev_index < 1 && (prev_index = length(ss.pts) - 1)
    return ss.pts[prev_index]
end

function _next_in_ring(ss::RelateSegmentString, index::Integer)
    next_index = index + 1
    next_index > length(ss.pts) && (next_index = 2)
    return ss.pts[next_index]
end

"""
    is_containing_segment(ss::RelateSegmentString, seg_index::Integer, pt)

Tests if a segment intersection point has that segment as its canonical
containing segment. Segments are half-closed, and contain their start point
but not the endpoint, except for the final segment in a non-closed segment
string, which contains its endpoint as well. This test ensures that vertices
are assigned to a unique segment in a segment string. In particular, this
avoids double-counting intersections which lie exactly at segment endpoints.
"""
function is_containing_segment(ss::RelateSegmentString, seg_index::Integer, pt)
    #-- intersection is at segment start vertex - process it
    _equals2(pt, ss.pts[seg_index]) && return true
    if _equals2(pt, ss.pts[seg_index + 1])
        is_final_segment = seg_index == length(ss.pts) - 1
        (is_closed(ss) || !is_final_segment) && return false
        #-- for final segment, process intersections with final endpoint
        return true
    end
    #-- intersection is interior - process it
    return true
end
