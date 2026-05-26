# # RelateNG geometry wrapper

"""
    RelateGeometry(geom; prepared = false, boundary_node_rule = Mod2BoundaryNodeRule())

Prepared-style wrapper that caches the geometry metadata needed by RelateNG.
"""
mutable struct RelateGeometry{G,R <: BoundaryNodeRule,E}
    geom::G
    prepared::Bool
    boundary_node_rule::R
    extent::E
    dimension::TopologicalDimension
    has_points::Bool
    has_lines::Bool
    has_areas::Bool
    all_linework_zero_length::Bool
    is_empty::Bool
    unique_points_cache::Union{Nothing,Set{Tuple{Float64,Float64}}}
    locator_cache::Union{Nothing,RelatePointLocator{G,R}}
    segment_strings_cache::Dict{Any,Any}
    prepared_edge_index_cache::Dict{DataType,Any}
end

function RelateGeometry(
    geom;
    prepared::Bool = false,
    boundary_node_rule::BoundaryNodeRule = Mod2BoundaryNodeRule(),
)
    dimension, has_points, has_lines, has_areas = _relate_analyze_dimensions(geom)
    is_empty = _relate_is_empty(geom)
    all_linework_zero_length = dimension == dim_line && _relate_all_linework_zero_length(geom)
    extent = _relate_cached_extent(geom)
    return RelateGeometry(
        geom,
        prepared,
        boundary_node_rule,
        extent,
        dimension,
        has_points,
        has_lines,
        has_areas,
        all_linework_zero_length,
        is_empty,
        nothing,
        nothing,
        Dict{Any,Any}(),
        Dict{DataType,Any}(),
    )
end

RelateGeometry(alg::RelateNG, geom) =
    RelateGeometry(geom; prepared = alg.prepared, boundary_node_rule = alg.boundary_node_rule)

function _relate_is_empty(geom)
    return _relate_is_empty(GI.trait(geom), geom)
end

_relate_is_empty(::GI.FeatureTrait, feature) =
    _relate_is_empty(GI.geometry(feature))

_relate_is_empty(::GI.FeatureCollectionTrait, fc) =
    all(_relate_is_empty, GI.getfeature(fc))

_relate_is_empty(::Nothing, iterable) =
    all(_relate_is_empty, iterable)

_relate_is_empty(::GI.AbstractGeometryTrait, geom) = GI.isempty(geom)

function _relate_cached_extent(geom)
    try
        return GI.extent(geom)
    catch
        return nothing
    end
end

function _relate_analyze_dimensions(geom)
    return _relate_analyze_dimensions(GI.trait(geom), geom)
end

_relate_analyze_dimensions(::GI.FeatureTrait, feature) =
    _relate_analyze_dimensions(GI.geometry(feature))

function _relate_analyze_dimensions(::GI.FeatureCollectionTrait, fc)
    return _relate_analyze_dimension_iterable(GI.getfeature(fc))
end

function _relate_analyze_dimensions(::Nothing, iterable)
    return _relate_analyze_dimension_iterable(iterable)
end

function _relate_analyze_dimension_iterable(iterable)
    dimension = dim_false
    has_points = false
    has_lines = false
    has_areas = false
    for geom in iterable
        child_dimension, child_points, child_lines, child_areas = _relate_analyze_dimensions(geom)
        dimension = max_dimension(dimension, child_dimension)
        has_points |= child_points
        has_lines |= child_lines
        has_areas |= child_areas
    end
    return dimension, has_points, has_lines, has_areas
end

_relate_analyze_dimensions(::GI.PointTrait, geom) =
    GI.isempty(geom) ? (dim_false, false, false, false) : (dim_point, true, false, false)

_relate_analyze_dimensions(::GI.MultiPointTrait, geom) =
    GI.isempty(geom) ? (dim_false, false, false, false) : (dim_point, true, false, false)

_relate_analyze_dimensions(::GI.AbstractCurveTrait, geom) =
    GI.isempty(geom) ? (dim_false, false, false, false) : (dim_line, false, true, false)

_relate_analyze_dimensions(::GI.PolygonTrait, geom) =
    GI.isempty(geom) ? (dim_false, false, false, false) : (dim_area, false, false, true)

_relate_analyze_dimensions(::GI.MultiPolygonTrait, geom) =
    GI.isempty(geom) ? (dim_false, false, false, false) : (dim_area, false, false, true)

function _relate_analyze_dimensions(::GI.AbstractGeometryTrait, geom)
    GI.isempty(geom) && return (dim_false, false, false, false)
    return _relate_analyze_dimension_iterable(GI.getgeom(geom))
end

function _relate_all_linework_zero_length(geom)
    return _relate_all_linework_zero_length(GI.trait(geom), geom)
end

_relate_all_linework_zero_length(::GI.FeatureTrait, feature) =
    _relate_all_linework_zero_length(GI.geometry(feature))

_relate_all_linework_zero_length(::GI.FeatureCollectionTrait, fc) =
    all(_relate_all_linework_zero_length, GI.getfeature(fc))

_relate_all_linework_zero_length(::Nothing, iterable) =
    all(_relate_all_linework_zero_length, iterable)

_relate_all_linework_zero_length(::GI.PointTrait, geom) = true
_relate_all_linework_zero_length(::GI.MultiPointTrait, geom) = true
_relate_all_linework_zero_length(::GI.PolygonTrait, geom) = true
_relate_all_linework_zero_length(::GI.MultiPolygonTrait, geom) = true

function _relate_all_linework_zero_length(::GI.AbstractCurveTrait, curve)
    GI.isempty(curve) && return true
    GI.npoint(curve) < 2 && return true

    first_point = _tuple_point(GI.getpoint(curve, 1))
    for point in GI.getpoint(curve)
        _tuple_point(point) == first_point || return false
    end
    return true
end

function _relate_all_linework_zero_length(::GI.AbstractGeometryTrait, geom)
    GI.isempty(geom) && return true
    return all(_relate_all_linework_zero_length, GI.getgeom(geom))
end

"""
    relate_dimension_real(relate_geometry)

Return the non-empty dimension, treating all-zero-length linework as points.
"""
function relate_dimension_real(relate_geometry::RelateGeometry)
    relate_geometry.is_empty && return dim_false
    if relate_geometry.dimension == dim_line && relate_geometry.all_linework_zero_length
        return dim_point
    end
    relate_geometry.has_areas && return dim_area
    relate_geometry.has_lines && return dim_line
    relate_geometry.has_points && return dim_point
    return dim_false
end

relate_has_dimension(relate_geometry::RelateGeometry, dimension::TopologicalDimension) =
    dimension == dim_point ? relate_geometry.has_points :
    dimension == dim_line ? relate_geometry.has_lines :
    dimension == dim_area ? relate_geometry.has_areas :
    dimension == dim_false ? relate_geometry.is_empty :
    false

relate_has_edges(relate_geometry::RelateGeometry) =
    relate_geometry.has_lines || relate_geometry.has_areas

relate_has_area_and_line(relate_geometry::RelateGeometry) =
    relate_geometry.has_areas && relate_geometry.has_lines

"""
    relate_has_boundary(relate_geometry)

Return whether linear components have boundary endpoints under the node rule.
"""
function relate_has_boundary(relate_geometry::RelateGeometry)
    relate_geometry.has_lines || return false
    endpoint_counts = Dict{Tuple{Float64,Float64},Int}()
    _relate_collect_line_endpoints!(
        endpoint_counts,
        GI.trait(relate_geometry.geom),
        relate_geometry.geom,
    )
    return any(count -> is_in_boundary(relate_geometry.boundary_node_rule, count), values(endpoint_counts))
end

relate_is_polygonal(relate_geometry::RelateGeometry) =
    _relate_is_polygonal(GI.trait(relate_geometry.geom), relate_geometry.geom)

function relate_is_self_noding_required(relate_geometry::RelateGeometry)
    relate_is_polygonal(relate_geometry) && return false
    relate_geometry.dimension == dim_point && return false
    if relate_geometry.has_areas && _relate_direct_child_count(relate_geometry.geom) == 1
        return false
    end
    (!relate_geometry.has_areas && !relate_geometry.has_lines) && return false
    return true
end

function _relate_direct_child_count(geom)
    trait = GI.trait(geom)
    trait isa GI.AbstractGeometryTrait || return 1
    return GI.ngeom(geom)
end

_relate_collect_line_endpoints!(counts, ::GI.PointTrait, geom) = counts
_relate_collect_line_endpoints!(counts, ::GI.MultiPointTrait, geom) = counts
_relate_collect_line_endpoints!(counts, ::GI.PolygonTrait, geom) = counts
_relate_collect_line_endpoints!(counts, ::GI.MultiPolygonTrait, geom) = counts

function _relate_collect_line_endpoints!(counts, ::GI.AbstractCurveTrait, curve)
    GI.isempty(curve) && return counts
    GI.npoint(curve) == 0 && return counts
    first_point = _tuple_point(GI.getpoint(curve, 1), Float64)
    last_point = _tuple_point(GI.getpoint(curve, GI.npoint(curve)), Float64)
    counts[first_point] = get(counts, first_point, 0) + 1
    counts[last_point] = get(counts, last_point, 0) + 1
    return counts
end

function _relate_collect_line_endpoints!(counts, ::GI.AbstractGeometryTrait, geom)
    GI.isempty(geom) && return counts
    for child in GI.getgeom(geom)
        _relate_collect_line_endpoints!(counts, GI.trait(child), child)
    end
    return counts
end

"""
    relate_point_locator(relate_geometry)

Return the cached `RelatePointLocator` for a wrapped geometry.
"""
function relate_point_locator(relate_geometry::RelateGeometry)
    if isnothing(relate_geometry.locator_cache)
        relate_geometry.locator_cache = RelatePointLocator(
            relate_geometry.geom;
            boundary_node_rule = relate_geometry.boundary_node_rule,
            prepared = relate_geometry.prepared,
        )
    end
    return relate_geometry.locator_cache
end

relate_locate(relate_geometry::RelateGeometry, point; kwargs...) =
    relate_locate(relate_point_locator(relate_geometry), point; kwargs...)

relate_locate_with_dim(relate_geometry::RelateGeometry, point; kwargs...) =
    relate_locate_with_dim(relate_point_locator(relate_geometry), point; kwargs...)

relate_locate_node_with_dim(relate_geometry::RelateGeometry, point, parent_polygonal = nothing; kwargs...) =
    relate_locate_node_with_dim(relate_point_locator(relate_geometry), point, parent_polygonal; kwargs...)

relate_locate_line_end_with_dim(relate_geometry::RelateGeometry, point; kwargs...) =
    relate_locate_line_end_with_dim(relate_point_locator(relate_geometry), point; kwargs...)

relate_is_node_in_area(relate_geometry::RelateGeometry, point, parent_polygonal = nothing; kwargs...) =
    relate_locate_node_with_dim(relate_geometry, point, parent_polygonal; kwargs...) ==
    DimensionLocation(dim_area, loc_interior)

relate_locate_area_vertex(relate_geometry::RelateGeometry, point; kwargs...) =
    relate_locate_node_with_dim(relate_geometry, point, nothing; kwargs...)

"""
    relate_unique_points(relate_geometry)

Return cached unique point coordinates from point components.
"""
function relate_unique_points(relate_geometry::RelateGeometry)
    if isnothing(relate_geometry.unique_points_cache)
        relate_geometry.unique_points_cache =
            Set(extracted.point for extracted in extract_ng_points(relate_geometry.geom))
    end
    return relate_geometry.unique_points_cache
end

"""
    relate_effective_points(relate_geometry)

Return point components not covered by higher-dimensional collection elements.
"""
function relate_effective_points(relate_geometry::RelateGeometry)
    points = extract_ng_points(relate_geometry.geom)
    dimension_value(relate_dimension_real(relate_geometry)) <= dimension_value(dim_point) && return points

    return filter(points) do extracted
        relate_locate_with_dim(relate_geometry, extracted.point).dimension == dim_point
    end
end

"""
    relate_segment_strings(relate_geometry; input_side = input_a, extent = nothing)

Extract and cache RelateNG-oriented segment strings for the wrapped geometry.
"""
function relate_segment_strings(
    relate_geometry::RelateGeometry,
    ::Type{T} = Float64;
    input_side::NGInputSide = input_a,
    extent = nothing,
) where {T}
    key = (T, input_side, extent)
    return get!(relate_geometry.segment_strings_cache, key) do
        extract_ng_segment_strings(
            relate_geometry.geom,
            T;
            input_side,
            extent,
            orient_rings = :relateng,
        )
    end
end

struct RelateSegmentRecord{S,E,X}
    segment::S
    segment_index::Int
    edge::E
    extent::X
end

struct RelatePreparedEdgeIndex{R,L,I}
    records::R
    lines::L
    index::I
end

"""
    relate_prepared_edge_index(relate_geometry, [T])

Return a cached NaturalIndex over A-side segment extents, or `nothing`.
"""
function relate_prepared_edge_index(relate_geometry::RelateGeometry, ::Type{T} = Float64) where {T}
    return get!(relate_geometry.prepared_edge_index_cache, T) do
        segments = relate_segment_strings(relate_geometry, T; input_side = input_a)
        records, lines, extents = _relate_segment_records(segments, T)
        isempty(records) && return nothing
        RelatePreparedEdgeIndex(records, lines, NaturalIndexing.NaturalIndex(extents))
    end
end

function _relate_segment_records(segments, ::Type{T}) where {T}
    records = Any[]
    lines = Any[]
    extents = Extents.Extent[]
    for segment in segments
        length(segment.points) < 2 && continue
        for i in 1:(length(segment.points) - 1)
            p1, p2 = segment.points[i], segment.points[i + 1]
            p1 == p2 && continue
            edge = (p1, p2)
            extent = ng_segment_extent(edge, T)
            push!(records, RelateSegmentRecord(segment, i, edge, extent))
            push!(lines, _lineedge(edge, T))
            push!(extents, extent)
        end
    end
    return records, lines, extents
end

function _relate_segment_lines(segments, ::Type{T}) where {T}
    lines = Any[]
    for segment in segments
        length(segment.points) < 2 && continue
        for i in 1:(length(segment.points) - 1)
            p1, p2 = segment.points[i], segment.points[i + 1]
            p1 == p2 && continue
            push!(lines, _lineedge((p1, p2), T))
        end
    end
    return lines
end
