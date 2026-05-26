# # Provenance-preserving geometry extraction for JTS NG
#
# Existing `flatten`/`eachedge` helpers remain the default traversal tools.
# This layer is only for NG paths that need operand, element, and ring context.

"""
    NGInputSide

Identifies whether extracted topology comes from operand A or B.
"""
@enum NGInputSide::Int8 input_a=0 input_b=1

"""
    NGRingRole

Role of a polygon ring in its parent polygonal component.
"""
@enum NGRingRole::Int8 ring_none=0 ring_shell=1 ring_hole=2

"""
    NGRingOrientation

Planar orientation of polygon ring coordinates before NG reorientation.
"""
@enum NGRingOrientation::Int8 ring_orientation_none=0 ring_clockwise=1 ring_counterclockwise=2

const NG_NO_RING_ID = -1

"""
    NGSegmentSource

Source metadata attached to an extracted line or area segment string.
"""
struct NGSegmentSource{G,P}
    input_side::NGInputSide
    source_dimension::TopologicalDimension
    element_id::Int
    ring_id::Int
    ring_role::NGRingRole
    source_orientation::NGRingOrientation
    depth_delta::Int8
    coordinates_reversed::Bool
    geometry::G
    parent_polygonal::P
end

"""
    NGSegmentString

Coordinate sequence plus provenance for later RelateNG or OverlayNG noding.
"""
struct NGSegmentString{T,S}
    points::Vector{Tuple{T,T}}
    source::S
    had_repeated_coordinates::Bool
    is_zero_length::Bool
end

"""
    NGPointSource

Source metadata attached to a point extracted from point-like input.
"""
struct NGPointSource{G}
    input_side::NGInputSide
    source_dimension::TopologicalDimension
    element_id::Int
    geometry::G
end

"""
    NGExtractedPoint

Point coordinate plus provenance for point-only and point-vs-geometry paths.
"""
struct NGExtractedPoint{T,S}
    point::Tuple{T,T}
    source::S
end

mutable struct NGExtractionState
    element_id::Int
end

NGExtractionState() = NGExtractionState(0)

function _next_ng_element_id!(state::NGExtractionState)
    state.element_id += 1
    return state.element_id
end

"""
    ng_source_dimension(geom)

Return the maximum topological dimension visible in a GeoInterface geometry.
"""
ng_source_dimension(geom) = _ng_source_dimension(GI.trait(geom), geom)

function _ng_source_dimension(::Nothing, iterable)
    dim = dim_false
    for geom in iterable
        dim = max_dimension(dim, ng_source_dimension(geom))
    end
    return dim
end

_ng_source_dimension(::GI.FeatureCollectionTrait, fc) =
    _ng_source_dimension(nothing, GI.getfeature(fc))
_ng_source_dimension(::GI.FeatureTrait, feature) = ng_source_dimension(GI.geometry(feature))
_ng_source_dimension(::GI.PointTrait, geom) = GI.isempty(geom) ? dim_false : dim_point
_ng_source_dimension(::GI.MultiPointTrait, geom) = GI.isempty(geom) ? dim_false : dim_point
_ng_source_dimension(::GI.AbstractCurveTrait, geom) = GI.isempty(geom) ? dim_false : dim_line
_ng_source_dimension(::GI.PolygonTrait, geom) = GI.isempty(geom) ? dim_false : dim_area

function _ng_source_dimension(::GI.AbstractGeometryTrait, geom)
    GI.isempty(geom) && return dim_false
    dim = dim_false
    for child in GI.getgeom(geom)
        dim = max_dimension(dim, ng_source_dimension(child))
    end
    return dim
end

"""
    ng_is_clockwise(points)

Return `true` if a ring coordinate sequence has clockwise planar orientation.
"""
function ng_is_clockwise(points::AbstractVector)
    length(points) < 2 && return false

    return _ng_orientation_sum(points) > 0.0
end

function _ng_orientation_sum(points::AbstractVector)
    isempty(points) && return 0.0

    total = 0.0
    prev = first(points)
    for point in points
        total += (point[1] - prev[1]) * (point[2] + prev[2])
        prev = point
    end
    return total
end

ng_is_clockwise(curve) = ng_is_clockwise(_ng_curve_points(curve, Float64))

function _ng_ring_orientation(points::AbstractVector)
    orientation_sum = _ng_orientation_sum(points)
    orientation_sum > 0.0 && return ring_clockwise
    orientation_sum < 0.0 && return ring_counterclockwise
    return ring_orientation_none
end

"""
    remove_repeated_ng_points(points)

Remove consecutive duplicate 2D points and report whether any were removed.
"""
function remove_repeated_ng_points(points::AbstractVector{<:Tuple})
    isempty(points) && return collect(points), false

    cleaned = [first(points)]
    had_repeated = false
    for point in Iterators.drop(points, 1)
        if point == last(cleaned)
            had_repeated = true
        else
            push!(cleaned, point)
        end
    end
    return cleaned, had_repeated
end

"""
    ng_is_zero_length(points)

Return `true` when a coordinate sequence has no non-degenerate segment.
"""
ng_is_zero_length(points::AbstractVector{<:Tuple}) = length(points) < 2

function _ng_curve_points(curve, ::Type{T}) where {T}
    return [_tuple_point(GI.getpoint(curve, i), T) for i in 1:GI.npoint(curve)]
end

_ng_extent_intersects(::Nothing, geom) = true
_ng_extent_intersects(extent::Extents.Extent, geom) = Extents.intersects(extent, GI.extent(geom))

function _ng_depth_delta(points::AbstractVector{<:Tuple}, ring_role::NGRingRole)
    ring_role == ring_none && return Int8(0)

    is_clockwise = ng_is_clockwise(points)
    is_canonical = ring_role == ring_shell ? is_clockwise : !is_clockwise
    return Int8(is_canonical ? 1 : -1)
end

function _ng_orient_ring_points(points::Vector{Tuple{T,T}}, ring_role::NGRingRole, orient_rings::Symbol) where {T}
    orient_rings == :source && return points, false
    orient_rings == :relateng ||
        throw(ArgumentError("orient_rings must be either :source or :relateng."))

    require_clockwise = ring_role == ring_shell
    needs_reverse = ng_is_clockwise(points) != require_clockwise
    needs_reverse || return points, false
    return reverse(points), true
end

function _ng_segment_string(
    points::Vector{Tuple{T,T}},
    source::NGSegmentSource;
    clean_repeated::Bool = true,
) where {T}
    cleaned, had_repeated = clean_repeated ? remove_repeated_ng_points(points) : (points, false)
    return NGSegmentString(cleaned, source, had_repeated, ng_is_zero_length(cleaned))
end

"""
    extract_ng_segment_strings(geom; input_side = input_a, T = Float64,
                               extent = nothing, orient_rings = :source)

Extract line and polygon-ring segment strings with NG provenance metadata.
"""
function extract_ng_segment_strings(
    geom,
    ::Type{T} = Float64;
    input_side::NGInputSide = input_a,
    extent = nothing,
    orient_rings::Symbol = :source,
) where {T}
    segments = NGSegmentString[]
    state = NGExtractionState()
    _extract_ng_segment_strings!(
        segments,
        state,
        GI.trait(geom),
        geom;
        input_side,
        T,
        extent,
        orient_rings,
        parent_polygonal = nothing,
    )
    return segments
end

function _extract_ng_segment_strings!(
    segments,
    state::NGExtractionState,
    ::Nothing,
    iterable;
    input_side,
    T,
    extent,
    orient_rings,
    parent_polygonal,
)
    for geom in iterable
        _extract_ng_segment_strings!(
            segments,
            state,
            GI.trait(geom),
            geom;
            input_side,
            T,
            extent,
            orient_rings,
            parent_polygonal,
        )
    end
    return segments
end

function _extract_ng_segment_strings!(
    segments,
    state::NGExtractionState,
    ::GI.FeatureCollectionTrait,
    fc;
    kwargs...,
)
    return _extract_ng_segment_strings!(segments, state, nothing, GI.getfeature(fc); kwargs...)
end

function _extract_ng_segment_strings!(
    segments,
    state::NGExtractionState,
    ::GI.FeatureTrait,
    feature;
    kwargs...,
)
    geom = GI.geometry(feature)
    return _extract_ng_segment_strings!(segments, state, GI.trait(geom), geom; kwargs...)
end

_extract_ng_segment_strings!(segments, state::NGExtractionState, ::GI.PointTrait, geom; kwargs...) =
    segments
_extract_ng_segment_strings!(segments, state::NGExtractionState, ::GI.MultiPointTrait, geom; kwargs...) =
    segments

function _extract_ng_segment_strings!(
    segments,
    state::NGExtractionState,
    ::GI.AbstractCurveTrait,
    curve;
    input_side,
    T,
    extent,
    orient_rings,
    parent_polygonal,
)
    GI.isempty(curve) && return segments
    _ng_extent_intersects(extent, curve) || return segments

    element_id = _next_ng_element_id!(state)
    points = _ng_curve_points(curve, T)
    source = NGSegmentSource(
        input_side,
        dim_line,
        element_id,
        NG_NO_RING_ID,
        ring_none,
        ring_orientation_none,
        Int8(0),
        false,
        curve,
        nothing,
    )
    push!(segments, _ng_segment_string(points, source))
    return segments
end

function _extract_ng_segment_strings!(
    segments,
    state::NGExtractionState,
    ::GI.PolygonTrait,
    polygon;
    input_side,
    T,
    extent,
    orient_rings,
    parent_polygonal,
)
    GI.isempty(polygon) && return segments
    _ng_extent_intersects(extent, polygon) || return segments

    element_id = _next_ng_element_id!(state)
    polygonal_parent = isnothing(parent_polygonal) ? polygon : parent_polygonal
    _add_ng_ring_segment!(
        segments,
        GI.getexterior(polygon),
        input_side,
        T,
        extent,
        orient_rings,
        polygon,
        polygonal_parent,
        element_id,
        0,
        ring_shell,
    )

    for (i, hole) in enumerate(GI.gethole(polygon))
        _add_ng_ring_segment!(
            segments,
            hole,
            input_side,
            T,
            extent,
            orient_rings,
            polygon,
            polygonal_parent,
            element_id,
            i,
            ring_hole,
        )
    end
    return segments
end

function _extract_ng_segment_strings!(
    segments,
    state::NGExtractionState,
    ::GI.MultiPolygonTrait,
    multipolygon;
    input_side,
    T,
    extent,
    orient_rings,
    parent_polygonal,
)
    GI.isempty(multipolygon) && return segments
    _ng_extent_intersects(extent, multipolygon) || return segments

    for polygon in GI.getgeom(multipolygon)
        _extract_ng_segment_strings!(
            segments,
            state,
            GI.trait(polygon),
            polygon;
            input_side,
            T,
            extent,
            orient_rings,
            parent_polygonal = multipolygon,
        )
    end
    return segments
end

function _extract_ng_segment_strings!(
    segments,
    state::NGExtractionState,
    ::GI.AbstractGeometryTrait,
    geom;
    input_side,
    T,
    extent,
    orient_rings,
    parent_polygonal,
)
    GI.isempty(geom) && return segments
    _ng_extent_intersects(extent, geom) || return segments

    for child in GI.getgeom(geom)
        _extract_ng_segment_strings!(
            segments,
            state,
            GI.trait(child),
            child;
            input_side,
            T,
            extent,
            orient_rings,
            parent_polygonal = nothing,
        )
    end
    return segments
end

function _add_ng_ring_segment!(
    segments,
    ring,
    input_side::NGInputSide,
    ::Type{T},
    extent,
    orient_rings::Symbol,
    polygon,
    parent_polygonal,
    element_id::Integer,
    ring_id::Integer,
    ring_role::NGRingRole,
) where {T}
    GI.isempty(ring) && return segments
    _ng_extent_intersects(extent, ring) || return segments

    source_points = _ng_curve_points(ring, T)
    source_orientation = _ng_ring_orientation(source_points)
    depth_delta = _ng_depth_delta(source_points, ring_role)
    points, coordinates_reversed = _ng_orient_ring_points(source_points, ring_role, orient_rings)
    source = NGSegmentSource(
        input_side,
        dim_area,
        element_id,
        ring_id,
        ring_role,
        source_orientation,
        depth_delta,
        coordinates_reversed,
        polygon,
        parent_polygonal,
    )
    push!(segments, _ng_segment_string(points, source))
    return segments
end

"""
    extract_ng_points(geom; input_side = input_a, T = Float64)

Extract point components with provenance, ignoring line and area components.
"""
function extract_ng_points(geom, ::Type{T} = Float64; input_side::NGInputSide = input_a) where {T}
    points = NGExtractedPoint[]
    state = NGExtractionState()
    _extract_ng_points!(points, state, GI.trait(geom), geom; input_side, T)
    return points
end

function _extract_ng_points!(points, state::NGExtractionState, ::Nothing, iterable; input_side, T)
    for geom in iterable
        _extract_ng_points!(points, state, GI.trait(geom), geom; input_side, T)
    end
    return points
end

_extract_ng_points!(points, state::NGExtractionState, ::GI.FeatureCollectionTrait, fc; kwargs...) =
    _extract_ng_points!(points, state, nothing, GI.getfeature(fc); kwargs...)

function _extract_ng_points!(points, state::NGExtractionState, ::GI.FeatureTrait, feature; kwargs...)
    geom = GI.geometry(feature)
    return _extract_ng_points!(points, state, GI.trait(geom), geom; kwargs...)
end

function _extract_ng_points!(points, state::NGExtractionState, ::GI.PointTrait, point; input_side, T)
    GI.isempty(point) && return points

    element_id = _next_ng_element_id!(state)
    source = NGPointSource(input_side, dim_point, element_id, point)
    push!(points, NGExtractedPoint(_tuple_point(point, T), source))
    return points
end

function _extract_ng_points!(points, state::NGExtractionState, ::GI.MultiPointTrait, multipoint; input_side, T)
    GI.isempty(multipoint) && return points

    for point in GI.getpoint(multipoint)
        element_id = _next_ng_element_id!(state)
        source = NGPointSource(input_side, dim_point, element_id, multipoint)
        push!(points, NGExtractedPoint(_tuple_point(point, T), source))
    end
    return points
end

_extract_ng_points!(points, state::NGExtractionState, ::GI.AbstractCurveTrait, geom; kwargs...) = points
_extract_ng_points!(points, state::NGExtractionState, ::GI.PolygonTrait, geom; kwargs...) = points

function _extract_ng_points!(points, state::NGExtractionState, ::GI.AbstractGeometryTrait, geom; input_side, T)
    GI.isempty(geom) && return points

    for child in GI.getgeom(geom)
        _extract_ng_points!(points, state, GI.trait(child), child; input_side, T)
    end
    return points
end
