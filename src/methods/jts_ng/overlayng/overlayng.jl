# # OverlayNG point-dispatch substrate
#
# Keep OverlayNG graph edges and labels out of `common/`.  Overlay labels carry
# operation-specific effective dimension/location and collapse state that should
# not be folded into RelateNG's node topology model.

@enum OverlayOpCode::Int8 begin
    overlay_intersection = 1
    overlay_union = 2
    overlay_difference = 3
    overlay_symdifference = 4
end

@enum OverlayLineState::Int8 begin
    overlay_not_part = 0
    overlay_line_part = 1
    overlay_boundary_part = 2
end

@enum OverlayCollapseRole::Int8 begin
    overlay_not_collapsed = 0
    overlay_collapsed = 1
end

"""
    OverlayInputGeometry(alg, geom)

OverlayNG input wrapper carrying dimension and point-location helpers.
"""
struct OverlayInputGeometry{G,L,C}
    geom::G
    dimension::TopologicalDimension
    locator::L
    segment_strings_cache::C
end

function OverlayInputGeometry(alg::OverlayNG, geom)
    dimension = ng_source_dimension(geom)
    locator = dimension_value(dimension) > dimension_value(dim_point) ?
        RelatePointLocator(geom) :
        nothing
    return OverlayInputGeometry(geom, dimension, locator, Dict{Any,Any}())
end

"""
    OverlayEdgeSourceInfo

OverlayNG source metadata attached to an extracted edge string.
"""
struct OverlayEdgeSourceInfo{G,P}
    input_side::NGInputSide
    source_dimension::TopologicalDimension
    element_id::Int
    ring_id::Int
    ring_role::NGRingRole
    source_orientation::NGRingOrientation
    depth_delta::Int8
    coordinates_reversed::Bool
    is_collapsed::Bool
    geometry::G
    parent_polygonal::P
end

"""
    OverlaySegmentString

OverlayNG edge coordinate sequence with overlay-specific source metadata.
"""
struct OverlaySegmentString{T,S}
    points::Vector{Tuple{T,T}}
    source::S
    had_repeated_coordinates::Bool
    is_zero_length::Bool
end

struct OverlaySegmentRecord{S,E,X}
    segment::S
    segment_index::Int
    edge_index::Int
    edge::E
    extent::X
end

"""
    OverlayEdgeKey

Direction-independent key for a noded OverlayNG edge.
"""
struct OverlayEdgeKey{T}
    p1::Tuple{T,T}
    p2::Tuple{T,T}
end

function OverlayEdgeKey(p1, p2, ::Type{T} = Float64) where {T}
    p1 = _tuple_point(p1, T)
    p2 = _tuple_point(p2, T)
    p1 <= p2 && return OverlayEdgeKey(p1, p2)
    return OverlayEdgeKey(p2, p1)
end

"""
    OverlayEdge

Merged noded edge with all coincident source contributions.
"""
mutable struct OverlayEdge{T,K}
    key::K
    points::Vector{Tuple{T,T}}
    sources::Vector{Any}
    source_directions::Vector{Bool}
    depth_delta::Int
    is_collapsed::Bool
end

"""
    OverlayInputLabel

Per-input OverlayNG label state for one merged edge.
"""
struct OverlayInputLabel
    dimension::TopologicalDimension
    on_location::TopologicalLocation
    left_location::TopologicalLocation
    right_location::TopologicalLocation
    line_state::OverlayLineState
    collapse_role::OverlayCollapseRole
end

"""
    OverlayLabel

Overlay-specific per-input label pair, kept separate from RelateNG topology.
"""
struct OverlayLabel
    input_a::OverlayInputLabel
    input_b::OverlayInputLabel
end

"""
    OverlayHalfEdge

Directed graph edge used by later OverlayNG node-star and ring phases.
"""
mutable struct OverlayHalfEdge{T,E,L}
    origin::Tuple{T,T}
    destination::Tuple{T,T}
    edge::E
    label::L
    angle::Float64
    sym::Any
    next::Any
    prev::Any
    result_area::Bool
    result_line::Bool
    visited::Bool
    ring_id::Int
end

"""
    OverlayGraph

Half-edge graph with angularly sorted node stars.
"""
mutable struct OverlayGraph
    half_edges::Vector{Any}
    node_stars::Dict{Any,Vector{Any}}
end

OverlayGraph() = OverlayGraph(Any[], Dict{Any,Vector{Any}}())

function OverlayEdgeSourceInfo(segment::NGSegmentString)
    source = segment.source
    return OverlayEdgeSourceInfo(
        source.input_side,
        source.source_dimension,
        source.element_id,
        source.ring_id,
        source.ring_role,
        source.source_orientation,
        source.depth_delta,
        source.coordinates_reversed,
        source.source_dimension == dim_area && segment.is_zero_length,
        source.geometry,
        source.parent_polygonal,
    )
end

function OverlaySegmentString(segment::NGSegmentString{T}) where {T}
    return OverlaySegmentString(
        segment.points,
        OverlayEdgeSourceInfo(segment),
        segment.had_repeated_coordinates,
        segment.is_zero_length,
    )
end

"""
    overlay_segment_strings(input, [T]; input_side = input_a, extent = nothing)

Extract and cache OverlayNG-oriented segment strings for graph overlay phases.
"""
function overlay_segment_strings(
    input::OverlayInputGeometry,
    ::Type{T} = Float64;
    input_side::NGInputSide = input_a,
    extent = nothing,
) where {T}
    key = (T, input_side, extent)
    return get!(input.segment_strings_cache, key) do
        map(
            OverlaySegmentString,
            extract_ng_segment_strings(input.geom, T; input_side, extent, orient_rings = :source),
        )
    end
end

function OverlayEdge(segment::OverlaySegmentString{T}) where {T}
    key = OverlayEdgeKey(first(segment.points), last(segment.points), T)
    is_forward = overlay_key_direction(key, segment)
    points = is_forward ? copy(segment.points) : reverse(segment.points)
    source = segment.source
    return OverlayEdge{T,typeof(key)}(
        key,
        points,
        Any[source],
        Bool[is_forward],
        overlay_depth_delta(source, is_forward),
        source.is_collapsed,
    )
end

function overlay_key_direction(key::OverlayEdgeKey, segment::OverlaySegmentString)
    return _tuple_point(first(segment.points), eltype(key.p1)) == key.p1 &&
        _tuple_point(last(segment.points), eltype(key.p2)) == key.p2
end

overlay_depth_delta(source::OverlayEdgeSourceInfo, is_forward::Bool) =
    is_forward ? Int(source.depth_delta) : -Int(source.depth_delta)

function overlay_add_source!(edge::OverlayEdge, segment::OverlaySegmentString)
    is_forward = overlay_key_direction(edge.key, segment)
    source = segment.source
    push!(edge.sources, source)
    push!(edge.source_directions, is_forward)
    edge.depth_delta += overlay_depth_delta(source, is_forward)
    edge.is_collapsed |= source.is_collapsed
    return edge
end

"""
    overlay_merge_edges(segments)

Merge coincident noded segment strings by direction-independent edge key.
"""
function overlay_merge_edges(segments)
    edges = OverlayEdge[]
    edge_indices = Dict{Any,Int}()
    for segment in segments
        length(segment.points) == 2 || throw(ArgumentError("Overlay edge merging requires noded two-point segment strings."))
        key = OverlayEdgeKey(first(segment.points), last(segment.points), eltype(first(segment.points)))
        edge_index = get(edge_indices, key, nothing)
        if isnothing(edge_index)
            push!(edges, OverlayEdge(segment))
            edge_indices[key] = length(edges)
        else
            overlay_add_source!(edges[edge_index], segment)
        end
    end
    return edges
end

function overlay_merge_edges(alg::OverlayNG, geom_a, geom_b, ::Type{T} = Float64; exact = True()) where {T}
    return overlay_merge_edges(overlay_node_segment_strings(alg, geom_a, geom_b, T; exact))
end

overlay_is_line_edge(edge::OverlayEdge) =
    any(source -> source.source_dimension == dim_line || source.is_collapsed, edge.sources)

overlay_is_boundary_edge(edge::OverlayEdge) =
    any(source -> source.source_dimension == dim_area && !source.is_collapsed, edge.sources)

function overlay_primary_ring_role(edge::OverlayEdge)
    any(source -> source.ring_role == ring_shell, edge.sources) && return ring_shell
    any(source -> source.ring_role == ring_hole, edge.sources) && return ring_hole
    return ring_none
end

"""
    overlay_label(edge)

Build the local OverlayNG label pair implied by an edge's coincident sources.
"""
function overlay_label(edge::OverlayEdge)
    return OverlayLabel(
        overlay_input_label(edge, input_a),
        overlay_input_label(edge, input_b),
    )
end

function overlay_input_label(edge::OverlayEdge, input_side::NGInputSide)
    label = overlay_empty_input_label()
    for (source, is_forward) in zip(edge.sources, edge.source_directions)
        source.input_side == input_side || continue
        label = overlay_merge_input_label(label, overlay_source_label(source, is_forward))
    end
    return label
end

overlay_empty_input_label() = OverlayInputLabel(
    dim_false,
    loc_exterior,
    loc_exterior,
    loc_exterior,
    overlay_not_part,
    overlay_not_collapsed,
)

function overlay_source_label(source::OverlayEdgeSourceInfo, is_forward::Bool)
    collapse_role = source.is_collapsed ? overlay_collapsed : overlay_not_collapsed
    if source.source_dimension == dim_area && !source.is_collapsed
        left_location, right_location = overlay_boundary_locations(source, is_forward)
        return OverlayInputLabel(
            dim_area,
            loc_interior,
            left_location,
            right_location,
            overlay_boundary_part,
            collapse_role,
        )
    elseif source.source_dimension == dim_line || source.is_collapsed
        return OverlayInputLabel(
            dim_line,
            loc_interior,
            loc_exterior,
            loc_exterior,
            overlay_line_part,
            collapse_role,
        )
    end
    return overlay_empty_input_label()
end

function overlay_boundary_locations(source::OverlayEdgeSourceInfo, is_forward::Bool)
    delta = overlay_depth_delta(source, is_forward)
    delta < 0 && return loc_interior, loc_exterior
    delta > 0 && return loc_exterior, loc_interior
    return loc_exterior, loc_exterior
end

function overlay_merge_input_label(a::OverlayInputLabel, b::OverlayInputLabel)
    return OverlayInputLabel(
        max_dimension(a.dimension, b.dimension),
        overlay_merge_location(a.on_location, b.on_location),
        overlay_merge_location(a.left_location, b.left_location),
        overlay_merge_location(a.right_location, b.right_location),
        overlay_merge_line_state(a.line_state, b.line_state),
        overlay_merge_collapse_role(a.collapse_role, b.collapse_role),
    )
end

function overlay_merge_location(a::TopologicalLocation, b::TopologicalLocation)
    (a == loc_interior || b == loc_interior) && return loc_interior
    (a == loc_boundary || b == loc_boundary) && return loc_boundary
    return loc_exterior
end

function overlay_merge_line_state(a::OverlayLineState, b::OverlayLineState)
    (a == overlay_boundary_part || b == overlay_boundary_part) && return overlay_boundary_part
    (a == overlay_line_part || b == overlay_line_part) && return overlay_line_part
    return overlay_not_part
end

overlay_merge_collapse_role(a::OverlayCollapseRole, b::OverlayCollapseRole) =
    (a == overlay_collapsed || b == overlay_collapsed) ? overlay_collapsed : overlay_not_collapsed

function OverlayHalfEdge(edge::OverlayEdge, origin, destination)
    origin = _tuple_point(origin, eltype(first(edge.points)))
    destination = _tuple_point(destination, eltype(first(edge.points)))
    return OverlayHalfEdge(
        origin,
        destination,
        edge,
        overlay_label(edge),
        overlay_half_edge_angle(origin, destination),
        nothing,
        nothing,
        nothing,
        false,
        false,
        false,
        0,
    )
end

overlay_half_edge_angle(origin, destination) =
    atan(destination[2] - origin[2], destination[1] - origin[1])

"""
    overlay_graph(edges)

Build a half-edge graph from merged OverlayNG edges.
"""
function overlay_graph(edges)
    graph = OverlayGraph()
    for edge in edges
        overlay_add_edge_pair!(graph, edge)
    end
    return graph
end

function overlay_add_edge_pair!(graph::OverlayGraph, edge::OverlayEdge)
    forward = OverlayHalfEdge(edge, first(edge.points), last(edge.points))
    reverse = OverlayHalfEdge(edge, last(edge.points), first(edge.points))
    forward.sym = reverse
    reverse.sym = forward
    push!(graph.half_edges, forward, reverse)
    overlay_insert_half_edge!(graph, forward)
    overlay_insert_half_edge!(graph, reverse)
    return graph
end

function overlay_insert_half_edge!(graph::OverlayGraph, half_edge::OverlayHalfEdge)
    star = get!(graph.node_stars, half_edge.origin) do
        Any[]
    end
    push!(star, half_edge)
    sort!(star, by = edge -> edge.angle)
    return half_edge
end

overlay_node_star(graph::OverlayGraph, point) =
    get(graph.node_stars, _tuple_point(point), Any[])

overlay_input_label(label::OverlayLabel, input_side::NGInputSide) =
    input_side == input_a ? label.input_a : label.input_b

function overlay_reverse_input_label(label::OverlayInputLabel)
    return OverlayInputLabel(
        label.dimension,
        label.on_location,
        label.right_location,
        label.left_location,
        label.line_state,
        label.collapse_role,
    )
end

function overlay_directed_label(half_edge::OverlayHalfEdge)
    is_forward = half_edge.origin == first(half_edge.edge.points)
    is_forward && return half_edge.label
    return OverlayLabel(
        overlay_reverse_input_label(half_edge.label.input_a),
        overlay_reverse_input_label(half_edge.label.input_b),
    )
end

function overlay_result_location(
    op::OverlayOpCode,
    location_a::TopologicalLocation,
    location_b::TopologicalLocation,
)
    in_a = location_a == loc_interior
    in_b = location_b == loc_interior
    if op == overlay_intersection
        return in_a && in_b
    elseif op == overlay_union
        return in_a || in_b
    elseif op == overlay_difference
        return in_a && !in_b
    elseif op == overlay_symdifference
        return xor(in_a, in_b)
    end
    throw(ArgumentError("Unknown OverlayNG operation code: $op"))
end

"""
    overlay_mark_result_edges!(graph, op)

Mark local result area and line half-edges from existing OverlayNG labels.
"""
function overlay_mark_result_edges!(graph::OverlayGraph, op::OverlayOpCode)
    overlay_clear_result_marks!(graph)
    overlay_mark_result_area_edges!(graph, op)
    overlay_unmark_duplicate_result_area_edges!(graph)
    overlay_mark_result_line_edges!(graph, op)
    return graph
end

function overlay_clear_result_marks!(graph::OverlayGraph)
    for half_edge in graph.half_edges
        half_edge.result_area = false
        half_edge.result_line = false
    end
    return graph
end

function overlay_mark_result_area_edges!(graph::OverlayGraph, op::OverlayOpCode)
    for half_edge in graph.half_edges
        overlay_is_boundary_edge(half_edge.edge) || continue
        label = overlay_directed_label(half_edge)
        half_edge.result_area = overlay_result_location(
            op,
            label.input_a.right_location,
            label.input_b.right_location,
        )
    end
    return graph
end

function overlay_unmark_duplicate_result_area_edges!(graph::OverlayGraph)
    for half_edge in graph.half_edges
        if half_edge.result_area && half_edge.sym.result_area
            half_edge.result_area = false
            half_edge.sym.result_area = false
        end
    end
    return graph
end

function overlay_mark_result_line_edges!(graph::OverlayGraph, op::OverlayOpCode)
    for half_edge in graph.half_edges
        half_edge.origin == first(half_edge.edge.points) || continue
        overlay_is_line_edge(half_edge.edge) || continue
        overlay_is_boundary_edge(half_edge.edge) && continue

        label = overlay_directed_label(half_edge)
        half_edge.result_line = overlay_result_location(
            op,
            label.input_a.on_location,
            label.input_b.on_location,
        )
    end
    return graph
end

"""
    overlay_node_segment_strings(alg, a, b, [T])

Split OverlayNG segment strings at all mutual and self intersections.
"""
function overlay_node_segment_strings(alg::OverlayNG, geom_a, geom_b, ::Type{T} = Float64; exact = True()) where {T}
    return overlay_node_segment_strings(
        alg,
        OverlayInputGeometry(alg, geom_a),
        OverlayInputGeometry(alg, geom_b),
        T;
        exact,
    )
end

function overlay_node_segment_strings(
    alg::OverlayNG,
    geom_a_input::OverlayInputGeometry,
    geom_b_input::OverlayInputGeometry,
    ::Type{T} = Float64;
    exact = True(),
) where {T}
    segments = Any[]
    append!(segments, overlay_segment_strings(geom_a_input, T; input_side = input_a))
    append!(segments, overlay_segment_strings(geom_b_input, T; input_side = input_b))
    return overlay_node_segment_strings(alg, segments, T; exact)
end

function overlay_node_segment_strings(
    alg::OverlayNG,
    segments,
    ::Type{T} = Float64;
    exact = True(),
) where {T}
    records = overlay_segment_records(segments, T)
    isempty(records) && return OverlaySegmentString[]

    split_points = _overlay_initial_split_points(records)
    overlay_add_intersection_split_points!(
        split_points,
        records,
        T;
        exact,
        precision_model = alg.precision_model,
    )

    noded = overlay_split_records(records, split_points, T)
    overlay_validate_fully_noded!(
        noded,
        T;
        exact,
        precision_model = alg.precision_model,
    )
    return noded
end

function overlay_segment_records(segments, ::Type{T} = Float64) where {T}
    records = OverlaySegmentRecord[]
    for (segment_index, segment) in enumerate(segments)
        segment.source.is_collapsed && continue
        length(segment.points) < 2 && continue
        for edge_index in 1:(length(segment.points) - 1)
            p1, p2 = segment.points[edge_index], segment.points[edge_index + 1]
            p1 == p2 && continue
            edge = (p1, p2)
            push!(
                records,
                OverlaySegmentRecord(
                    segment,
                    segment_index,
                    edge_index,
                    edge,
                    ng_segment_extent(edge, T),
                ),
            )
        end
    end
    return records
end

_overlay_initial_split_points(records) =
    [Any[record.edge[1], record.edge[2]] for record in records]

function overlay_add_intersection_split_points!(
    split_points,
    records,
    ::Type{T};
    exact,
    precision_model = nothing,
) where {T}
    extents = getproperty.(records, :extent)
    index = NaturalIndexing.NaturalIndex(extents)
    for (i, record_a) in enumerate(records)
        candidate_indices = SpatialTreeInterface.query(index, record_a.extent)
        for j in candidate_indices
            j <= i && continue
            record_b = records[j]
            ng_segments_maybe_intersect(record_a.edge, record_b.edge, T) || continue
            intersection = ng_segment_intersection(
                record_a.edge,
                record_b.edge,
                T;
                exact,
                precision_model,
            )
            ng_has_intersection(intersection) || continue
            for point in ng_intersection_points(intersection)
                push!(split_points[i], point)
                push!(split_points[j], point)
            end
        end
    end
    return split_points
end

function overlay_split_records(records, split_points, ::Type{T} = Float64) where {T}
    noded = OverlaySegmentString[]
    for (record, points) in zip(records, split_points)
        ordered_points = overlay_unique_ordered_split_points(record.edge, points)
        for i in 1:(length(ordered_points) - 1)
            p1, p2 = ordered_points[i], ordered_points[i + 1]
            p1 == p2 && continue
            push!(
                noded,
                OverlaySegmentString(
                    Tuple{T,T}[_tuple_point(p1, T), _tuple_point(p2, T)],
                    record.segment.source,
                    record.segment.had_repeated_coordinates,
                    false,
                ),
            )
        end
    end
    return noded
end

function overlay_unique_ordered_split_points(edge, points)
    sorted_points = sort(collect(points); by = point -> overlay_segment_fraction(edge, point))
    unique_points = Any[]
    for point in sorted_points
        point in unique_points && continue
        push!(unique_points, point)
    end
    return unique_points
end

function overlay_segment_fraction((p1, p2), point)
    dx = p2[1] - p1[1]
    dy = p2[2] - p1[2]
    if abs(dx) >= abs(dy)
        dx == 0 && return zero(dx)
        return (point[1] - p1[1]) / dx
    else
        dy == 0 && return zero(dy)
        return (point[2] - p1[2]) / dy
    end
end

function overlay_validate_fully_noded!(
    segments,
    ::Type{T} = Float64;
    exact = True(),
    precision_model = nothing,
) where {T}
    overlay_is_fully_noded(segments, T; exact, precision_model) && return segments
    throw(ArgumentError("OverlayNG noder produced linework that is not fully noded."))
end

function overlay_is_fully_noded(
    segments,
    ::Type{T} = Float64;
    exact = True(),
    precision_model = nothing,
) where {T}
    records = overlay_segment_records(segments, T)
    length(records) <= 1 && return true

    extents = getproperty.(records, :extent)
    index = NaturalIndexing.NaturalIndex(extents)
    for (i, record_a) in enumerate(records)
        candidate_indices = SpatialTreeInterface.query(index, record_a.extent)
        for j in candidate_indices
            j <= i && continue
            record_b = records[j]
            ng_segments_maybe_intersect(record_a.edge, record_b.edge, T) || continue
            intersection = ng_segment_intersection(
                record_a.edge,
                record_b.edge,
                T;
                exact,
                precision_model,
            )
            ng_has_intersection(intersection) || continue
            for point in ng_intersection_points(intersection)
                (_overlay_is_edge_endpoint(record_a.edge, point, T) &&
                 _overlay_is_edge_endpoint(record_b.edge, point, T)) || return false
            end
        end
    end
    return true
end

function _overlay_is_edge_endpoint((p1, p2), point, ::Type{T}) where {T}
    point = _tuple_point(point, T)
    return point == _tuple_point(p1, T) || point == _tuple_point(p2, T)
end

overlay(alg::OverlayNG, op::OverlayOpCode, geom_a, geom_b, ::Type{T} = Float64; target = nothing) where {T <: AbstractFloat} =
    overlay(alg, op, OverlayInputGeometry(alg, geom_a), OverlayInputGeometry(alg, geom_b), T; target)

function overlay(
    alg::OverlayNG,
    op::OverlayOpCode,
    input_a::OverlayInputGeometry,
    input_b::OverlayInputGeometry,
    ::Type{T} = Float64;
    target = nothing,
) where {T <: AbstractFloat}
    if overlay_has_point_dispatch(input_a, input_b)
        return overlay_compute_point_dispatch(alg, op, input_a, input_b, T; target)
    end
    throw(ArgumentError("OverlayNG edge overlay is not implemented yet."))
end

intersection(alg::OverlayNG, geom_a, geom_b, ::Type{T} = Float64; target = nothing, kwargs...) where {T <: AbstractFloat} =
    overlay(alg, overlay_intersection, geom_a, geom_b, T; target)

union(alg::OverlayNG, geom_a, geom_b, ::Type{T} = Float64; target = nothing, kwargs...) where {T <: AbstractFloat} =
    overlay(alg, overlay_union, geom_a, geom_b, T; target)

difference(alg::OverlayNG, geom_a, geom_b, ::Type{T} = Float64; target = nothing, kwargs...) where {T <: AbstractFloat} =
    overlay(alg, overlay_difference, geom_a, geom_b, T; target)

symdifference(alg::OverlayNG, geom_a, geom_b, ::Type{T} = Float64; target = nothing, kwargs...) where {T <: AbstractFloat} =
    overlay(alg, overlay_symdifference, geom_a, geom_b, T; target)

overlay_is_pointlike(input::OverlayInputGeometry) =
    input.dimension == dim_false || input.dimension == dim_point

overlay_has_point_dispatch(input_a::OverlayInputGeometry, input_b::OverlayInputGeometry) =
    overlay_is_pointlike(input_a) || overlay_is_pointlike(input_b)

function overlay_compute_point_dispatch(
    alg::OverlayNG,
    op::OverlayOpCode,
    input_a::OverlayInputGeometry,
    input_b::OverlayInputGeometry,
    ::Type{T};
    target = nothing,
) where {T}
    if overlay_is_pointlike(input_a) && overlay_is_pointlike(input_b)
        return overlay_compute_point_point(alg, op, input_a, input_b, T; target)
    elseif overlay_is_pointlike(input_a)
        return overlay_compute_point_nonpoint(alg, op, input_a, input_b, true, T; target)
    else
        return overlay_compute_point_nonpoint(alg, op, input_b, input_a, false, T; target)
    end
end

function overlay_compute_point_point(
    alg::OverlayNG,
    op::OverlayOpCode,
    input_a::OverlayInputGeometry,
    input_b::OverlayInputGeometry,
    ::Type{T};
    target = nothing,
) where {T}
    points_a = overlay_unique_points(input_a, T)
    points_b = overlay_unique_points(input_b, T)
    set_a = Set(points_a)
    set_b = Set(points_b)

    points = if op == overlay_intersection
        [point for point in points_a if point in set_b]
    elseif op == overlay_union
        overlay_union_points(points_a, points_b)
    elseif op == overlay_difference
        [point for point in points_a if !(point in set_b)]
    elseif op == overlay_symdifference
        overlay_union_points(
            [point for point in points_a if !(point in set_b)],
            [point for point in points_b if !(point in set_a)],
        )
    else
        throw(ArgumentError("Unknown OverlayNG operation code: $op"))
    end
    return overlay_filter_results(alg, target, overlay_point_geometries(points))
end

function overlay_compute_point_nonpoint(
    alg::OverlayNG,
    op::OverlayOpCode,
    point_input::OverlayInputGeometry,
    nonpoint_input::OverlayInputGeometry,
    point_is_a::Bool,
    ::Type{T};
    target = nothing,
) where {T}
    covered_points, exterior_points = overlay_partition_points(point_input, nonpoint_input, T)

    if op == overlay_intersection
        return overlay_filter_results(alg, target, overlay_point_geometries(covered_points))
    elseif op == overlay_union || op == overlay_symdifference
        return overlay_filter_results(
            alg,
            target,
            Any[nonpoint_input.geom, overlay_point_geometries(exterior_points)...],
        )
    elseif op == overlay_difference
        if point_is_a
            return overlay_filter_results(alg, target, overlay_point_geometries(exterior_points))
        else
            return overlay_filter_results(alg, target, Any[nonpoint_input.geom])
        end
    end
    throw(ArgumentError("Unknown OverlayNG operation code: $op"))
end

function overlay_partition_points(point_input::OverlayInputGeometry, nonpoint_input::OverlayInputGeometry, ::Type{T}) where {T}
    covered_points = Any[]
    exterior_points = Any[]
    for point in overlay_unique_points(point_input, T)
        target_location = relate_locate_with_dim(nonpoint_input.locator, point)
        if target_location.location == loc_exterior
            push!(exterior_points, point)
        else
            push!(covered_points, point)
        end
    end
    return covered_points, exterior_points
end

function overlay_unique_points(input::OverlayInputGeometry, ::Type{T}) where {T}
    seen = Set{Any}()
    points = Any[]
    for extracted in extract_ng_points(input.geom, T)
        point = extracted.point
        point in seen && continue
        push!(seen, point)
        push!(points, point)
    end
    return points
end

function overlay_union_points(points_a, points_b)
    seen = Set{Any}()
    points = Any[]
    for point in Iterators.flatten((points_a, points_b))
        point in seen && continue
        push!(seen, point)
        push!(points, point)
    end
    return points
end

overlay_point_geometries(points) = Any[GI.Point(point[1], point[2]) for point in points]

function overlay_filter_results(alg::OverlayNG, target, geometries)
    target_trait = isnothing(target) ? nothing : TraitTarget(target)
    results = Any[]
    for geom in geometries
        overlay_accepts_result(alg, target_trait, geom) && push!(results, geom)
    end
    return results
end

function overlay_accepts_result(alg::OverlayNG, target, geom)
    if alg.area_result_only && ng_source_dimension(geom) != dim_area
        return false
    end
    return overlay_accepts_target(target, geom)
end

overlay_accepts_target(::Nothing, geom) = true
overlay_accepts_target(::TraitTarget{Nothing}, geom) = true
overlay_accepts_target(::TraitTarget{Target}, geom) where {Target} =
    GI.trait(geom) isa Target
