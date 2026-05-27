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
    OverlaySnapNode

JTS-style node location on a snap-rounding segment string.
"""
mutable struct OverlaySnapNode{T}
    point::Tuple{T,T}
    segment_index::Int
end

"""
    OverlaySnapSegmentString

Mutable noding wrapper used by the fixed-precision snap-rounding path.
"""
mutable struct OverlaySnapSegmentString{T,S}
    points::Vector{Tuple{T,T}}
    source::S
    had_repeated_coordinates::Bool
    nodes::Vector{OverlaySnapNode{T}}
end

"""
    OverlayHotPixel

Partially-open snap-rounding pixel centered on a rounded grid coordinate.
"""
mutable struct OverlayHotPixel{T}
    coordinate::Tuple{T,T}
    scale::T
    hpx::T
    hpy::T
    is_node::Bool
end

"""
    OverlayHotPixelIndex

Unique hot-pixel registry with JTS node/non-node pixel state.
"""
mutable struct OverlayHotPixelIndex{T,P}
    precision_model::P
    scale::T
    pixels::Vector{OverlayHotPixel{T}}
    pixel_indices::Dict{Tuple{T,T},Int}
end

Base.copy(node::OverlaySnapNode) = OverlaySnapNode(node.point, node.segment_index)

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

function OverlayEdgeKey(points::AbstractVector{<:Tuple}, ::Type{T} = Float64) where {T}
    length(points) >= 2 || throw(ArgumentError("Overlay edge keys require at least two points."))
    if overlay_edge_direction(points)
        return OverlayEdgeKey(_tuple_point(points[1], T), _tuple_point(points[2], T))
    end
    return OverlayEdgeKey(_tuple_point(points[end], T), _tuple_point(points[end - 1], T))
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
    on_location::Union{Nothing,TopologicalLocation}
    left_location::Union{Nothing,TopologicalLocation}
    right_location::Union{Nothing,TopologicalLocation}
    line_state::OverlayLineState
    collapse_role::OverlayCollapseRole
    ring_role::NGRingRole
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
    OverlayResultRing

Closed result-area ring built from marked OverlayNG half-edges.
"""
mutable struct OverlayResultRing{T}
    points::Vector{Tuple{T,T}}
    is_hole::Bool
    holes::Vector{Any}
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
    is_forward::Bool
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
    node_edges::Dict{Any,Any}
end

OverlayGraph() = OverlayGraph(Any[], Dict{Any,Vector{Any}}(), Dict{Any,Any}())

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
    key = OverlayEdgeKey(segment.points, T)
    source = segment.source
    return OverlayEdge{T,typeof(key)}(
        key,
        copy(segment.points),
        Any[source],
        Bool[true],
        overlay_depth_delta(source, true),
        source.is_collapsed,
    )
end

function overlay_key_direction(key::OverlayEdgeKey, segment::OverlaySegmentString)
    return _tuple_point(first(segment.points), eltype(key.p1)) == key.p1 &&
        _tuple_point(segment.points[2], eltype(key.p2)) == key.p2
end

overlay_depth_delta(source::OverlayEdgeSourceInfo, is_forward::Bool) =
    is_forward ? Int(source.depth_delta) : -Int(source.depth_delta)

function overlay_add_source!(edge::OverlayEdge, segment::OverlaySegmentString)
    is_forward = overlay_relative_direction(edge, segment)
    source = segment.source
    push!(edge.sources, source)
    push!(edge.source_directions, is_forward)
    edge.depth_delta += overlay_depth_delta(source, is_forward)
    edge.is_collapsed |= source.is_collapsed
    return edge
end

function overlay_relative_direction(edge::OverlayEdge, segment::OverlaySegmentString)
    length(edge.points) == length(segment.points) ||
        throw(ArgumentError("OverlayNG cannot merge differently noded coincident edges."))
    return _tuple_point(segment.points[1], eltype(first(edge.points))) == edge.points[1] &&
        _tuple_point(segment.points[2], eltype(first(edge.points))) == edge.points[2]
end

function overlay_edge_direction(points::AbstractVector)
    length(points) >= 2 || throw(ArgumentError("OverlayNG edge direction requires at least two points."))
    cmp = points[1] == points[end] ? 0 : (points[1] < points[end] ? -1 : 1)
    if cmp == 0
        cmp = points[2] == points[end - 1] ? 0 : (points[2] < points[end - 1] ? -1 : 1)
    end
    cmp == 0 && throw(ArgumentError("OverlayNG edge direction cannot be determined because endpoints are equal."))
    return cmp == -1
end

"""
    overlay_merge_edges(segments)

Merge coincident noded segment strings by direction-independent edge key.
"""
function overlay_merge_edges(segments)
    edges = OverlayEdge[]
    edge_indices = Dict{Any,Int}()
    for segment in segments
        length(segment.points) >= 2 || continue
        key = OverlayEdgeKey(segment.points, eltype(first(segment.points)))
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
    if overlay_is_area_collapse(edge, input_side)
        return overlay_area_collapse_line_label(overlay_primary_ring_role(edge, input_side))
    end

    label = overlay_empty_input_label()
    for (source, is_forward) in zip(edge.sources, edge.source_directions)
        source.input_side == input_side || continue
        label = overlay_merge_input_label(label, overlay_source_label(source, is_forward))
    end
    return label
end

function overlay_is_area_collapse(edge::OverlayEdge, input_side::NGInputSide)
    has_area_source = false
    depth_delta = 0
    for (source, is_forward) in zip(edge.sources, edge.source_directions)
        source.input_side == input_side || continue
        source.source_dimension == dim_area || continue
        has_area_source = true
        depth_delta += overlay_depth_delta(source, is_forward)
    end
    return has_area_source && depth_delta == 0
end

function overlay_primary_ring_role(edge::OverlayEdge, input_side::NGInputSide)
    any(source -> source.input_side == input_side && source.ring_role == ring_shell, edge.sources) && return ring_shell
    any(source -> source.input_side == input_side && source.ring_role == ring_hole, edge.sources) && return ring_hole
    return ring_none
end

overlay_area_collapse_line_label(ring_role::NGRingRole) = OverlayInputLabel(
    dim_line,
    nothing,
    nothing,
    nothing,
    overlay_line_part,
    overlay_collapsed,
    ring_role,
)

overlay_empty_input_label() = OverlayInputLabel(
    dim_false,
    nothing,
    nothing,
    nothing,
    overlay_not_part,
    overlay_not_collapsed,
    ring_none,
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
            source.ring_role,
        )
    elseif source.source_dimension == dim_line || source.is_collapsed
        return OverlayInputLabel(
            dim_line,
            nothing,
            nothing,
            nothing,
            overlay_line_part,
            collapse_role,
            source.is_collapsed ? source.ring_role : ring_none,
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
        overlay_merge_ring_role(a.ring_role, b.ring_role),
    )
end

overlay_merge_location(::Nothing, ::Nothing) = nothing
overlay_merge_location(::Nothing, b::Union{Nothing,TopologicalLocation}) = b
overlay_merge_location(a::Union{Nothing,TopologicalLocation}, ::Nothing) = a

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

function overlay_merge_ring_role(a::NGRingRole, b::NGRingRole)
    (a == ring_shell || b == ring_shell) && return ring_shell
    (a == ring_hole || b == ring_hole) && return ring_hole
    return ring_none
end

function OverlayHalfEdge(edge::OverlayEdge, origin, destination)
    is_forward = origin == first(edge.points) && origin != last(edge.points)
    direction_point = is_forward ? edge.points[2] : edge.points[end - 1]
    return OverlayHalfEdge(edge, origin, destination, direction_point, is_forward)
end

function OverlayHalfEdge(edge::OverlayEdge, origin, destination, direction_point, is_forward::Bool)
    origin = _tuple_point(origin, eltype(first(edge.points)))
    destination = _tuple_point(destination, eltype(first(edge.points)))
    direction_point = _tuple_point(direction_point, eltype(first(edge.points)))
    return OverlayHalfEdge(
        origin,
        destination,
        edge,
        overlay_label(edge),
        overlay_half_edge_angle(origin, direction_point),
        is_forward,
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
    forward = OverlayHalfEdge(edge, first(edge.points), last(edge.points), edge.points[2], true)
    reverse = OverlayHalfEdge(edge, last(edge.points), first(edge.points), edge.points[end - 1], false)
    forward.sym = reverse
    reverse.sym = forward
    push!(graph.half_edges, forward, reverse)
    overlay_insert_half_edge!(graph, forward)
    overlay_insert_half_edge!(graph, reverse)
    return graph
end

function overlay_insert_half_edge!(graph::OverlayGraph, half_edge::OverlayHalfEdge)
    get!(graph.node_edges, half_edge.origin) do
        half_edge
    end
    star = get!(graph.node_stars, half_edge.origin) do
        Any[]
    end
    push!(star, half_edge)
    sort!(star, by = edge -> overlay_jts_angle(edge.angle))
    return half_edge
end

overlay_jts_angle(angle::Real) =
    angle < 0 ? angle + 2π : angle

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
        label.ring_role,
    )
end

function overlay_directed_label(half_edge::OverlayHalfEdge)
    half_edge.is_forward && return half_edge.label
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
function overlay_mark_result_edges!(
    graph::OverlayGraph,
    op::OverlayOpCode;
    strict::Bool = false,
    input_area_side = nothing,
)
    overlay_clear_result_marks!(graph)
    overlay_mark_result_area_edges!(graph, op)
    overlay_unmark_duplicate_result_area_edges!(graph)
    has_result_area = any(half_edge -> half_edge.result_area, graph.half_edges)
    overlay_mark_result_line_edges!(
        graph,
        op;
        strict,
        has_result_area,
        input_area_side,
    )
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
        label = overlay_directed_label(half_edge)
        overlay_is_boundary_either(label) || continue
        half_edge.result_area = overlay_result_location(
            op,
            overlay_boundary_or_line_location(label.input_a, side_right),
            overlay_boundary_or_line_location(label.input_b, side_right),
        )
    end
    return graph
end

function overlay_boundary_or_line_location(label::OverlayInputLabel, side::SidePosition)
    location = if overlay_is_boundary_input(label)
        side == side_left ? label.left_location : label.right_location
    else
        label.on_location
    end
    isnothing(location) && throw(ArgumentError("OverlayNG label location is unknown after labelling."))
    return location
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

function overlay_mark_result_line_edges!(
    graph::OverlayGraph,
    op::OverlayOpCode;
    strict::Bool = false,
    has_result_area::Bool = any(half_edge -> half_edge.result_area, graph.half_edges),
    input_area_side = nothing,
)
    for half_edge in graph.half_edges
        half_edge.is_forward || continue
        overlay_half_edge_in_result_either(half_edge) && continue

        label = overlay_directed_label(half_edge)
        overlay_is_result_line(
            label,
            op,
            strict;
            has_result_area,
            input_area_side,
        ) || continue
        half_edge.result_line = true
        half_edge.sym.result_line = true
    end
    return graph
end

overlay_half_edge_in_result(half_edge::OverlayHalfEdge) =
    half_edge.result_area || half_edge.result_line

overlay_half_edge_in_result_either(half_edge::OverlayHalfEdge) =
    overlay_half_edge_in_result(half_edge) || overlay_half_edge_in_result(half_edge.sym)

overlay_is_boundary_input(label::OverlayInputLabel) =
    label.line_state == overlay_boundary_part && label.collapse_role == overlay_not_collapsed

overlay_is_line_input(label::OverlayInputLabel) =
    label.line_state == overlay_line_part && label.collapse_role == overlay_not_collapsed

overlay_is_collapse_input(label::OverlayInputLabel) =
    label.collapse_role == overlay_collapsed

overlay_is_boundary_singleton(label::OverlayLabel) =
    (overlay_is_boundary_input(label.input_a) && label.input_b.dimension == dim_false) ||
    (overlay_is_boundary_input(label.input_b) && label.input_a.dimension == dim_false)

overlay_is_boundary_both(label::OverlayLabel) =
    overlay_is_boundary_input(label.input_a) && overlay_is_boundary_input(label.input_b)

overlay_is_boundary_either(label::OverlayLabel) =
    overlay_is_boundary_input(label.input_a) || overlay_is_boundary_input(label.input_b)

overlay_is_boundary_collapse(label::OverlayLabel) =
    !overlay_is_line(label) && !overlay_is_boundary_both(label)

overlay_is_boundary_touch(label::OverlayLabel) =
    overlay_is_boundary_both(label) &&
    label.input_a.right_location != label.input_b.right_location

overlay_is_line(label::OverlayLabel) =
    overlay_is_line_input(label.input_a) || overlay_is_line_input(label.input_b)

overlay_is_interior_collapse(label::OverlayLabel) =
    (overlay_is_collapse_input(label.input_a) && label.input_a.on_location == loc_interior) ||
    (overlay_is_collapse_input(label.input_b) && label.input_b.on_location == loc_interior)

overlay_is_collapse_and_not_part_interior(label::OverlayLabel) =
    (overlay_is_collapse_input(label.input_a) &&
     label.input_b.dimension == dim_false &&
     label.input_b.on_location == loc_interior) ||
    (overlay_is_collapse_input(label.input_b) &&
     label.input_a.dimension == dim_false &&
     label.input_a.on_location == loc_interior)

function overlay_is_result_line(
    label::OverlayLabel,
    op::OverlayOpCode,
    strict::Bool;
    has_result_area::Bool,
    input_area_side,
)
    overlay_is_boundary_singleton(label) && return false
    strict && overlay_is_boundary_collapse(label) && return false
    overlay_is_interior_collapse(label) && return false

    if op != overlay_intersection
        overlay_is_collapse_and_not_part_interior(label) && return false
        if has_result_area && !isnothing(input_area_side)
            overlay_is_line_in_area(label, input_area_side) && return false
        end
    end

    !strict && op == overlay_intersection && overlay_is_boundary_touch(label) && return true

    return overlay_result_location(
        op,
        overlay_effective_line_location(label.input_a),
        overlay_effective_line_location(label.input_b),
    )
end

function overlay_effective_line_location(label::OverlayInputLabel)
    if overlay_is_collapse_input(label) || overlay_is_line_input(label)
        return loc_interior
    end
    isnothing(label.on_location) && throw(ArgumentError("OverlayNG line location is unknown after labelling."))
    return label.on_location
end

function overlay_is_line_in_area(label::OverlayLabel, input_side::NGInputSide)
    input_label = overlay_input_label(label, input_side)
    return input_label.on_location == loc_interior
end

overlay_has_edges(input::OverlayInputGeometry) =
    dimension_value(input.dimension) > dimension_value(dim_point)

"""
    overlay_compute_labelling!(graph, input_a, input_b)

Run the JTS OverlayLabeller pass order over a noded overlay graph.
"""
function overlay_compute_labelling!(
    graph::OverlayGraph,
    input_a_geom::OverlayInputGeometry,
    input_b_geom::OverlayInputGeometry,
)
    overlay_label_area_node_edges!(graph, input_a_geom, input_b_geom)
    overlay_label_connected_linear_edges!(graph, input_a_geom, input_b_geom)
    overlay_label_collapsed_edges!(graph)
    overlay_label_connected_linear_edges!(graph, input_a_geom, input_b_geom)
    overlay_label_disconnected_edges!(graph, input_a_geom, input_b_geom)
    return graph
end

function overlay_label_area_node_edges!(
    graph::OverlayGraph,
    input_a_geom::OverlayInputGeometry,
    input_b_geom::OverlayInputGeometry,
)
    for node_edge in values(graph.node_edges)
        input_a_geom.dimension == dim_area && overlay_propagate_area_locations!(graph, node_edge, input_a)
        input_b_geom.dimension == dim_area && overlay_propagate_area_locations!(graph, node_edge, input_b)
    end
    return graph
end

function overlay_propagate_area_locations!(
    graph::OverlayGraph,
    node_edge::OverlayHalfEdge,
    input_side::NGInputSide,
)
    star = overlay_node_star(graph, node_edge.origin)
    length(star) <= 1 && return star

    node_edge_index = findfirst(candidate -> candidate === node_edge, star)
    isnothing(node_edge_index) && return star

    start_index = nothing
    for offset in 0:(length(star) - 1)
        candidate_index = mod1(node_edge_index + offset, length(star))
        half_edge = star[candidate_index]
        if overlay_is_boundary_input(overlay_input_label(overlay_directed_label(half_edge), input_side))
            start_index = candidate_index
            break
        end
    end
    isnothing(start_index) && return star

    start_label = overlay_input_label(overlay_directed_label(star[start_index]), input_side)
    current_location = start_label.left_location
    isnothing(current_location) && throw(ArgumentError("OverlayNG boundary edge has no left location."))

    for offset in 1:(length(star) - 1)
        half_edge = star[mod1(start_index + offset, length(star))]
        input_label = overlay_input_label(overlay_directed_label(half_edge), input_side)
        if overlay_is_boundary_input(input_label)
            input_label.right_location == current_location ||
                throw(ArgumentError("OverlayNG side location conflict during area propagation."))
            isnothing(input_label.left_location) &&
                throw(ArgumentError("OverlayNG boundary edge has no left location."))
            current_location = input_label.left_location
        elseif overlay_is_line_location_unknown(input_label)
            overlay_set_input_line_location!(half_edge, input_side, current_location)
        end
    end
    return star
end

function overlay_label_connected_linear_edges!(
    graph::OverlayGraph,
    input_a_geom::OverlayInputGeometry,
    input_b_geom::OverlayInputGeometry,
)
    overlay_propagate_linear_locations!(graph, input_a; input_is_line = input_a_geom.dimension == dim_line)
    overlay_has_edges(input_b_geom) &&
        overlay_propagate_linear_locations!(graph, input_b; input_is_line = input_b_geom.dimension == dim_line)
    return graph
end

function overlay_propagate_linear_locations!(
    graph::OverlayGraph,
    input_side::NGInputSide;
    input_is_line::Bool,
)
    queue = Any[]
    for half_edge in graph.half_edges
        input_label = overlay_input_label(half_edge.label, input_side)
        if overlay_is_linear_input(input_label) && !overlay_is_line_location_unknown(input_label)
            push!(queue, half_edge)
        end
    end

    while !isempty(queue)
        half_edge = popfirst!(queue)
        input_label = overlay_input_label(half_edge.label, input_side)
        line_location = input_label.on_location
        isnothing(line_location) && continue
        input_is_line && line_location != loc_exterior && continue

        star = overlay_node_star(graph, half_edge.origin)
        start_index = findfirst(candidate -> candidate === half_edge, star)
        isnothing(start_index) && continue
        for offset in 1:(length(star) - 1)
            candidate = star[mod1(start_index + offset, length(star))]
            candidate_label = overlay_input_label(candidate.label, input_side)
            overlay_is_line_location_unknown(candidate_label) || continue
            overlay_set_input_line_location!(candidate, input_side, line_location)
            pushfirst!(queue, candidate.sym)
        end
    end
    return graph
end

overlay_is_linear_input(label::OverlayInputLabel) =
    overlay_is_line_input(label) || overlay_is_collapse_input(label)

overlay_is_line_location_unknown(label::OverlayInputLabel) =
    isnothing(label.on_location)

function overlay_label_collapsed_edges!(graph::OverlayGraph)
    for half_edge in graph.half_edges
        half_edge.is_forward || continue
        overlay_label_collapsed_edge!(half_edge, input_a)
        overlay_label_collapsed_edge!(half_edge, input_b)
    end
    return graph
end

function overlay_label_collapsed_edge!(half_edge::OverlayHalfEdge, input_side::NGInputSide)
    input_label = overlay_input_label(half_edge.label, input_side)
    overlay_is_line_location_unknown(input_label) || return half_edge
    overlay_is_collapse_input(input_label) || return half_edge
    location = input_label.ring_role == ring_hole ? loc_interior : loc_exterior
    overlay_set_input_line_location!(half_edge, input_side, location)
    return half_edge
end

function overlay_set_input_line_location!(
    half_edge::OverlayHalfEdge,
    input_side::NGInputSide,
    location::TopologicalLocation,
)
    label = overlay_replace_input_label(
        half_edge.label,
        input_side,
        overlay_with_line_location(overlay_input_label(half_edge.label, input_side), location),
    )
    half_edge.label = label
    half_edge.sym.label = label
    return half_edge
end

function overlay_replace_input_label(
    label::OverlayLabel,
    input_side::NGInputSide,
    input_label::OverlayInputLabel,
)
    input_side == input_a && return OverlayLabel(input_label, label.input_b)
    return OverlayLabel(label.input_a, input_label)
end

function overlay_with_line_location(label::OverlayInputLabel, location::TopologicalLocation)
    return OverlayInputLabel(
        label.dimension,
        location,
        label.left_location,
        label.right_location,
        label.line_state,
        label.collapse_role,
        label.ring_role,
    )
end

function overlay_with_all_locations(label::OverlayInputLabel, location::TopologicalLocation)
    return OverlayInputLabel(
        label.dimension,
        location,
        location,
        location,
        label.line_state,
        label.collapse_role,
        label.ring_role,
    )
end

function overlay_label_disconnected_edges!(
    graph::OverlayGraph,
    input_a_geom::OverlayInputGeometry,
    input_b_geom::OverlayInputGeometry,
)
    for half_edge in graph.half_edges
        half_edge.is_forward || continue
        overlay_label_disconnected_input!(half_edge, input_a, input_a_geom)
        overlay_label_disconnected_input!(half_edge, input_b, input_b_geom)
    end
    return graph
end

function overlay_label_disconnected_input!(
    half_edge::OverlayHalfEdge,
    input_side::NGInputSide,
    input::OverlayInputGeometry,
)
    input_label = overlay_input_label(half_edge.label, input_side)
    overlay_is_line_location_unknown(input_label) || return half_edge
    location = overlay_disconnected_edge_location(half_edge, input)
    new_label = overlay_replace_input_label(
        half_edge.label,
        input_side,
        overlay_with_all_locations(input_label, location),
    )
    half_edge.label = new_label
    half_edge.sym.label = new_label
    return half_edge
end

function overlay_disconnected_edge_location(half_edge::OverlayHalfEdge, input::OverlayInputGeometry)
    input.dimension == dim_area || return loc_exterior
    origin_location = relate_locate_with_dim(input.locator, half_edge.origin).location
    destination_location = relate_locate_with_dim(input.locator, half_edge.destination).location
    return origin_location != loc_exterior && destination_location != loc_exterior ? loc_interior : loc_exterior
end

function overlay_result_area_half_edges(graph::OverlayGraph)
    return [half_edge for half_edge in graph.half_edges if half_edge.result_area]
end

"""
    overlay_extract_result_polygons(graph)

Build polygon geometries from marked result-area half-edges.
"""
function overlay_extract_result_polygons(graph::OverlayGraph)
    overlay_link_result_area_edges!(graph)
    rings = overlay_result_area_rings!(graph)
    return overlay_result_ring_polygons(rings)
end

function overlay_link_result_area_edges!(graph::OverlayGraph)
    for half_edge in graph.half_edges
        half_edge.next = nothing
        half_edge.prev = nothing
        half_edge.visited = false
        half_edge.ring_id = 0
    end

    for half_edge in overlay_result_area_half_edges(graph)
        next_edge = overlay_next_result_area_edge(graph, half_edge)
        half_edge.next = next_edge
        next_edge.prev = half_edge
    end
    return graph
end

function overlay_next_result_area_edge(graph::OverlayGraph, half_edge::OverlayHalfEdge)
    star = overlay_node_star(graph, half_edge.destination)
    sym_index = findfirst(candidate -> candidate === half_edge.sym, star)
    isnothing(sym_index) && throw(ArgumentError("OverlayNG graph is missing a symmetric half-edge at a result node."))

    for offset in 1:length(star)
        candidate = star[mod1(sym_index + offset, length(star))]
        candidate.result_area && return candidate
    end
    throw(ArgumentError("OverlayNG result area edge has no outgoing continuation."))
end

function overlay_result_area_rings!(graph::OverlayGraph)
    rings = OverlayResultRing[]
    ring_id = 0
    for half_edge in overlay_result_area_half_edges(graph)
        half_edge.visited && continue
        ring_id += 1
        ring = overlay_result_area_ring!(half_edge, ring_id)
        append!(rings, overlay_minimal_result_rings(ring.points))
    end
    return rings
end

function overlay_result_area_ring!(start_edge::OverlayHalfEdge, ring_id::Integer)
    points = typeof(start_edge.origin)[]
    half_edge = start_edge
    while true
        half_edge.visited && throw(ArgumentError("OverlayNG result ring visited an edge twice."))
        half_edge.visited = true
        half_edge.ring_id = ring_id
        overlay_add_coordinate_list!(points, overlay_half_edge_points(half_edge))

        half_edge = half_edge.next
        isnothing(half_edge) && throw(ArgumentError("OverlayNG result ring has an unlinked edge."))
        half_edge === start_edge && break
    end
    last(points) == first(points) || push!(points, first(points))
    return OverlayResultRing(points, !ng_is_clockwise(points), Any[])
end

function overlay_minimal_result_rings(points)
    rings = OverlayResultRing[]
    length(points) < 4 && return rings

    open_points = collect(points[1:(end - 1)])
    extracted_loops = Any[]
    while true
        loop_range = overlay_repeated_vertex_loop_range(open_points)
        isnothing(loop_range) && break

        start_index, end_index = loop_range
        loop_points = collect(open_points[start_index:end_index])
        push!(extracted_loops, loop_points)
        open_points = vcat(open_points[1:start_index], open_points[(end_index + 1):end])
    end

    for loop_points in extracted_loops
        overlay_push_result_ring!(rings, loop_points)
    end
    overlay_push_result_ring!(rings, open_points)
    return rings
end

function overlay_repeated_vertex_loop_range(points)
    seen = Dict{Any,Int}()
    for (index, point) in enumerate(points)
        previous_index = get(seen, point, nothing)
        if !isnothing(previous_index) && index > previous_index + 1
            return previous_index, index
        end
        seen[point] = index
    end
    return nothing
end

function overlay_push_result_ring!(rings, open_points)
    length(open_points) >= 3 || return rings
    closed_points = collect(open_points)
    first(closed_points) == last(closed_points) || push!(closed_points, first(closed_points))
    length(closed_points) >= 4 || return rings
    _ng_orientation_sum(closed_points) == 0.0 && return rings
    push!(rings, OverlayResultRing(closed_points, !ng_is_clockwise(closed_points), Any[]))
    return rings
end

function overlay_result_ring_polygons(rings)
    shells = [ring for ring in rings if !ring.is_hole]
    holes = [ring for ring in rings if ring.is_hole]

    for hole in holes
        shell = overlay_find_containing_shell(hole, shells)
        isnothing(shell) && throw(ArgumentError("OverlayNG polygon extraction found a hole with no containing shell."))
        push!(shell.holes, hole)
    end

    return Any[GI.Polygon(Any[shell.points, getproperty.(shell.holes, :points)...]) for shell in shells]
end

function overlay_find_containing_shell(hole::OverlayResultRing, shells)
    containing_shell = nothing
    containing_area = Inf
    for shell in shells
        overlay_ring_contains_ring(shell, hole) || continue
        shell_area = abs(_ng_orientation_sum(shell.points))
        if shell_area < containing_area
            containing_shell = shell
            containing_area = shell_area
        end
    end
    return containing_shell
end

function overlay_ring_contains_ring(shell::OverlayResultRing, hole::OverlayResultRing)
    locator = RelatePointLocator(GI.Polygon([shell.points]))
    for point in Iterators.drop(hole.points, 1)
        location = relate_locate_with_dim(locator, point).location
        location == loc_interior && return true
        location == loc_exterior && return false
    end
    return false
end

"""
    overlay_extract_result_lines(graph)

Build noded line geometries from marked result-line half-edges.
"""
function overlay_extract_result_lines(graph::OverlayGraph)
    lines = Any[]
    for half_edge in graph.half_edges
        half_edge.result_line || continue
        half_edge.visited && continue
        push!(lines, GI.LineString(overlay_half_edge_points(half_edge)))
        half_edge.visited = true
        half_edge.sym.visited = true
    end
    return lines
end

function overlay_half_edge_points(half_edge::OverlayHalfEdge)
    if half_edge.is_forward
        return copy(half_edge.edge.points)
    end
    return reverse(half_edge.edge.points)
end

"""
    overlay_extract_intersection_points(graph; strict = false)

Build point results for non-point intersection nodes not already in the result.
"""
function overlay_extract_intersection_points(graph::OverlayGraph; strict::Bool = false)
    points = Any[]
    for point in sort(collect(keys(graph.node_stars)))
        star = graph.node_stars[point]
        overlay_is_result_intersection_point(star; strict) || continue
        push!(points, GI.Point(point[1], point[2]))
    end
    return points
end

function overlay_is_result_intersection_point(star; strict::Bool = false)
    is_edge_of_a = false
    is_edge_of_b = false
    for half_edge in star
        overlay_half_edge_in_result(half_edge) && return false
        label = half_edge.label
        is_edge_of_a |= overlay_is_edge_of(label, input_a; strict)
        is_edge_of_b |= overlay_is_edge_of(label, input_b; strict)
    end
    return is_edge_of_a && is_edge_of_b
end

function overlay_is_edge_of(label::OverlayLabel, input_side::NGInputSide; strict::Bool = false)
    strict && overlay_is_boundary_collapse(label) && return false
    input_label = overlay_input_label(label, input_side)
    return overlay_is_boundary_input(input_label) || overlay_is_line_input(input_label)
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
    if alg.precision_model isa FixedPrecisionModel
        return overlay_snapround_node_segment_strings(alg.precision_model, segments, T; exact)
    end

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

const OVERLAY_SNAPROUNDING_NEARNESS_FACTOR = 100

function OverlaySnapSegmentString(segment::OverlaySegmentString{T}) where {T}
    return OverlaySnapSegmentString(
        copy(segment.points),
        segment.source,
        segment.had_repeated_coordinates,
        OverlaySnapNode{T}[],
    )
end

function OverlayHotPixel(point, scale::T) where {T}
    point = _tuple_point(point, T)
    hpx = scale == one(T) ? point[1] : _jts_precision_round(point[1] * scale)
    hpy = scale == one(T) ? point[2] : _jts_precision_round(point[2] * scale)
    return OverlayHotPixel(point, scale, hpx, hpy, false)
end

function OverlayHotPixelIndex(precision_model::FixedPrecisionModel, ::Type{T} = Float64) where {T}
    scale = T(precision_model.scale)
    return OverlayHotPixelIndex(precision_model, scale, OverlayHotPixel{T}[], Dict{Tuple{T,T},Int}())
end

"""
    overlay_snapround_node_segment_strings(model, segments, [T])

Fixed-precision OverlayNG noder following JTS `SnapRoundingNoder` phase order.
"""
function overlay_snapround_node_segment_strings(
    precision_model::FixedPrecisionModel,
    segments,
    ::Type{T} = Float64;
    exact = True(),
) where {T}
    snap_strings = Any[
        OverlaySnapSegmentString(segment) for segment in segments if !segment.source.is_collapsed
    ]
    isempty(snap_strings) && return OverlaySegmentString[]

    pixel_index = OverlayHotPixelIndex(precision_model, T)
    overlay_snapround_add_intersection_pixels!(pixel_index, snap_strings, T; exact)
    overlay_snapround_add_vertex_pixels!(pixel_index, snap_strings)
    snapped = overlay_snapround_compute_snaps(pixel_index, precision_model, snap_strings, T; exact)

    noded = OverlaySegmentString[]
    for snap_string in snapped
        append!(noded, overlay_snap_noded_substrings(snap_string, T))
    end
    overlay_validate_fully_noded!(noded, T; exact, precision_model = nothing)
    return noded
end

function overlay_snapround_add_intersection_pixels!(
    pixel_index::OverlayHotPixelIndex,
    snap_strings,
    ::Type{T};
    exact,
) where {T}
    records = overlay_snap_segment_records(snap_strings, T)
    isempty(records) && return pixel_index

    nearness_tol = inv(T(pixel_index.scale)) / T(OVERLAY_SNAPROUNDING_NEARNESS_FACTOR)
    extents = getproperty.(records, :extent)
    index = NaturalIndexing.NaturalIndex(extents)
    intersection_points = Tuple{T,T}[]
    for (i, record_a) in enumerate(records)
        query_extent = overlay_expand_extent(record_a.extent, nearness_tol)
        candidate_indices = SpatialTreeInterface.query(index, query_extent)
        for j in candidate_indices
            j <= i && continue
            record_b = records[j]
            Extents.intersects(query_extent, record_b.extent) || continue
            overlay_process_snaprounding_intersections!(
                intersection_points,
                records,
                i,
                j,
                T;
                exact,
                nearness_tol,
            )
        end
    end
    overlay_hot_pixel_index_add_nodes!(pixel_index, intersection_points)
    return pixel_index
end

function overlay_snap_segment_records(snap_strings, ::Type{T}) where {T}
    records = OverlaySegmentRecord[]
    for (segment_index, snap_string) in enumerate(snap_strings)
        length(snap_string.points) < 2 && continue
        for edge_index in 1:(length(snap_string.points) - 1)
            p1, p2 = snap_string.points[edge_index], snap_string.points[edge_index + 1]
            p1 == p2 && continue
            edge = (p1, p2)
            push!(
                records,
                OverlaySegmentRecord(
                    snap_string,
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

function overlay_expand_extent(extent::Extents.Extent, amount)
    return Extents.Extent(
        X = (extent.X[1] - amount, extent.X[2] + amount),
        Y = (extent.Y[1] - amount, extent.Y[2] + amount),
    )
end

function overlay_process_snaprounding_intersections!(
    intersection_points,
    records,
    i::Integer,
    j::Integer,
    ::Type{T};
    exact,
    nearness_tol,
) where {T}
    record_a = records[i]
    record_b = records[j]
    if record_a.segment === record_b.segment && record_a.edge_index == record_b.edge_index
        return intersection_points
    end

    intersection = ng_segment_intersection(record_a.edge, record_b.edge, T; exact, precision_model = nothing)
    if overlay_snap_has_interior_intersection(intersection)
        for point in ng_intersection_points(intersection)
            point = _tuple_point(point, T)
            push!(intersection_points, point)
            overlay_snap_add_intersection!(record_a.segment, point, record_a.edge_index)
            overlay_snap_add_intersection!(record_b.segment, point, record_b.edge_index)
        end
        return intersection_points
    end

    overlay_snap_process_near_vertex!(intersection_points, record_a.edge[1], record_b, nearness_tol, T)
    overlay_snap_process_near_vertex!(intersection_points, record_a.edge[2], record_b, nearness_tol, T)
    overlay_snap_process_near_vertex!(intersection_points, record_b.edge[1], record_a, nearness_tol, T)
    overlay_snap_process_near_vertex!(intersection_points, record_b.edge[2], record_a, nearness_tol, T)
    return intersection_points
end

overlay_snap_has_interior_intersection(intersection::NGSegmentIntersection) =
    intersection.orientation == line_cross || intersection.orientation == line_over

function overlay_snap_process_near_vertex!(
    intersection_points,
    point,
    record,
    nearness_tol,
    ::Type{T},
) where {T}
    p0, p1 = record.edge
    overlay_point_distance(point, p0, T) < nearness_tol && return intersection_points
    overlay_point_distance(point, p1, T) < nearness_tol && return intersection_points
    overlay_point_segment_distance(point, p0, p1, T) < nearness_tol || return intersection_points

    point = _tuple_point(point, T)
    push!(intersection_points, point)
    overlay_snap_add_intersection!(record.segment, point, record.edge_index)
    return intersection_points
end

function overlay_point_distance(a, b, ::Type{T}) where {T}
    a = _tuple_point(a, T)
    b = _tuple_point(b, T)
    return hypot(a[1] - b[1], a[2] - b[2])
end

function overlay_point_segment_distance(point, p0, p1, ::Type{T}) where {T}
    point = _tuple_point(point, T)
    p0 = _tuple_point(p0, T)
    p1 = _tuple_point(p1, T)
    dx = p1[1] - p0[1]
    dy = p1[2] - p0[2]
    len2 = dx * dx + dy * dy
    iszero(len2) && return overlay_point_distance(point, p0, T)
    frac = clamp(((point[1] - p0[1]) * dx + (point[2] - p0[2]) * dy) / len2, zero(T), one(T))
    closest = (p0[1] + frac * dx, p0[2] + frac * dy)
    return overlay_point_distance(point, closest, T)
end

function overlay_snapround_add_vertex_pixels!(pixel_index::OverlayHotPixelIndex, snap_strings)
    for snap_string in snap_strings
        for point in snap_string.points
            overlay_hot_pixel_index_add!(pixel_index, point)
        end
    end
    return pixel_index
end

function overlay_snapround_compute_snaps(
    pixel_index::OverlayHotPixelIndex,
    precision_model::FixedPrecisionModel,
    snap_strings,
    ::Type{T};
    exact,
) where {T}
    snapped = Any[]
    for snap_string in snap_strings
        snapped_string = overlay_snapround_compute_segment_snaps(
            pixel_index,
            precision_model,
            snap_string,
            T;
            exact,
        )
        isnothing(snapped_string) || push!(snapped, snapped_string)
    end
    for snap_string in snapped
        overlay_snapround_add_vertex_node_snaps!(pixel_index, snap_string)
    end
    return snapped
end

function overlay_snapround_compute_segment_snaps(
    pixel_index::OverlayHotPixelIndex,
    precision_model::FixedPrecisionModel,
    snap_string::OverlaySnapSegmentString,
    ::Type{T};
    exact,
) where {T}
    noded_points = overlay_snap_noded_coordinates(snap_string, T)
    rounded_points = overlay_round_coordinate_list(precision_model, noded_points, T)
    length(rounded_points) <= 1 && return nothing

    snapped_string = OverlaySnapSegmentString(
        rounded_points,
        snap_string.source,
        snap_string.had_repeated_coordinates,
        OverlaySnapNode{T}[],
    )
    snapped_index = 1
    for i in 1:(length(noded_points) - 1)
        current_snap = snapped_string.points[snapped_index]
        next_round = apply_ng_precision(precision_model, noded_points[i + 1], T)
        if next_round == current_snap
            continue
        end
        overlay_snapround_snap_segment!(
            pixel_index,
            snapped_string,
            noded_points[i],
            noded_points[i + 1],
            snapped_index;
            exact,
        )
        snapped_index += 1
    end
    return snapped_string
end

function overlay_round_coordinate_list(precision_model::FixedPrecisionModel, points, ::Type{T}) where {T}
    rounded = Tuple{T,T}[]
    for point in points
        rounded_point = apply_ng_precision(precision_model, point, T)
        if isempty(rounded) || rounded_point != last(rounded)
            push!(rounded, rounded_point)
        end
    end
    return rounded
end

function overlay_snapround_snap_segment!(
    pixel_index::OverlayHotPixelIndex,
    snap_string::OverlaySnapSegmentString,
    p0,
    p1,
    segment_index::Integer;
    exact,
)
    overlay_hot_pixel_index_query(pixel_index, p0, p1) do hot_pixel
        if !hot_pixel.is_node
            (overlay_hot_pixel_intersects(hot_pixel, p0) || overlay_hot_pixel_intersects(hot_pixel, p1)) && return
        end
        if overlay_hot_pixel_intersects(hot_pixel, p0, p1; exact)
            overlay_snap_add_intersection!(snap_string, hot_pixel.coordinate, segment_index)
            hot_pixel.is_node = true
        end
    end
    return snap_string
end

function overlay_snapround_add_vertex_node_snaps!(
    pixel_index::OverlayHotPixelIndex,
    snap_string::OverlaySnapSegmentString,
)
    length(snap_string.points) <= 2 && return snap_string
    for coordinate_index in 2:(length(snap_string.points) - 1)
        point = snap_string.points[coordinate_index]
        overlay_hot_pixel_index_query(pixel_index, point, point) do hot_pixel
            if hot_pixel.is_node && hot_pixel.coordinate == point
                overlay_snap_add_intersection!(snap_string, point, coordinate_index)
            end
        end
    end
    return snap_string
end

function overlay_snap_add_intersection!(
    snap_string::OverlaySnapSegmentString{T},
    point,
    segment_index::Integer,
) where {T}
    point = _tuple_point(point, T)
    normalized_index = Int(segment_index)
    next_index = normalized_index + 1
    if next_index <= length(snap_string.points) && point == snap_string.points[next_index]
        normalized_index = next_index
    end
    any(node -> node.segment_index == normalized_index && node.point == point, snap_string.nodes) &&
        return snap_string
    push!(snap_string.nodes, OverlaySnapNode(point, normalized_index))
    return snap_string
end

function overlay_snap_noded_coordinates(snap_string::OverlaySnapSegmentString{T}, ::Type{T}) where {T}
    nodes = overlay_snap_sorted_nodes(snap_string; add_collapsed_nodes = false)
    coordinates = Tuple{T,T}[]
    for i in 1:(length(nodes) - 1)
        overlay_add_coordinate_list!(coordinates, overlay_snap_split_edge_points(snap_string, nodes[i], nodes[i + 1]))
    end
    return coordinates
end

function overlay_snap_noded_substrings(snap_string::OverlaySnapSegmentString{T}, ::Type{T}) where {T}
    nodes = overlay_snap_sorted_nodes(snap_string; add_collapsed_nodes = true)
    substrings = OverlaySegmentString[]
    for i in 1:(length(nodes) - 1)
        points = overlay_snap_split_edge_points(snap_string, nodes[i], nodes[i + 1])
        overlay_edge_points_collapsed(points) && continue
        push!(
            substrings,
            OverlaySegmentString(
                Tuple{T,T}[_tuple_point(point, T) for point in points],
                snap_string.source,
                snap_string.had_repeated_coordinates,
                false,
            ),
        )
    end
    return substrings
end

function overlay_snap_sorted_nodes(
    snap_string::OverlaySnapSegmentString{T};
    add_collapsed_nodes::Bool,
) where {T}
    nodes = OverlaySnapNode{T}[copy(node) for node in snap_string.nodes]
    push!(nodes, OverlaySnapNode(first(snap_string.points), 1))
    push!(nodes, OverlaySnapNode(last(snap_string.points), length(snap_string.points)))
    if add_collapsed_nodes
        overlay_snap_add_collapsed_nodes!(nodes, snap_string)
    end
    sort!(nodes; lt = (a, b) -> overlay_snap_compare_nodes(snap_string, a, b) < 0)
    return overlay_snap_unique_nodes(snap_string, nodes)
end

function overlay_snap_add_collapsed_nodes!(nodes, snap_string::OverlaySnapSegmentString{T}) where {T}
    sorted_nodes = overlay_snap_unique_nodes(
        snap_string,
        sort!(
            OverlaySnapNode{T}[copy(node) for node in nodes];
            lt = (a, b) -> overlay_snap_compare_nodes(snap_string, a, b) < 0,
        ),
    )
    for i in 1:(length(sorted_nodes) - 1)
        collapse_index = overlay_snap_inserted_collapse_index(snap_string, sorted_nodes[i], sorted_nodes[i + 1])
        isnothing(collapse_index) || push!(nodes, OverlaySnapNode(snap_string.points[collapse_index], collapse_index))
    end
    for i in 1:(length(snap_string.points) - 2)
        snap_string.points[i] == snap_string.points[i + 2] &&
            push!(nodes, OverlaySnapNode(snap_string.points[i + 1], i + 1))
    end
    return nodes
end

function overlay_snap_inserted_collapse_index(snap_string, node_a::OverlaySnapNode, node_b::OverlaySnapNode)
    node_a.point == node_b.point || return nothing
    vertices_between = node_b.segment_index - node_a.segment_index
    overlay_snap_node_is_interior(snap_string, node_b) || (vertices_between -= 1)
    vertices_between == 1 && return node_a.segment_index + 1
    return nothing
end

function overlay_snap_unique_nodes(snap_string, nodes)
    unique_nodes = typeof(nodes)()
    for node in nodes
        if isempty(unique_nodes) || overlay_snap_compare_nodes(snap_string, last(unique_nodes), node) != 0
            push!(unique_nodes, node)
        end
    end
    return unique_nodes
end

function overlay_snap_compare_nodes(snap_string, a::OverlaySnapNode, b::OverlaySnapNode)
    a.segment_index < b.segment_index && return -1
    a.segment_index > b.segment_index && return 1
    a.point == b.point && return 0

    a_interior = overlay_snap_node_is_interior(snap_string, a)
    b_interior = overlay_snap_node_is_interior(snap_string, b)
    !a_interior && return -1
    !b_interior && return 1
    return overlay_segment_point_compare(overlay_segment_octant(snap_string.points, a.segment_index), a.point, b.point)
end

function overlay_snap_node_is_interior(snap_string, node::OverlaySnapNode)
    node.segment_index <= length(snap_string.points) || return true
    return node.point != snap_string.points[node.segment_index]
end

function overlay_segment_octant(points, segment_index::Integer)
    segment_index >= length(points) && return -1
    p0, p1 = points[segment_index], points[segment_index + 1]
    dx = p1[1] - p0[1]
    dy = p1[2] - p0[2]
    dx == 0 && dy == 0 && return 0
    abs_dx = abs(dx)
    abs_dy = abs(dy)
    if dx >= 0
        return dy >= 0 ? (abs_dx >= abs_dy ? 0 : 1) : (abs_dx >= abs_dy ? 7 : 6)
    end
    return dy >= 0 ? (abs_dx >= abs_dy ? 3 : 2) : (abs_dx >= abs_dy ? 4 : 5)
end

function overlay_segment_point_compare(octant::Integer, p0, p1)
    p0 == p1 && return 0
    x_sign = overlay_relative_sign(p0[1], p1[1])
    y_sign = overlay_relative_sign(p0[2], p1[2])
    if octant == 0
        return overlay_compare_signs(x_sign, y_sign)
    elseif octant == 1
        return overlay_compare_signs(y_sign, x_sign)
    elseif octant == 2
        return overlay_compare_signs(y_sign, -x_sign)
    elseif octant == 3
        return overlay_compare_signs(-x_sign, y_sign)
    elseif octant == 4
        return overlay_compare_signs(-x_sign, -y_sign)
    elseif octant == 5
        return overlay_compare_signs(-y_sign, -x_sign)
    elseif octant == 6
        return overlay_compare_signs(-y_sign, x_sign)
    elseif octant == 7
        return overlay_compare_signs(x_sign, -y_sign)
    end
    return 0
end

overlay_relative_sign(a, b) = a < b ? -1 : (a > b ? 1 : 0)

function overlay_compare_signs(primary::Integer, secondary::Integer)
    primary != 0 && return primary
    return secondary
end

function overlay_snap_split_edge_points(snap_string, node_a::OverlaySnapNode, node_b::OverlaySnapNode)
    npoints = node_b.segment_index - node_a.segment_index + 2
    npoints == 2 && return [node_a.point, node_b.point]

    last_segment_start = snap_string.points[node_b.segment_index]
    use_end_node = overlay_snap_node_is_interior(snap_string, node_b) || node_b.point != last_segment_start
    points = Any[node_a.point]
    for point_index in (node_a.segment_index + 1):node_b.segment_index
        push!(points, snap_string.points[point_index])
    end
    use_end_node && push!(points, node_b.point)
    return points
end

function overlay_add_coordinate_list!(coordinates, points)
    for point in points
        if isempty(coordinates) || point != last(coordinates)
            push!(coordinates, point)
        end
    end
    return coordinates
end

function overlay_edge_points_collapsed(points)
    length(points) < 2 && return true
    points[1] == points[2] && return true
    length(points) > 2 && points[end] == points[end - 1] && return true
    return false
end

function overlay_hot_pixel_index_add_nodes!(index::OverlayHotPixelIndex, points)
    for point in points
        overlay_hot_pixel_index_add!(index, point).is_node = true
    end
    return index
end

function overlay_hot_pixel_index_add!(index::OverlayHotPixelIndex{T}, point) where {T}
    rounded_point = apply_ng_precision(index.precision_model, point, T)
    pixel_index = get(index.pixel_indices, rounded_point, nothing)
    if !isnothing(pixel_index)
        pixel = index.pixels[pixel_index]
        pixel.is_node = true
        return pixel
    end

    pixel = OverlayHotPixel(rounded_point, index.scale)
    push!(index.pixels, pixel)
    index.pixel_indices[rounded_point] = length(index.pixels)
    return pixel
end

function overlay_hot_pixel_index_query(visitor, index::OverlayHotPixelIndex, p0, p1)
    width = inv(index.scale)
    extent = overlay_expand_extent(ng_segment_extent((p0, p1), typeof(index.scale)), width)
    for pixel in index.pixels
        overlay_point_in_extent(pixel.coordinate, extent) && visitor(pixel)
    end
    return nothing
end

function overlay_point_in_extent(point, extent::Extents.Extent)
    return extent.X[1] <= point[1] <= extent.X[2] &&
        extent.Y[1] <= point[2] <= extent.Y[2]
end

function overlay_hot_pixel_intersects(pixel::OverlayHotPixel{T}, point) where {T}
    x = T(GI.x(point)) * pixel.scale
    y = T(GI.y(point)) * pixel.scale
    tolerance = one(T) / 2
    x >= pixel.hpx + tolerance && return false
    x < pixel.hpx - tolerance && return false
    y >= pixel.hpy + tolerance && return false
    y < pixel.hpy - tolerance && return false
    return true
end

function overlay_hot_pixel_intersects(pixel::OverlayHotPixel{T}, p0, p1; exact = True()) where {T}
    return overlay_hot_pixel_intersects_scaled(
        pixel,
        T(GI.x(p0)) * pixel.scale,
        T(GI.y(p0)) * pixel.scale,
        T(GI.x(p1)) * pixel.scale,
        T(GI.y(p1)) * pixel.scale;
        exact,
    )
end

function overlay_hot_pixel_intersects_scaled(
    pixel::OverlayHotPixel{T},
    p0x,
    p0y,
    p1x,
    p1y;
    exact,
) where {T}
    px, py, qx, qy = p0x, p0y, p1x, p1y
    if px > qx
        px, py, qx, qy = p1x, p1y, p0x, p0y
    end

    tolerance = one(T) / 2
    maxx = pixel.hpx + tolerance
    min(px, qx) >= maxx && return false
    minx = pixel.hpx - tolerance
    max(px, qx) < minx && return false
    maxy = pixel.hpy + tolerance
    min(py, qy) >= maxy && return false
    miny = pixel.hpy - tolerance
    max(py, qy) < miny && return false

    px == qx && return true
    py == qy && return true

    orient_ul = overlay_hot_pixel_orientation(px, py, qx, qy, minx, maxy; exact)
    if orient_ul == 0
        py < qy && return false
        return true
    end

    orient_ur = overlay_hot_pixel_orientation(px, py, qx, qy, maxx, maxy; exact)
    if orient_ur == 0
        py > qy && return false
        return true
    end
    orient_ul != orient_ur && return true

    orient_ll = overlay_hot_pixel_orientation(px, py, qx, qy, minx, miny; exact)
    orient_ll == 0 && return true
    orient_ll != orient_ul && return true

    orient_lr = overlay_hot_pixel_orientation(px, py, qx, qy, maxx, miny; exact)
    if orient_lr == 0
        py < qy && return false
        return true
    end
    orient_ll != orient_lr && return true
    orient_lr != orient_ur && return true
    return false
end

function overlay_hot_pixel_orientation(px, py, qx, qy, x, y; exact)
    # JTS evaluates exact corner hits against computed noding coordinates.  The
    # local line intersector can leave those coordinates a few ulps off the
    # corner, so normalize near-zero determinants back to the intended corner hit.
    det = (qx - px) * (y - py) - (qy - py) * (x - px)
    scale = max(abs(qx - px), abs(qy - py), abs(x - px), abs(y - py), 1)
    abs(det) <= 16 * eps(typeof(det)) * scale * scale && return 0
    orientation = ng_orientation((px, py), (qx, qy), (x, y); exact)
    return orientation < 0 ? -1 : (orientation > 0 ? 1 : 0)
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
    return overlay_compute_edge_overlay(alg, op, input_a, input_b, T; target)
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

function overlay_compute_edge_overlay(
    alg::OverlayNG,
    op::OverlayOpCode,
    input_a::OverlayInputGeometry,
    input_b::OverlayInputGeometry,
    ::Type{T};
    target = nothing,
) where {T}
    noded = overlay_node_segment_strings(alg, input_a, input_b, T)
    graph = overlay_graph(overlay_merge_edges(noded))
    overlay_compute_labelling!(graph, input_a, input_b)
    overlay_mark_result_edges!(
        graph,
        op;
        strict = alg.strict,
        input_area_side = overlay_input_area_side(input_a, input_b),
    )

    polygons = overlay_extract_result_polygons(graph)
    results = copy(polygons)
    if !alg.area_result_only
        lines = Any[]
        if overlay_allows_result_lines(op, alg.strict, !isempty(polygons))
            lines = overlay_extract_result_lines(graph)
            append!(results, lines)
        end

        has_result_components = !isempty(polygons) || !isempty(lines)
        if op == overlay_intersection &&
           overlay_allows_result_intersection_points(alg.strict, has_result_components)
            append!(results, overlay_extract_intersection_points(graph; strict = alg.strict))
        end
    end
    return overlay_filter_results(alg, target, results)
end

overlay_allows_result_lines(op::OverlayOpCode, strict::Bool, has_result_area::Bool) =
    !has_result_area || !strict || op == overlay_symdifference || op == overlay_union

overlay_allows_result_intersection_points(strict::Bool, has_result_components::Bool) =
    !has_result_components || !strict

function overlay_input_area_side(a_input::OverlayInputGeometry, b_input::OverlayInputGeometry)
    a_input.dimension == dim_area && return input_a
    b_input.dimension == dim_area && return input_b
    return nothing
end

function overlay_compute_point_point(
    alg::OverlayNG,
    op::OverlayOpCode,
    input_a::OverlayInputGeometry,
    input_b::OverlayInputGeometry,
    ::Type{T};
    target = nothing,
) where {T}
    points_a = overlay_unique_points(input_a, T; precision_model = alg.precision_model)
    points_b = overlay_unique_points(input_b, T; precision_model = alg.precision_model)
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
    if op == overlay_intersection
        covered_points, _ = overlay_partition_points(
            point_input,
            nonpoint_input,
            T;
            precision_model = alg.precision_model,
        )
        return overlay_filter_results(alg, target, overlay_point_geometries(covered_points))
    elseif op == overlay_union || op == overlay_symdifference
        nonpoint_geometries = overlay_point_dispatch_nonpoint_geometries(alg, nonpoint_input, T)
        _, exterior_points = overlay_partition_points_against_geometry(
            point_input,
            overlay_components_geometry(nonpoint_geometries),
            T;
            precision_model = alg.precision_model,
        )
        return overlay_filter_results(
            alg,
            target,
            Any[nonpoint_geometries..., overlay_point_geometries(exterior_points)...],
        )
    elseif op == overlay_difference
        if point_is_a
            _, exterior_points = overlay_partition_points(
                point_input,
                nonpoint_input,
                T;
                precision_model = alg.precision_model,
            )
            return overlay_filter_results(alg, target, overlay_point_geometries(exterior_points))
        else
            return overlay_filter_results(
                alg,
                target,
                overlay_point_dispatch_nonpoint_geometries(alg, nonpoint_input, T),
            )
        end
    end
    throw(ArgumentError("Unknown OverlayNG operation code: $op"))
end

function overlay_point_dispatch_nonpoint_geometries(
    alg::OverlayNG,
    input::OverlayInputGeometry,
    ::Type{T},
) where {T}
    if overlay_has_fixed_precision(alg.precision_model)
        return overlay_precision_geometries(input.geom, alg.precision_model, T)
    end

    input.dimension == dim_line || return Any[input.geom]

    raw_segments = overlay_segment_strings(input, T; input_side = input_a)
    raw_keys = overlay_record_segment_keys(overlay_segment_records(raw_segments, T), T)
    noded_segments = overlay_node_segment_strings(alg, raw_segments, T)
    noded_keys = overlay_segment_string_keys(noded_segments, T)
    raw_keys == noded_keys && return Any[input.geom]

    geometries = Any[]
    seen = Set{Any}()
    for segment in noded_segments
        length(segment.points) < 2 && continue
        edge = (first(segment.points), last(segment.points))
        key = overlay_segment_key(edge, T)
        key in seen && continue
        push!(seen, key)
        push!(geometries, GI.LineString(copy(segment.points)))
    end
    return geometries
end

overlay_has_fixed_precision(model) = false
overlay_has_fixed_precision(::FixedPrecisionModel) = true

function overlay_record_segment_keys(records, ::Type{T}) where {T}
    keys = Set{Any}()
    for record in records
        push!(keys, overlay_segment_key(record.edge, T))
    end
    return keys
end

function overlay_segment_string_keys(segments, ::Type{T}) where {T}
    keys = Set{Any}()
    for segment in segments
        length(segment.points) < 2 && continue
        push!(keys, overlay_segment_key((first(segment.points), last(segment.points)), T))
    end
    return keys
end

function overlay_segment_key((p1, p2), ::Type{T}) where {T}
    edge = (_tuple_point(p1, T), _tuple_point(p2, T))
    return min(edge, (edge[2], edge[1]))
end

function overlay_partition_points(
    point_input::OverlayInputGeometry,
    nonpoint_input::OverlayInputGeometry,
    ::Type{T};
    precision_model = nothing,
) where {T}
    return overlay_partition_points_with_locator(
        point_input,
        nonpoint_input.locator,
        T;
        precision_model,
    )
end

function overlay_partition_points_against_geometry(
    point_input::OverlayInputGeometry,
    nonpoint_geom,
    ::Type{T};
    precision_model = nothing,
) where {T}
    return overlay_partition_points_with_locator(
        point_input,
        RelatePointLocator(nonpoint_geom),
        T;
        precision_model,
    )
end

function overlay_partition_points_with_locator(
    point_input::OverlayInputGeometry,
    locator::RelatePointLocator,
    ::Type{T};
    precision_model = nothing,
) where {T}
    covered_points = Any[]
    exterior_points = Any[]
    for point in overlay_unique_points(point_input, T; precision_model)
        target_location = relate_locate_with_dim(locator, point)
        if target_location.location == loc_exterior
            push!(exterior_points, point)
        else
            push!(covered_points, point)
        end
    end
    return covered_points, exterior_points
end

overlay_components_geometry(geometries) =
    length(geometries) == 1 ? only(geometries) : GI.GeometryCollection(geometries)

function overlay_unique_points(
    input::OverlayInputGeometry,
    ::Type{T};
    precision_model = nothing,
) where {T}
    seen = Set{Any}()
    points = Any[]
    for extracted in extract_ng_points(input.geom, T)
        point = apply_ng_precision(precision_model, extracted.point, T)
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

function overlay_precision_geometries(geom, precision_model, ::Type{T}) where {T}
    return overlay_precision_geometries(GI.trait(geom), geom, precision_model, T)
end

function overlay_precision_geometries(::GI.PointTrait, geom, precision_model, ::Type{T}) where {T}
    GI.isempty(geom) && return Any[]
    point = apply_ng_precision(precision_model, tuples(geom), T)
    return Any[GI.Point(point[1], point[2])]
end

function overlay_precision_geometries(::GI.LineStringTrait, geom, precision_model, ::Type{T}) where {T}
    points = overlay_precision_curve_points(geom, precision_model, T)
    length(points) < 2 && return Any[]
    return Any[GI.LineString(points)]
end

function overlay_precision_geometries(::GI.LinearRingTrait, geom, precision_model, ::Type{T}) where {T}
    points = overlay_precision_curve_points(geom, precision_model, T)
    length(points) < 2 && return Any[]
    return Any[GI.LineString(points)]
end

function overlay_precision_geometries(::GI.PolygonTrait, geom, precision_model, ::Type{T}) where {T}
    rings = Any[]
    for ring in GI.getring(geom)
        points = overlay_precision_ring_points(ring, precision_model, T)
        length(points) >= 4 || continue
        push!(rings, points)
    end
    isempty(rings) && return Any[]
    return Any[GI.Polygon(rings)]
end

function overlay_precision_geometries(trait::GI.AbstractGeometryTrait, geom, precision_model, ::Type{T}) where {T}
    if trait isa GI.MultiPointTrait ||
       trait isa GI.MultiLineStringTrait ||
       trait isa GI.MultiPolygonTrait ||
       trait isa GI.GeometryCollectionTrait
        return reduce(
            vcat,
            (overlay_precision_geometries(child, precision_model, T) for child in GI.getgeom(geom));
            init = Any[],
        )
    end
    return Any[geom]
end

function overlay_precision_curve_points(curve, precision_model, ::Type{T}) where {T}
    points = [apply_ng_precision(precision_model, tuples(point), T) for point in GI.getpoint(curve)]
    cleaned, _ = remove_repeated_ng_points(points)
    return cleaned
end

function overlay_precision_ring_points(ring, precision_model, ::Type{T}) where {T}
    points = overlay_precision_curve_points(ring, precision_model, T)
    isempty(points) && return points
    first(points) == last(points) || push!(points, first(points))
    cleaned, _ = remove_repeated_ng_points(points)
    first(cleaned) == last(cleaned) || push!(cleaned, first(cleaned))
    return cleaned
end

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
