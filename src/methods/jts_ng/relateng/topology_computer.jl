# # RelateNG topology computer

"""
    RelateNodeSection

Incident edge neighbourhood for a RelateNG node on one input segment string.
"""
struct RelateNodeSection{T,P}
    input_side::NGInputSide
    dimension::TopologicalDimension
    element_id::Int
    ring_id::Int
    parent_polygonal::P
    is_node_at_vertex::Bool
    node_point::Tuple{T,T}
    previous_vertex::Union{Nothing,Tuple{T,T}}
    next_vertex::Union{Nothing,Tuple{T,T}}
end

"""
    RelateTopologyComputer(predicate, geom_a, geom_b)

Event sink that updates a RelateNG predicate and accumulates edge node sections.
"""
mutable struct RelateTopologyComputer{P,A,B}
    predicate::P
    geom_a::A
    geom_b::B
    node_sections::Dict{Any,Any}
end

function RelateTopologyComputer(
    predicate::TopologyPredicate,
    geom_a::RelateGeometry,
    geom_b::RelateGeometry,
)
    relate_init_predicate!(predicate, geom_a, geom_b)
    computer = RelateTopologyComputer(predicate, geom_a, geom_b, Dict{Any,Any}())
    relate_init_exterior_dimensions!(computer)
    return computer
end

function RelateTopologyComputer(alg::RelateNG, predicate::TopologyPredicate, a, b)
    geom_a = RelateGeometry(
        a;
        prepared = alg.prepared,
        boundary_node_rule = alg.boundary_node_rule,
    )
    geom_b = RelateGeometry(
        b;
        prepared = false,
        boundary_node_rule = alg.boundary_node_rule,
    )
    return RelateTopologyComputer(predicate, geom_a, geom_b)
end

relate_geometry(computer::RelateTopologyComputer, input_side::NGInputSide) =
    input_side == input_a ? computer.geom_a : computer.geom_b

relate_opposite_side(input_side::NGInputSide) =
    input_side == input_a ? input_b : input_a

relate_is_area_area(computer::RelateTopologyComputer) =
    computer.geom_a.dimension == dim_area && computer.geom_b.dimension == dim_area

function relate_is_self_noding_required(computer::RelateTopologyComputer)
    require_self_noding(computer.predicate) || return false
    relate_is_self_noding_required(computer.geom_a) && return true
    relate_has_area_and_line(computer.geom_b) && return true
    return false
end

relate_is_exterior_check_required(computer::RelateTopologyComputer, input_side::NGInputSide) =
    require_exterior_check(computer.predicate, input_side)

relate_is_result_known(computer::RelateTopologyComputer) =
    predicate_is_known(computer.predicate)

relate_result(computer::RelateTopologyComputer) =
    predicate_value(computer.predicate)

function relate_finish!(computer::RelateTopologyComputer)
    relate_finish!(computer.predicate)
    return computer
end

"""
    relate_update_dimension!(computer, loc_a, loc_b, dimension)

Record that the A/B DE-9IM cell is at least `dimension`.
"""
function relate_update_dimension!(
    computer::RelateTopologyComputer,
    loc_a::TopologicalLocation,
    loc_b::TopologicalLocation,
    dimension::TopologicalDimension,
)
    relate_update_dimension!(computer.predicate, loc_a, loc_b, dimension)
    return computer
end

function relate_update_dimension!(
    computer::RelateTopologyComputer,
    input_side::NGInputSide,
    source_location::TopologicalLocation,
    target_location::TopologicalLocation,
    dimension::TopologicalDimension,
)
    if input_side == input_a
        relate_update_dimension!(computer, source_location, target_location, dimension)
    else
        relate_update_dimension!(computer, target_location, source_location, dimension)
    end
    return computer
end

"""
    relate_init_exterior_dimensions!(computer)

Seed DE-9IM facts implied by operand dimensions before detailed events run.
"""
function relate_init_exterior_dimensions!(computer::RelateTopologyComputer)
    dim_a = relate_dimension_real(computer.geom_a)
    dim_b = relate_dimension_real(computer.geom_b)

    if dim_a == dim_point && dim_b == dim_line
        relate_update_dimension!(computer, loc_exterior, loc_interior, dim_line)
    elseif dim_a == dim_line && dim_b == dim_point
        relate_update_dimension!(computer, loc_interior, loc_exterior, dim_line)
    elseif dim_a == dim_point && dim_b == dim_area
        relate_update_dimension!(computer, loc_exterior, loc_interior, dim_area)
        relate_update_dimension!(computer, loc_exterior, loc_boundary, dim_line)
    elseif dim_a == dim_area && dim_b == dim_point
        relate_update_dimension!(computer, loc_interior, loc_exterior, dim_area)
        relate_update_dimension!(computer, loc_boundary, loc_exterior, dim_line)
    elseif dim_a == dim_line && dim_b == dim_area
        relate_update_dimension!(computer, loc_exterior, loc_interior, dim_area)
    elseif dim_a == dim_area && dim_b == dim_line
        relate_update_dimension!(computer, loc_interior, loc_exterior, dim_area)
    elseif dim_a == dim_false || dim_b == dim_false
        dim_a != dim_false && relate_init_exterior_empty!(computer, input_a)
        dim_b != dim_false && relate_init_exterior_empty!(computer, input_b)
    end
    return computer
end

function relate_init_exterior_empty!(
    computer::RelateTopologyComputer,
    non_empty_side::NGInputSide,
)
    geometry = relate_geometry(computer, non_empty_side)
    dim = geometry.dimension
    if dim == dim_point
        relate_update_dimension!(computer, non_empty_side, loc_interior, loc_exterior, dim_point)
    elseif dim == dim_line
        if relate_has_boundary(geometry)
            relate_update_dimension!(computer, non_empty_side, loc_boundary, loc_exterior, dim_point)
        end
        relate_update_dimension!(computer, non_empty_side, loc_interior, loc_exterior, dim_line)
    elseif dim == dim_area
        relate_update_dimension!(computer, non_empty_side, loc_boundary, loc_exterior, dim_line)
        relate_update_dimension!(computer, non_empty_side, loc_interior, loc_exterior, dim_area)
    end
    return computer
end

"""
    relate_add_point_on_point_interior!(computer, point)

Record a point component shared by both inputs.
"""
function relate_add_point_on_point_interior!(computer::RelateTopologyComputer, point = nothing)
    relate_update_dimension!(computer, loc_interior, loc_interior, dim_point)
    return computer
end

function relate_add_point_on_point_exterior!(
    computer::RelateTopologyComputer,
    input_side::NGInputSide,
    point = nothing,
)
    relate_update_dimension!(computer, input_side, loc_interior, loc_exterior, dim_point)
    return computer
end

function relate_add_point_on_geometry!(
    computer::RelateTopologyComputer,
    point_side::NGInputSide,
    target::DimensionLocation,
    point = nothing,
)
    return relate_add_point_on_geometry!(
        computer,
        point_side,
        target.location,
        target.dimension,
        point,
    )
end

"""
    relate_add_point_on_geometry!(computer, point_side, loc_target, dim_target, point)

Record a point component location against the opposite geometry.
"""
function relate_add_point_on_geometry!(
    computer::RelateTopologyComputer,
    point_side::NGInputSide,
    target_location::TopologicalLocation,
    target_dimension::TopologicalDimension,
    point = nothing,
)
    relate_update_dimension!(computer, point_side, loc_interior, target_location, dim_point)

    target = relate_geometry(computer, relate_opposite_side(point_side))
    target.is_empty && return computer

    target_dimension == dim_point && return computer
    target_dimension == dim_line && return computer
    if target_dimension == dim_area
        relate_update_dimension!(computer, point_side, loc_exterior, loc_interior, dim_area)
        relate_update_dimension!(computer, point_side, loc_exterior, loc_boundary, dim_line)
        return computer
    end
    throw(ArgumentError("Unknown target dimension: $target_dimension"))
end

function relate_add_line_end_on_geometry!(
    computer::RelateTopologyComputer,
    line_side::NGInputSide,
    line_end_location::TopologicalLocation,
    target::DimensionLocation,
    point = nothing,
)
    return relate_add_line_end_on_geometry!(
        computer,
        line_side,
        line_end_location,
        target.location,
        target.dimension,
        point,
    )
end

"""
    relate_add_line_end_on_geometry!(computer, line_side, loc_line_end, loc_target, dim_target, point)

Record topology implied by a significant line endpoint.
"""
function relate_add_line_end_on_geometry!(
    computer::RelateTopologyComputer,
    line_side::NGInputSide,
    line_end_location::TopologicalLocation,
    target_location::TopologicalLocation,
    target_dimension::TopologicalDimension,
    point = nothing,
)
    relate_update_dimension!(computer, line_side, line_end_location, target_location, dim_point)

    target = relate_geometry(computer, relate_opposite_side(line_side))
    target.is_empty && return computer

    if target_dimension == dim_point
        return computer
    elseif target_dimension == dim_line
        relate_add_line_end_on_line!(computer, line_side, line_end_location, target_location, point)
    elseif target_dimension == dim_area
        relate_add_line_end_on_area!(computer, line_side, line_end_location, target_location, point)
    else
        throw(ArgumentError("Unknown target dimension: $target_dimension"))
    end
    return computer
end

function relate_add_line_end_on_line!(
    computer::RelateTopologyComputer,
    line_side::NGInputSide,
    line_end_location::TopologicalLocation,
    target_location::TopologicalLocation,
    point = nothing,
)
    if target_location == loc_exterior
        relate_update_dimension!(computer, line_side, loc_interior, loc_exterior, dim_line)
    end
    return computer
end

function relate_add_line_end_on_area!(
    computer::RelateTopologyComputer,
    line_side::NGInputSide,
    line_end_location::TopologicalLocation,
    area_location::TopologicalLocation,
    point = nothing,
)
    if area_location != loc_boundary
        relate_update_dimension!(computer, line_side, loc_interior, area_location, dim_line)
        relate_update_dimension!(computer, line_side, loc_exterior, area_location, dim_area)
    end
    return computer
end

function relate_add_area_vertex!(
    computer::RelateTopologyComputer,
    area_side::NGInputSide,
    area_location::TopologicalLocation,
    target::DimensionLocation,
    point = nothing,
)
    return relate_add_area_vertex!(
        computer,
        area_side,
        area_location,
        target.location,
        target.dimension,
        point,
    )
end

"""
    relate_add_area_vertex!(computer, area_side, loc_area, loc_target, dim_target, point)

Record topology implied by an area vertex against the opposite geometry.
"""
function relate_add_area_vertex!(
    computer::RelateTopologyComputer,
    area_side::NGInputSide,
    area_location::TopologicalLocation,
    target_location::TopologicalLocation,
    target_dimension::TopologicalDimension,
    point = nothing,
)
    if target_location == loc_exterior
        relate_update_dimension!(computer, area_side, loc_interior, loc_exterior, dim_area)
        if area_location == loc_boundary
            relate_update_dimension!(computer, area_side, loc_boundary, loc_exterior, dim_line)
            relate_update_dimension!(computer, area_side, loc_exterior, loc_exterior, dim_area)
        end
        return computer
    end

    if target_dimension == dim_point
        relate_add_area_vertex_on_point!(computer, area_side, area_location, point)
    elseif target_dimension == dim_line
        relate_add_area_vertex_on_line!(computer, area_side, area_location, target_location, point)
    elseif target_dimension == dim_area
        relate_add_area_vertex_on_area!(computer, area_side, area_location, target_location, point)
    else
        throw(ArgumentError("Unknown target dimension: $target_dimension"))
    end
    return computer
end

function relate_add_area_vertex_on_point!(
    computer::RelateTopologyComputer,
    area_side::NGInputSide,
    area_location::TopologicalLocation,
    point = nothing,
)
    relate_update_dimension!(computer, area_side, area_location, loc_interior, dim_point)
    relate_update_dimension!(computer, area_side, loc_interior, loc_exterior, dim_area)
    if area_location == loc_boundary
        relate_update_dimension!(computer, area_side, loc_boundary, loc_exterior, dim_line)
        relate_update_dimension!(computer, area_side, loc_exterior, loc_exterior, dim_area)
    end
    return computer
end

function relate_add_area_vertex_on_line!(
    computer::RelateTopologyComputer,
    area_side::NGInputSide,
    area_location::TopologicalLocation,
    target_location::TopologicalLocation,
    point = nothing,
)
    relate_update_dimension!(computer, area_side, area_location, target_location, dim_point)
    if area_location == loc_interior
        relate_update_dimension!(computer, area_side, loc_interior, loc_exterior, dim_area)
    end
    return computer
end

function relate_add_area_vertex_on_area!(
    computer::RelateTopologyComputer,
    area_side::NGInputSide,
    area_location::TopologicalLocation,
    target_location::TopologicalLocation,
    point = nothing,
)
    if target_location == loc_boundary
        if area_location == loc_boundary
            relate_update_dimension!(computer, area_side, loc_boundary, loc_boundary, dim_point)
        else
            relate_update_dimension!(computer, area_side, loc_interior, loc_interior, dim_area)
            relate_update_dimension!(computer, area_side, loc_interior, loc_boundary, dim_line)
            relate_update_dimension!(computer, area_side, loc_interior, loc_exterior, dim_area)
        end
    else
        relate_update_dimension!(computer, area_side, loc_interior, target_location, dim_area)
        if area_location == loc_boundary
            relate_update_dimension!(computer, area_side, loc_boundary, target_location, dim_line)
            relate_update_dimension!(computer, area_side, loc_exterior, target_location, dim_area)
        end
    end
    return computer
end

function RelateNodeSection(
    segment::NGSegmentString{T},
    segment_index::Integer,
    node_point,
) where {T}
    1 <= segment_index < length(segment.points) ||
        throw(BoundsError(segment.points, segment_index))

    point = _tuple_point(node_point, T)
    source = segment.source
    previous_vertex = _relate_previous_vertex(segment.points, segment_index, point)
    next_vertex = _relate_next_vertex(segment.points, segment_index, point)
    is_node_at_vertex =
        point == segment.points[segment_index] ||
        point == segment.points[segment_index + 1]
    return RelateNodeSection(
        source.input_side,
        source.source_dimension,
        source.element_id,
        source.ring_id,
        source.parent_polygonal,
        is_node_at_vertex,
        point,
        previous_vertex,
        next_vertex,
    )
end

function _relate_previous_vertex(points, segment_index::Integer, point)
    segment_start = points[segment_index]
    segment_start != point && return segment_start
    segment_index > 1 && return points[segment_index - 1]
    _relate_is_closed_points(points) && length(points) > 2 && return points[end - 1]
    return nothing
end

function _relate_next_vertex(points, segment_index::Integer, point)
    segment_end = points[segment_index + 1]
    segment_end != point && return segment_end
    segment_index < length(points) - 1 && return points[segment_index + 2]
    _relate_is_closed_points(points) && length(points) > 2 && return points[2]
    return nothing
end

_relate_is_closed_points(points) = length(points) > 1 && first(points) == last(points)

relate_is_same_geometry(a::RelateNodeSection, b::RelateNodeSection) =
    a.input_side == b.input_side

relate_is_area_area(a::RelateNodeSection, b::RelateNodeSection) =
    a.dimension == dim_area && b.dimension == dim_area

relate_is_proper(section::RelateNodeSection) = !section.is_node_at_vertex
relate_is_proper(a::RelateNodeSection, b::RelateNodeSection) =
    relate_is_proper(a) && relate_is_proper(b)

function relate_node_sections(computer::RelateTopologyComputer, point)
    node_sections = get(computer.node_sections, _tuple_point(point), nothing)
    isnothing(node_sections) && return Any[]
    return node_sections.sections
end

"""
    relate_add_intersection!(computer, a, b)

Record two node sections created by an edge intersection.
"""
function relate_add_intersection!(
    computer::RelateTopologyComputer,
    a::RelateNodeSection,
    b::RelateNodeSection,
)
    if !relate_is_same_geometry(a, b)
        relate_update_intersection_ab!(computer, a, b)
    end
    relate_add_node_sections!(computer, a, b)
    return computer
end

function relate_update_intersection_ab!(
    computer::RelateTopologyComputer,
    a::RelateNodeSection,
    b::RelateNodeSection,
)
    if relate_is_area_area(a, b)
        relate_update_area_area_cross!(computer, a, b)
    end
    relate_update_node_location!(computer, a, b)
    return computer
end

function relate_update_area_area_cross!(
    computer::RelateTopologyComputer,
    a::RelateNodeSection,
    b::RelateNodeSection,
)
    if relate_is_proper(a, b)
        relate_update_dimension!(computer, loc_interior, loc_interior, dim_area)
    end
    return computer
end

function relate_update_node_location!(
    computer::RelateTopologyComputer,
    a::RelateNodeSection,
    b::RelateNodeSection,
)
    section_a = a.input_side == input_a ? a : b
    section_b = a.input_side == input_b ? a : b
    point = section_a.node_point
    loc_a = relate_locate_node_with_dim(
        computer.geom_a,
        point,
        section_a.parent_polygonal,
    ).location
    loc_b = relate_locate_node_with_dim(
        computer.geom_b,
        point,
        section_b.parent_polygonal,
    ).location
    relate_update_dimension!(computer, loc_a, loc_b, dim_point)
    return computer
end

function relate_add_node_sections!(
    computer::RelateTopologyComputer,
    a::RelateNodeSection,
    b::RelateNodeSection,
)
    node_sections = get!(computer.node_sections, a.node_point) do
        RelateNodeSections(a.node_point)
    end
    relate_add_node_section!(node_sections, a)
    relate_add_node_section!(node_sections, b)
    return computer
end

"""
    relate_evaluate_nodes!(computer)

Evaluate accumulated node sections and update side/on topology dimensions.
"""
function relate_evaluate_nodes!(computer::RelateTopologyComputer)
    for node_sections in values(computer.node_sections)
        relate_has_interaction_ab(node_sections) || continue
        node = relate_create_node(node_sections)
        point = node_sections.point
        is_area_interior_a = relate_is_node_in_area(
            computer.geom_a,
            point,
            relate_get_polygonal(node_sections, input_a),
        )
        is_area_interior_b = relate_is_node_in_area(
            computer.geom_b,
            point,
            relate_get_polygonal(node_sections, input_b),
        )
        relate_finish_node!(node, is_area_interior_a, is_area_interior_b)
        relate_evaluate_node_edges!(computer, node)
        relate_is_result_known(computer) && return computer
    end
    return computer
end
