# # RelateNG local node topology

const RELATE_DIM_UNKNOWN = nothing
const RELATE_LOC_UNKNOWN = nothing
const RelateMaybeDimension = Union{Nothing,TopologicalDimension}
const RelateMaybeLocation = Union{Nothing,TopologicalLocation}

"""
    RelateNodeSections(point)

Node-section accumulator for all incident RelateNG edge neighbourhoods at a point.
"""
struct RelateNodeSections{T}
    point::Tuple{T,T}
    sections::Vector{Any}
end

RelateNodeSections(point::Tuple{T,T}) where {T} = RelateNodeSections(point, Any[])

function relate_add_node_section!(node_sections::RelateNodeSections, section::RelateNodeSection)
    push!(node_sections.sections, section)
    return node_sections
end

Base.length(node_sections::RelateNodeSections) = length(node_sections.sections)
Base.iterate(node_sections::RelateNodeSections, state...) = iterate(node_sections.sections, state...)

function relate_has_interaction_ab(node_sections::RelateNodeSections)
    has_a = false
    has_b = false
    for section in node_sections.sections
        has_a |= section.input_side == input_a
        has_b |= section.input_side == input_b
        has_a && has_b && return true
    end
    return false
end

function relate_get_polygonal(node_sections::RelateNodeSections, input_side::NGInputSide)
    for section in node_sections.sections
        section.input_side == input_side || continue
        isnothing(section.parent_polygonal) || return section.parent_polygonal
    end
    return nothing
end

function relate_prepare_sections(node_sections::RelateNodeSections)
    return sort(
        copy(node_sections.sections);
        by = section -> (
            Int(section.input_side),
            dimension_value(section.dimension),
            section.element_id,
            section.ring_id,
            _relate_point_sort_key(section.previous_vertex),
            _relate_point_sort_key(section.next_vertex),
        ),
    )
end

_relate_point_sort_key(::Nothing) = (-Inf, -Inf)
_relate_point_sort_key(point::Tuple) = point

relate_is_area(section::RelateNodeSection) = section.dimension == dim_area
relate_is_shell(section::RelateNodeSection) = section.ring_id == 0
relate_is_same_polygon(a::RelateNodeSection, b::RelateNodeSection) =
    a.input_side == b.input_side && a.element_id == b.element_id

"""
    RelateNode

Angularly ordered edge star around a topology node.
"""
mutable struct RelateNode{T}
    point::Tuple{T,T}
    edges::Vector{Any}
end

RelateNode(point::Tuple{T,T}) where {T} = RelateNode(point, Any[])

"""
    RelateEdge

Directed edge label around a `RelateNode`, carrying A/B side and on locations.
"""
mutable struct RelateEdge{T}
    node_point::Tuple{T,T}
    direction_point::Tuple{T,T}
    dim_a::RelateMaybeDimension
    left_a::RelateMaybeLocation
    on_a::RelateMaybeLocation
    right_a::RelateMaybeLocation
    dim_b::RelateMaybeDimension
    left_b::RelateMaybeLocation
    on_b::RelateMaybeLocation
    right_b::RelateMaybeLocation
end

function RelateEdge(node_point::Tuple{T,T}, direction_point::Tuple{T,T}) where {T}
    return RelateEdge(
        node_point,
        direction_point,
        RELATE_DIM_UNKNOWN,
        RELATE_LOC_UNKNOWN,
        RELATE_LOC_UNKNOWN,
        RELATE_LOC_UNKNOWN,
        RELATE_DIM_UNKNOWN,
        RELATE_LOC_UNKNOWN,
        RELATE_LOC_UNKNOWN,
        RELATE_LOC_UNKNOWN,
    )
end

function relate_create_edge(
    node_point,
    direction_point,
    input_side::NGInputSide,
    dimension::TopologicalDimension,
    is_forward::Bool,
)
    edge = RelateEdge(node_point, direction_point)
    relate_merge_edge!(edge, input_side, dimension, is_forward)
    return edge
end

function relate_edge_angle(edge::RelateEdge)
    return _relate_angle(edge.node_point, edge.direction_point)
end

function _relate_angle(origin, point)
    angle = atan(point[2] - origin[2], point[1] - origin[1])
    return angle < 0 ? angle + 2pi : angle
end

function relate_edge_location(edge::RelateEdge, input_side::NGInputSide, side::SidePosition)
    if input_side == input_a
        side == side_left && return edge.left_a
        side == side_on && return edge.on_a
        side == side_right && return edge.right_a
    else
        side == side_left && return edge.left_b
        side == side_on && return edge.on_b
        side == side_right && return edge.right_b
    end
    error("Unknown side position: $side")
end

function relate_set_edge_location!(
    edge::RelateEdge,
    input_side::NGInputSide,
    side::SidePosition,
    location::TopologicalLocation,
)
    if input_side == input_a
        side == side_left && (edge.left_a = location; return edge)
        side == side_on && (edge.on_a = location; return edge)
        side == side_right && (edge.right_a = location; return edge)
    else
        side == side_left && (edge.left_b = location; return edge)
        side == side_on && (edge.on_b = location; return edge)
        side == side_right && (edge.right_b = location; return edge)
    end
    error("Unknown side position: $side")
end

function relate_edge_dimension(edge::RelateEdge, input_side::NGInputSide)
    return input_side == input_a ? edge.dim_a : edge.dim_b
end

function relate_set_edge_dimension!(
    edge::RelateEdge,
    input_side::NGInputSide,
    dimension::TopologicalDimension,
)
    if input_side == input_a
        edge.dim_a = dimension
    else
        edge.dim_b = dimension
    end
    return edge
end

relate_is_known(edge::RelateEdge, input_side::NGInputSide) =
    !isnothing(relate_edge_dimension(edge, input_side))

relate_is_known(edge::RelateEdge, input_side::NGInputSide, side::SidePosition) =
    !isnothing(relate_edge_location(edge, input_side, side))

function relate_set_line_locations!(edge::RelateEdge, input_side::NGInputSide)
    relate_set_edge_dimension!(edge, input_side, dim_line)
    relate_set_edge_location!(edge, input_side, side_left, loc_exterior)
    relate_set_edge_location!(edge, input_side, side_on, loc_interior)
    relate_set_edge_location!(edge, input_side, side_right, loc_exterior)
    return edge
end

function relate_set_area_locations!(
    edge::RelateEdge,
    input_side::NGInputSide,
    is_forward::Bool,
)
    relate_set_edge_dimension!(edge, input_side, dim_area)
    relate_set_edge_location!(edge, input_side, side_left, is_forward ? loc_exterior : loc_interior)
    relate_set_edge_location!(edge, input_side, side_on, loc_boundary)
    relate_set_edge_location!(edge, input_side, side_right, is_forward ? loc_interior : loc_exterior)
    return edge
end

function relate_merge_edge!(
    edge::RelateEdge,
    input_side::NGInputSide,
    dimension::TopologicalDimension,
    is_forward::Bool = false,
)
    if !relate_is_known(edge, input_side)
        dimension == dim_area ?
            relate_set_area_locations!(edge, input_side, is_forward) :
            relate_set_line_locations!(edge, input_side)
        return edge
    end

    edge_location = dimension == dim_area ? loc_boundary : loc_interior
    left_location = dimension == dim_area ? (is_forward ? loc_exterior : loc_interior) : loc_exterior
    right_location = dimension == dim_area ? (is_forward ? loc_interior : loc_exterior) : loc_exterior

    if dimension == dim_area && relate_edge_dimension(edge, input_side) == dim_line
        relate_set_edge_dimension!(edge, input_side, dim_area)
        relate_set_edge_location!(edge, input_side, side_on, loc_boundary)
    else
        relate_merge_on_location!(edge, input_side, edge_location)
    end
    relate_merge_side_location!(edge, input_side, side_left, left_location)
    relate_merge_side_location!(edge, input_side, side_right, right_location)
    return edge
end

function relate_merge_on_location!(
    edge::RelateEdge,
    input_side::NGInputSide,
    location::TopologicalLocation,
)
    current = relate_edge_location(edge, input_side, side_on)
    isnothing(current) && relate_set_edge_location!(edge, input_side, side_on, location)
    return edge
end

function relate_merge_side_location!(
    edge::RelateEdge,
    input_side::NGInputSide,
    side::SidePosition,
    location::TopologicalLocation,
)
    current = relate_edge_location(edge, input_side, side)
    if current != loc_interior
        relate_set_edge_location!(edge, input_side, side, location)
    end
    return edge
end

function relate_set_area_interior!(edge::RelateEdge, input_side::NGInputSide)
    relate_set_edge_location!(edge, input_side, side_left, loc_interior)
    relate_set_edge_location!(edge, input_side, side_on, loc_interior)
    relate_set_edge_location!(edge, input_side, side_right, loc_interior)
    return edge
end

function relate_set_unknown_locations!(
    edge::RelateEdge,
    input_side::NGInputSide,
    location::TopologicalLocation,
)
    for side in (side_left, side_on, side_right)
        relate_is_known(edge, input_side, side) ||
            relate_set_edge_location!(edge, input_side, side, location)
    end
    return edge
end

function relate_add_edges!(node::RelateNode, section::RelateNodeSection)
    if section.dimension == dim_line
        relate_add_line_edge!(node, section.input_side, section.previous_vertex)
        relate_add_line_edge!(node, section.input_side, section.next_vertex)
    elseif section.dimension == dim_area
        edge0 = relate_add_area_edge!(node, section.input_side, section.previous_vertex, false)
        edge1 = relate_add_area_edge!(node, section.input_side, section.next_vertex, true)
        if !isnothing(edge0) && !isnothing(edge1)
            index0 = findfirst(==(edge0), node.edges)
            index1 = findfirst(==(edge1), node.edges)
            relate_update_edges_in_area!(node, section.input_side, index0, index1)
            relate_update_if_area_prev!(node, section.input_side, index0)
            relate_update_if_area_next!(node, section.input_side, index1)
        end
    end
    return node
end

function relate_add_line_edge!(node::RelateNode, input_side::NGInputSide, direction_point)
    return relate_add_edge!(node, input_side, direction_point, dim_line, false)
end

function relate_add_area_edge!(
    node::RelateNode,
    input_side::NGInputSide,
    direction_point,
    is_forward::Bool,
)
    return relate_add_edge!(node, input_side, direction_point, dim_area, is_forward)
end

function relate_add_edge!(
    node::RelateNode,
    input_side::NGInputSide,
    direction_point::Nothing,
    dimension::TopologicalDimension,
    is_forward::Bool,
)
    return nothing
end

function relate_add_edge!(
    node::RelateNode{T},
    input_side::NGInputSide,
    direction_point,
    dimension::TopologicalDimension,
    is_forward::Bool,
) where {T}
    direction = _tuple_point(direction_point, T)
    direction == node.point && return nothing

    angle = _relate_angle(node.point, direction)
    for (i, edge) in pairs(node.edges)
        edge_angle = relate_edge_angle(edge)
        if isapprox(edge_angle, angle; atol = eps(T), rtol = zero(T))
            relate_merge_edge!(edge, input_side, dimension, is_forward)
            return edge
        elseif edge_angle > angle
            edge = relate_create_edge(node.point, direction, input_side, dimension, is_forward)
            insert!(node.edges, i, edge)
            return edge
        end
    end

    edge = relate_create_edge(node.point, direction, input_side, dimension, is_forward)
    push!(node.edges, edge)
    return edge
end

function relate_update_edges_in_area!(
    node::RelateNode,
    input_side::NGInputSide,
    index_from::Integer,
    index_to::Integer,
)
    index = _relate_next_index(node.edges, index_from)
    while index != index_to
        relate_set_area_interior!(node.edges[index], input_side)
        index = _relate_next_index(node.edges, index)
    end
    return node
end

function relate_update_if_area_prev!(node::RelateNode, input_side::NGInputSide, index::Integer)
    prev_edge = node.edges[_relate_prev_index(node.edges, index)]
    if relate_edge_location(prev_edge, input_side, side_left) == loc_interior
        relate_set_area_interior!(node.edges[index], input_side)
    end
    return node
end

function relate_update_if_area_next!(node::RelateNode, input_side::NGInputSide, index::Integer)
    next_edge = node.edges[_relate_next_index(node.edges, index)]
    if relate_edge_location(next_edge, input_side, side_right) == loc_interior
        relate_set_area_interior!(node.edges[index], input_side)
    end
    return node
end

_relate_prev_index(list, index) = index > firstindex(list) ? index - 1 : lastindex(list)
_relate_next_index(list, index) = index == lastindex(list) ? firstindex(list) : index + 1

function relate_find_known_edge_index(node::RelateNode, input_side::NGInputSide)
    return findfirst(edge -> relate_is_known(edge, input_side), node.edges)
end

function relate_finish_node!(
    node::RelateNode,
    input_side::NGInputSide,
    is_area_interior::Bool,
)
    if is_area_interior
        foreach(edge -> relate_set_area_interior!(edge, input_side), node.edges)
        return node
    end

    start_index = relate_find_known_edge_index(node, input_side)
    isnothing(start_index) && return node
    current_location = relate_edge_location(node.edges[start_index], input_side, side_left)
    isnothing(current_location) && return node

    index = _relate_next_index(node.edges, start_index)
    while index != start_index
        edge = node.edges[index]
        relate_set_unknown_locations!(edge, input_side, current_location)
        current_location = relate_edge_location(edge, input_side, side_left)
        index = _relate_next_index(node.edges, index)
    end
    return node
end

function relate_finish_node!(
    node::RelateNode,
    is_area_interior_a::Bool,
    is_area_interior_b::Bool,
)
    relate_finish_node!(node, input_a, is_area_interior_a)
    relate_finish_node!(node, input_b, is_area_interior_b)
    return node
end

function relate_create_node(node_sections::RelateNodeSections)
    node = RelateNode(node_sections.point)
    sections = relate_prepare_sections(node_sections)
    i = firstindex(sections)
    while i <= lastindex(sections)
        section = sections[i]
        if relate_is_area(section) && _relate_has_multiple_polygon_sections(sections, i)
            polygon_sections, next_index = _relate_collect_polygon_sections(sections, i)
            for converted_section in relate_convert_polygon_sections(polygon_sections)
                relate_add_edges!(node, converted_section)
            end
            i = next_index
        else
            relate_add_edges!(node, section)
            i += 1
        end
    end
    return node
end

function _relate_has_multiple_polygon_sections(sections, index::Integer)
    index < lastindex(sections) || return false
    return relate_is_same_polygon(sections[index], sections[index + 1])
end

function _relate_collect_polygon_sections(sections, index::Integer)
    polygon_sections = Any[]
    first_section = sections[index]
    while index <= lastindex(sections) && relate_is_same_polygon(first_section, sections[index])
        push!(polygon_sections, sections[index])
        index += 1
    end
    return polygon_sections, index
end

"""
    relate_convert_polygon_sections(sections)

Convert same-polygon ring sections at a node from minimal-ring to maximal-ring form.
"""
function relate_convert_polygon_sections(poly_sections)
    sections = relate_extract_unique_sections(
        sort(
            collect(poly_sections);
            by = section -> _relate_section_angle(section.previous_vertex, section.node_point),
        ),
    )
    length(sections) == 1 && return sections

    shell_index = findfirst(relate_is_shell, sections)
    isnothing(shell_index) && return relate_convert_hole_sections(sections)

    converted_sections = Any[]
    next_shell_index = shell_index
    while true
        next_shell_index = relate_convert_shell_and_holes!(
            converted_sections,
            sections,
            next_shell_index,
        )
        next_shell_index == shell_index && break
    end
    return converted_sections
end

function relate_convert_shell_and_holes!(converted_sections, sections, shell_index::Integer)
    shell_section = sections[shell_index]
    in_vertex = shell_section.previous_vertex
    i = _relate_next_wrapped_index(sections, shell_index)

    while !relate_is_shell(sections[i])
        hole_section = sections[i]
        push!(
            converted_sections,
            relate_create_section(shell_section, in_vertex, hole_section.next_vertex),
        )
        in_vertex = hole_section.previous_vertex
        i = _relate_next_wrapped_index(sections, i)
    end

    push!(
        converted_sections,
        relate_create_section(shell_section, in_vertex, shell_section.next_vertex),
    )
    return i
end

function relate_convert_hole_sections(sections)
    converted_sections = Any[]
    copy_section = first(sections)
    for i in eachindex(sections)
        next_index = _relate_next_wrapped_index(sections, i)
        push!(
            converted_sections,
            relate_create_section(
                copy_section,
                sections[i].previous_vertex,
                sections[next_index].next_vertex,
            ),
        )
    end
    return converted_sections
end

function relate_create_section(section::RelateNodeSection, previous_vertex, next_vertex)
    return RelateNodeSection(
        section.input_side,
        dim_area,
        section.element_id,
        0,
        section.parent_polygonal,
        section.is_node_at_vertex,
        section.node_point,
        previous_vertex,
        next_vertex,
    )
end

function relate_extract_unique_sections(sections)
    isempty(sections) && return sections
    unique_sections = Any[first(sections)]
    last_unique = first(sections)
    for section in Iterators.drop(sections, 1)
        _relate_section_key(section) == _relate_section_key(last_unique) && continue
        push!(unique_sections, section)
        last_unique = section
    end
    return unique_sections
end

function _relate_section_key(section::RelateNodeSection)
    return (
        section.input_side,
        section.dimension,
        section.element_id,
        section.ring_id,
        section.previous_vertex,
        section.node_point,
        section.next_vertex,
    )
end

function _relate_section_angle(previous_vertex, node_point)
    isnothing(previous_vertex) && return -Inf
    return _relate_angle(node_point, previous_vertex)
end

_relate_next_wrapped_index(list, index::Integer) =
    index == lastindex(list) ? firstindex(list) : index + 1

function relate_evaluate_node_edges!(computer::RelateTopologyComputer, node::RelateNode)
    for edge in node.edges
        if relate_is_area_area(computer)
            relate_update_edge_side_dimension!(computer, edge, side_left, dim_area)
            relate_update_edge_side_dimension!(computer, edge, side_right, dim_area)
        end
        relate_update_edge_side_dimension!(computer, edge, side_on, dim_line)
        relate_is_result_known(computer) && return computer
    end
    return computer
end

function relate_update_edge_side_dimension!(
    computer::RelateTopologyComputer,
    edge::RelateEdge,
    side::SidePosition,
    dimension::TopologicalDimension,
)
    loc_a = relate_edge_location(edge, input_a, side)
    loc_b = relate_edge_location(edge, input_b, side)
    (isnothing(loc_a) || isnothing(loc_b)) && return computer
    relate_update_dimension!(computer, loc_a, loc_b, dimension)
    return computer
end
