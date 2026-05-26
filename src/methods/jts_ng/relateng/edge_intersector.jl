# # RelateNG edge intersection events

"""
    relate_compute_edge_intersections!(computer, [T])

Intersect A/B segment strings and feed node-section events to the topology computer.
"""
function relate_compute_edge_intersections!(
    computer::RelateTopologyComputer,
    ::Type{T} = Float64;
    exact = True(),
) where {T}
    _relate_extents_intersect(computer.geom_a.extent, computer.geom_b.extent) || return computer

    segments_a = relate_segment_strings(computer.geom_a, T; input_side = input_a)
    segments_b = relate_segment_strings(computer.geom_b, T; input_side = input_b)
    for segment_a in segments_a
        relate_process_segment_pairs!(computer, segment_a, segments_b, T; exact)
        relate_is_result_known(computer) && return computer
    end
    return computer
end

_relate_extents_intersect(::Nothing, extent_b) = true
_relate_extents_intersect(extent_a, ::Nothing) = true
_relate_extents_intersect(extent_a, extent_b) = Extents.intersects(extent_a, extent_b)

function relate_process_segment_pairs!(
    computer::RelateTopologyComputer,
    segment_a::NGSegmentString,
    segments_b,
    ::Type{T};
    exact,
) where {T}
    length(segment_a.points) < 2 && return computer
    for segment_b in segments_b
        length(segment_b.points) < 2 && continue
        relate_process_segment_pairs!(computer, segment_a, segment_b, T; exact)
        relate_is_result_known(computer) && return computer
    end
    return computer
end

function relate_process_segment_pairs!(
    computer::RelateTopologyComputer,
    segment_a::NGSegmentString,
    segment_b::NGSegmentString,
    ::Type{T};
    exact,
) where {T}
    for index_a in 1:(length(segment_a.points) - 1)
        edge_a = (segment_a.points[index_a], segment_a.points[index_a + 1])
        for index_b in 1:(length(segment_b.points) - 1)
            edge_b = (segment_b.points[index_b], segment_b.points[index_b + 1])
            ng_segments_maybe_intersect(edge_a, edge_b, T) || continue
            relate_process_segment_intersection!(
                computer,
                segment_a,
                index_a,
                edge_a,
                segment_b,
                index_b,
                edge_b,
                T;
                exact,
            )
            relate_is_result_known(computer) && return computer
        end
    end
    return computer
end

function relate_process_segment_intersection!(
    computer::RelateTopologyComputer,
    segment_a::NGSegmentString,
    index_a::Integer,
    edge_a,
    segment_b::NGSegmentString,
    index_b::Integer,
    edge_b,
    ::Type{T};
    exact,
) where {T}
    intersection = ng_segment_intersection(edge_a, edge_b, T; exact)
    ng_has_intersection(intersection) || return computer

    if intersection.orientation == line_over
        relate_update_dimension!(computer, loc_interior, loc_interior, dim_line)
    end

    is_proper = intersection.orientation == line_cross
    for point in ng_intersection_points(intersection)
        if is_proper ||
                (ng_is_containing_segment(segment_a, index_a, point) &&
                 ng_is_containing_segment(segment_b, index_b, point))
            section_a = RelateNodeSection(segment_a, index_a, point)
            section_b = RelateNodeSection(segment_b, index_b, point)
            relate_add_intersection!(computer, section_a, section_b)
        end
        relate_is_result_known(computer) && return computer
    end
    return computer
end
