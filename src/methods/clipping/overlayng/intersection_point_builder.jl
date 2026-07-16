# # IntersectionPointBuilder — result points (port of JTS `IntersectionPointBuilder`)
#
# Phase 2b of the OverlayNG port (design doc §3). A port of
# `operation/overlayng/IntersectionPointBuilder.java`: extracts `Point` results
# from an Intersection of non-point inputs, at nodes incident on edges of BOTH
# inputs where none of the incident edges is itself in the result (i.e. the
# inputs meet at an isolated point). `isAllowCollapseLines` is always true here
# (non-strict), so the boundary-collapse guard in `isEdgeOf` is elided.
#
# Returns the result points as coordinate tuples; the driver wraps them. Nothing
# here is exported.

function _build_points(g::OverlayGraph)
    edges = g.edges
    points = Tuple{Float64, Float64}[]
    for ne in graph_node_edges(g)
        if _is_result_point(edges, ne)
            push!(points, node_point(g.arr, he_origin(edges, ne)))
        end
    end
    return points
end

# Port of `isResultPoint`: a node incident on edges of both inputs, none of them
# in the result.
function _is_result_point(edges, node_edge::Integer)
    is_edge_of_a = false
    is_edge_of_b = false
    e = Int32(node_edge)
    while true
        oe_in_result(edges, e) && return false
        lbl = oe_label(edges, e)
        is_edge_of_a |= _is_edge_of(lbl, 0)
        is_edge_of_b |= _is_edge_of(lbl, 1)
        e = he_onext(edges, e)
        e == node_edge && break
    end
    return is_edge_of_a && is_edge_of_b
end

@inline _is_edge_of(lbl::OverlayLabel, i::Integer) = is_boundary(lbl, i) || is_line(lbl, i)
