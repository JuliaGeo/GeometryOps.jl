# # LineBuilder — result lines from the overlay graph (port of JTS `LineBuilder`)
#
# Phase 2b of the OverlayNG port (design doc §3). A port of
# `operation/overlayng/LineBuilder.java`, restricted to the raw-edge extraction
# JTS actually uses (`addResultLines`); the merged path (`addResultLinesMerged`)
# is marked NOT USED in JTS and is skipped. Strict-mode branches are dropped —
# this engine runs the original (non-strict) JTS semantics, so
# `isAllowMixedResult` and `isAllowCollapseLines` are always true.
#
# Each `OverlayEdge` is a single segment, so every result line is one two-point
# `LineString` (the raw noded-edge output; JTS emits the same edges, just as
# possibly-longer chains). Everything here is internal — nothing is exported.

# Build the result lines from the graph (port of `getLines`: `markResultLines`
# then `addResultLines`).
function _build_lines(m::Manifold, g::OverlayGraph, input, has_result_area::Bool,
        op::_OverlayOpCode; exact)
    edges = g.edges
    area_index = _input_area_index(input)
    #-- markResultLines
    for i in eachindex(edges)
        oe_in_result_either(edges, i) && continue
        if _is_result_line(oe_label(edges, i), op, has_result_area, area_index)
            oe_mark_in_result_line!(edges, i)
        end
    end
    #-- addResultLines (raw noded edges)
    lines = Vector{typeof(_edge_line(g, 1))}()
    for i in eachindex(edges)
        oe_in_result_line(edges, i) || continue
        oe_is_visited(edges, i) && continue
        push!(lines, _edge_line(g, i))
        oe_mark_visited_both!(edges, i)
    end
    return lines
end

# The `LineString` of one result edge (a single segment, in the half-edge's
# direction).
_edge_line(g::OverlayGraph, i::Integer) =
    GI.LineString([node_point(g.arr, he_origin(g.edges, i)),
                   node_point(g.arr, he_dest(g.edges, i))])

# Port of `isResultLine`. `is_allow_collapse_lines` / `is_allow_mixed_result`
# are always true here (non-strict), so those guards are elided.
function _is_result_line(lbl::OverlayLabel, op::_OverlayOpCode, has_result_area::Bool,
        area_index::Integer)
    #-- a boundary of a single geometry is only in the result as part of an area
    is_boundary_singleton(lbl) && return false
    #-- a collapse interior to its parent area (narrow gore / hole spike)
    is_interior_collapse(lbl) && return false

    if op != OVERLAY_INTERSECTION
        #-- collapsed edge in the other area's interior
        is_collapse_and_not_part_interior(lbl) && return false
        #-- a line inside the result area is subsumed by it
        if has_result_area && area_index >= 0 && is_line_in_area(lbl, area_index)
            return false
        end
    end

    #-- touching area boundaries produce a line for Intersection (mixed result)
    op == OVERLAY_INTERSECTION && is_boundary_touch(lbl) && return true

    #-- otherwise, the op boolean logic over the effective line locations
    a_loc = _effective_location(lbl, 0)
    b_loc = _effective_location(lbl, 1)
    return _is_result_of_op(op, a_loc, b_loc)
end

# Port of `effectiveLocation`: line and collapse edges report INTERIOR so the op
# logic can include them where warranted.
@inline function _effective_location(lbl::OverlayLabel, gi::Integer)
    is_collapse(lbl, gi) && return LOC_INTERIOR
    is_line(lbl, gi) && return LOC_INTERIOR
    return get_line_location(lbl, gi)
end
