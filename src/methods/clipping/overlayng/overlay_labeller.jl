# # OverlayLabeller — the five-pass edge labelling (port of JTS `OverlayLabeller`)
#
# Phase 2b of the OverlayNG port (design doc `2026-07-16-overlayng-noding-substrate.md`,
# §3). A faithful port of `operation/overlayng/OverlayLabeller.java`, operating on
# the phase-2a index-based `OverlayGraph`. The single shared `OverlayLabel` between a
# symmetric pair means every `set_location_*!` here is seen by both half-edges —
# JTS relies on this, and so does this port.
#
# Also hosts the overlay op-code enum and `_is_result_of_op` (JTS
# `OverlayNG.isResultOfOp`, the boolean core used by both the labeller and the line
# builder), and the topology-error type (JTS `TopologyException`).
#
# Everything here is internal to GeometryOps — nothing is exported.

# ## Overlay operation codes and the op boolean core

# The four set-theoretic overlay operations (port of the `OverlayNG.INTERSECTION`
# / `UNION` / `DIFFERENCE` / `SYMDIFFERENCE` op codes).
@enum _OverlayOpCode::UInt8 begin
    OVERLAY_INTERSECTION
    OVERLAY_UNION
    OVERLAY_DIFFERENCE
    OVERLAY_SYMDIFFERENCE
end

"""
    _OverlayTopologyError(msg)

A robustness/topology error raised by the overlay engine (port of JTS
`TopologyException`). Signals an inconsistency the builder could not resolve
(e.g. a side-location conflict during area propagation, or a ring that cannot
be closed).
"""
struct _OverlayTopologyError <: Exception
    msg::String
end
Base.showerror(io::IO, e::_OverlayTopologyError) = print(io, "OverlayTopologyError: ", e.msg)

# Port of `OverlayNG.isResultOfOp`: whether a point with the given per-input
# locations lies in the result of `op`. `LOC_BOUNDARY` counts as `LOC_INTERIOR`.
@inline function _is_result_of_op(op::_OverlayOpCode, loc0::Integer, loc1::Integer)
    loc0 == LOC_BOUNDARY && (loc0 = LOC_INTERIOR)
    loc1 == LOC_BOUNDARY && (loc1 = LOC_INTERIOR)
    if op == OVERLAY_INTERSECTION
        return loc0 == LOC_INTERIOR && loc1 == LOC_INTERIOR
    elseif op == OVERLAY_UNION
        return loc0 == LOC_INTERIOR || loc1 == LOC_INTERIOR
    elseif op == OVERLAY_DIFFERENCE
        return loc0 == LOC_INTERIOR && loc1 != LOC_INTERIOR
    else # OVERLAY_SYMDIFFERENCE
        return (loc0 == LOC_INTERIOR && loc1 != LOC_INTERIOR) ||
               (loc0 != LOC_INTERIOR && loc1 == LOC_INTERIOR)
    end
end

# ## computeLabelling — the five passes, in order (port of `computeLabelling`)

function _compute_labelling!(g::OverlayGraph, input)
    edges = g.edges
    for ne in graph_node_edges(g)
        _propagate_area_locations!(edges, input, ne, 0)
        _input_has_edges(input, 1) && _propagate_area_locations!(edges, input, ne, 1)
    end
    _label_connected_linear_edges!(edges, input)
    _label_collapsed_edges!(edges)
    _label_connected_linear_edges!(edges, input)
    _label_disconnected_edges!(g, input)
    return nothing
end

# ### Pass 1: area-node propagation (port of `propagateAreaLocations`)

# Scans a node's star CCW, propagating side labels for one area input to every
# edge whose location for that input is still unknown. A side-location conflict
# between two boundary edges is a topology error (port of the JTS check).
function _propagate_area_locations!(edges, input, node_edge::Integer, gi::Integer)
    _input_is_area(input, gi) || return nothing
    #-- one-edge node: nothing to propagate (dangling edge)
    he_degree(edges, node_edge) == 1 && return nothing

    e_start = _find_propagation_start_edge(edges, node_edge, gi)
    e_start == 0 && return nothing

    curr_loc = oe_get_location(edges, e_start, gi, POS_LEFT)
    e = he_onext(edges, e_start)
    while e != e_start
        label = oe_label(edges, e)
        if !is_boundary(label, gi)
            #-- non-boundary edge: its location relative to this area is now known
            set_location_line!(label, gi, curr_loc)
        else
            loc_right = oe_get_location(edges, e, gi, POS_RIGHT)
            loc_right == curr_loc ||
                throw(_OverlayTopologyError("side location conflict: arg $gi"))
            loc_left = oe_get_location(edges, e, gi, POS_LEFT)
            #-- loc_left == LOC_NONE should never happen for a boundary edge
            curr_loc = loc_left
        end
        e = he_onext(edges, e)
    end
    return nothing
end

# Port of `findPropagationStartEdge`: a boundary edge for `gi` in the node's
# star, or `0` if none.
function _find_propagation_start_edge(edges, node_edge::Integer, gi::Integer)
    e = Int32(node_edge)
    while true
        is_boundary(oe_label(edges, e), gi) && return e
        e = he_onext(edges, e)
        e == node_edge && break
    end
    return Int32(0)
end

# ### Pass 2/4: connected linear propagation (port of `labelConnectedLinearEdges`)

function _label_connected_linear_edges!(edges, input)
    _propagate_linear_locations!(edges, input, 0)
    _input_has_edges(input, 1) && _propagate_linear_locations!(edges, input, 1)
    return nothing
end

# BFS over linear (line or collapse) edges with a known location, propagating
# that location to connected unknown edges (port of `propagateLinearLocations`).
function _propagate_linear_locations!(edges, input, gi::Integer)
    stack = Int32[]
    for i in eachindex(edges)
        lbl = oe_label(edges, i)
        if is_linear(lbl, gi) && !is_line_location_unknown(lbl, gi)
            push!(stack, Int32(i))
        end
    end
    isempty(stack) && return nothing

    is_input_line = _input_is_line(input, gi)
    while !isempty(stack)
        e_node = pop!(stack)
        _propagate_linear_at_node!(edges, e_node, gi, is_input_line, stack)
    end
    return nothing
end

# Port of `propagateLinearLocationAtNode`. Line parents propagate EXTERIOR only.
function _propagate_linear_at_node!(edges, e_node::Integer, gi::Integer,
        is_input_line::Bool, stack::Vector{Int32})
    line_loc = get_line_location(oe_label(edges, e_node), gi)
    #-- a Line parent only propagates EXTERIOR locations
    is_input_line && line_loc != LOC_EXTERIOR && return nothing

    e = he_onext(edges, e_node)
    while e != e_node
        label = oe_label(edges, e)
        if is_line_location_unknown(label, gi)
            set_location_line!(label, gi, line_loc)
            #-- continue the traversal from the far node (don't re-add e itself)
            push!(stack, he_sym(edges, e))
        end
        e = he_onext(edges, e)
    end
    return nothing
end

# ### Pass 3: collapsed-edge ring-role labelling (port of `labelCollapsedEdges`)

function _label_collapsed_edges!(edges)
    for i in eachindex(edges)
        label = oe_label(edges, i)
        is_line_location_unknown(label, 0) && _label_collapsed_edge!(label, 0)
        is_line_location_unknown(label, 1) && _label_collapsed_edge!(label, 1)
    end
    return nothing
end

function _label_collapsed_edge!(label::OverlayLabel, gi::Integer)
    is_collapse(label, gi) || return nothing
    #-- disconnected collapsed edge: label from its parent ring role (shell/hole)
    set_location_collapse!(label, gi)
    return nothing
end

# ### Pass 5: disconnected-edge PIP labelling (port of `labelDisconnectedEdges`)

function _label_disconnected_edges!(g::OverlayGraph, input)
    edges = g.edges
    for i in eachindex(edges)
        label = oe_label(edges, i)
        is_line_location_unknown(label, 0) && _label_disconnected_edge!(g, input, i, 0)
        is_line_location_unknown(label, 1) && _label_disconnected_edge!(g, input, i, 1)
    end
    return nothing
end

# Locates a disconnected edge against the ORIGINAL input area (design §3
# amendment 7: never against the reduced/collapsed linework), using both
# endpoints for robustness (port of `labelDisconnectedEdge` +
# `locateEdgeBothEnds`).
function _label_disconnected_edge!(g::OverlayGraph, input, i::Integer, gi::Integer)
    label = oe_label(g.edges, i)
    if !_input_is_area(input, gi)
        #-- non-area target: a disconnected edge must be EXTERIOR
        set_location_all!(label, gi, LOC_EXTERIOR)
        return nothing
    end
    loc_orig = _input_locate_in_area(input, gi, node_point(g.arr, he_origin(g.edges, i)))
    loc_dest = _input_locate_in_area(input, gi, node_point(g.arr, he_dest(g.edges, i)))
    is_int = loc_orig != LOC_EXTERIOR && loc_dest != LOC_EXTERIOR
    set_location_all!(label, gi, is_int ? LOC_INTERIOR : LOC_EXTERIOR)
    return nothing
end

# ## Result-area marking (ports of `markResultAreaEdges` / `unmarkDuplicate…`)

function _mark_result_area_edges!(g::OverlayGraph, op::_OverlayOpCode)
    edges = g.edges
    for i in eachindex(edges)
        _mark_in_result_area!(edges, i, op)
    end
    return nothing
end

# Port of `markInResultArea`: mark an edge whose right-side (boundary) or line
# location makes it part of the result-area boundary under `op`.
@inline function _mark_in_result_area!(edges, i::Integer, op::_OverlayOpCode)
    label = oe_label(edges, i)
    if is_boundary_either(label) && _is_result_of_op(op,
            oe_get_location_boundary_or_line(edges, i, 0, POS_RIGHT),
            oe_get_location_boundary_or_line(edges, i, 1, POS_RIGHT))
        oe_mark_in_result_area!(edges, i)
    end
    return nothing
end

# Port of `unmarkDuplicateEdgesFromResultArea`: an edge whose sym is also in the
# result area cancels (merges edge-adjacent result areas per polygon validity).
function _unmark_duplicate_edges_from_result_area!(g::OverlayGraph)
    edges = g.edges
    for i in eachindex(edges)
        oe_in_result_area_both(edges, i) && oe_unmark_from_result_area_both!(edges, i)
    end
    return nothing
end
