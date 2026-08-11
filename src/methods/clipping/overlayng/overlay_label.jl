# # OverlayLabel — the per-edge topological label (port of JTS `OverlayLabel`)
#
# Phase 2a of the OverlayNG port (design doc `2026-07-16-overlayng-noding-substrate.md`,
# §3). A faithful port of `operation/overlayng/OverlayLabel.java`: a mutable plain
# struct (NOT bit-packed), one instance shared between the two oppositely-oriented
# `OverlayEdge`s of a symmetric pair. Orientation-sensitive accessors take the
# containing edge's `is_forward` flag and swap Left/Right when it is false.
#
# A label records, for each of the (up to) two input geometries (index `0` = A,
# `1` = B, matching JTS), one of four states via its `dim` field:
#   - Boundary  (`DIM_A`)       : an area-boundary edge; `loc_left`/`loc_right` set,
#                                 `loc_line = LOC_INTERIOR`, `is_hole` = ring role.
#   - Collapse  (`DIM_COLLAPSE`): two or more coincident area edges summed to a
#                                 zero depth delta; `loc_line` filled later.
#   - Line      (`DIM_L`)       : a linear edge; only `loc_line` is meaningful.
#   - NotPart   (`DIM_NOT_PART`): the edge is not part of this input.
#
# Everything here is internal to GeometryOps — nothing is exported.

# ## Edge-state dimension codes (port of `OverlayLabel.DIM_*`)
#
# `DIM_BOUNDARY` / `DIM_LINE` reuse the DE9IM `DIM_A` (2) / `DIM_L` (1) values,
# which JTS's `OverlayLabel` constants numerically coincide with (JTS
# `DIM_BOUNDARY == Dimension.A`, `DIM_LINE == Dimension.L`); only the two overlay-
# specific states need fresh values.
const DIM_NOT_PART = Int8(-1)    # JTS OverlayLabel.DIM_NOT_PART == DIM_UNKNOWN == Dimension.FALSE
const DIM_COLLAPSE = Int8(3)     # JTS OverlayLabel.DIM_COLLAPSE

# `Position` codes (`POS_ON`/`POS_LEFT`/`POS_RIGHT`) and location codes
# (`LOC_INTERIOR`/`LOC_BOUNDARY`/`LOC_EXTERIOR`/`LOC_NONE`) are the RelateNG ports
# already in scope (relate_node.jl / de9im.jl); `LOC_NONE` is JTS `LOC_UNKNOWN`.

"""
    OverlayLabel

The topological label of one edge of the overlay graph (port of JTS
`OverlayLabel`). Mutable plain fields, one instance shared between a symmetric
`OverlayEdge` pair; the `is_forward`-parameterized accessors swap Left/Right for
the reverse half-edge. Index `0` selects input A, index `1` input B.
"""
mutable struct OverlayLabel
    a_dim      :: Int8
    a_is_hole  :: Bool
    a_loc_left :: Int8
    a_loc_right:: Int8
    a_loc_line :: Int8

    b_dim      :: Int8
    b_is_hole  :: Bool
    b_loc_left :: Int8
    b_loc_right:: Int8
    b_loc_line :: Int8
end

# Uninitialized label: both inputs NotPart, all locations unknown.
OverlayLabel() = OverlayLabel(DIM_NOT_PART, false, LOC_NONE, LOC_NONE, LOC_NONE,
                              DIM_NOT_PART, false, LOC_NONE, LOC_NONE, LOC_NONE)

Base.copy(l::OverlayLabel) = OverlayLabel(l.a_dim, l.a_is_hole, l.a_loc_left, l.a_loc_right, l.a_loc_line,
                                          l.b_dim, l.b_is_hole, l.b_loc_left, l.b_loc_right, l.b_loc_line)

# ## Initializers (port of `initBoundary`/`initCollapse`/`initLine`/`initNotPart`)

# Area-boundary edge: side locations known, `loc_line = LOC_INTERIOR`.
function init_boundary!(l::OverlayLabel, index::Integer, loc_left, loc_right, is_hole::Bool)
    if index == 0
        l.a_dim = DIM_A; l.a_is_hole = is_hole
        l.a_loc_left = loc_left; l.a_loc_right = loc_right; l.a_loc_line = LOC_INTERIOR
    else
        l.b_dim = DIM_A; l.b_is_hole = is_hole
        l.b_loc_left = loc_left; l.b_loc_right = loc_right; l.b_loc_line = LOC_INTERIOR
    end
    return l
end

# Collapsed area edge: location unknown, resolved later from the graph topology.
function init_collapse!(l::OverlayLabel, index::Integer, is_hole::Bool)
    if index == 0
        l.a_dim = DIM_COLLAPSE; l.a_is_hole = is_hole
    else
        l.b_dim = DIM_COLLAPSE; l.b_is_hole = is_hole
    end
    return l
end

# Line edge: only the line location is meaningful, initialized unknown.
function init_line!(l::OverlayLabel, index::Integer)
    if index == 0
        l.a_dim = DIM_L; l.a_loc_line = LOC_NONE
    else
        l.b_dim = DIM_L; l.b_loc_line = LOC_NONE
    end
    return l
end

# Not part of this input (locations assumed already unknown).
function init_not_part!(l::OverlayLabel, index::Integer)
    index == 0 ? (l.a_dim = DIM_NOT_PART) : (l.b_dim = DIM_NOT_PART)
    return l
end

# ## Location setters (used during label propagation, phase 2b)

function set_location_line!(l::OverlayLabel, index::Integer, loc)
    index == 0 ? (l.a_loc_line = loc) : (l.b_loc_line = loc)
    return l
end

function set_location_all!(l::OverlayLabel, index::Integer, loc)
    if index == 0
        l.a_loc_line = loc; l.a_loc_left = loc; l.a_loc_right = loc
    else
        l.b_loc_line = loc; l.b_loc_left = loc; l.b_loc_right = loc
    end
    return l
end

# A collapsed edge with no boundary information takes its parent ring role:
# a hole collapse is INTERIOR, a shell collapse is EXTERIOR.
function set_location_collapse!(l::OverlayLabel, index::Integer)
    loc = is_hole(l, index) ? LOC_INTERIOR : LOC_EXTERIOR
    index == 0 ? (l.a_loc_line = loc) : (l.b_loc_line = loc)
    return l
end

# ## State predicates (port of the `OverlayLabel` boolean surface)

dimension(l::OverlayLabel, index::Integer) = index == 0 ? l.a_dim : l.b_dim

is_line(l::OverlayLabel) = l.a_dim == DIM_L || l.b_dim == DIM_L
is_line(l::OverlayLabel, index::Integer) = (index == 0 ? l.a_dim : l.b_dim) == DIM_L
is_linear(l::OverlayLabel, index::Integer) =
    (d = index == 0 ? l.a_dim : l.b_dim; d == DIM_L || d == DIM_COLLAPSE)
is_known(l::OverlayLabel, index::Integer) = (index == 0 ? l.a_dim : l.b_dim) != DIM_NOT_PART
is_not_part(l::OverlayLabel, index::Integer) = (index == 0 ? l.a_dim : l.b_dim) == DIM_NOT_PART

is_boundary_either(l::OverlayLabel) = l.a_dim == DIM_A || l.b_dim == DIM_A
is_boundary_both(l::OverlayLabel)   = l.a_dim == DIM_A && l.b_dim == DIM_A
is_boundary(l::OverlayLabel, index::Integer) = (index == 0 ? l.a_dim : l.b_dim) == DIM_A

# A collapse coincident with the other input's (non-collapsed) boundary.
is_boundary_collapse(l::OverlayLabel) = is_line(l) ? false : !is_boundary_both(l)

# Two areas touching along their common boundary (opposite Right locations).
is_boundary_touch(l::OverlayLabel) =
    is_boundary_both(l) &&
    get_location(l, 0, POS_RIGHT, true) != get_location(l, 1, POS_RIGHT, true)

is_boundary_singleton(l::OverlayLabel) =
    (l.a_dim == DIM_A && l.b_dim == DIM_NOT_PART) ||
    (l.b_dim == DIM_A && l.a_dim == DIM_NOT_PART)

is_line_location_unknown(l::OverlayLabel, index::Integer) =
    (index == 0 ? l.a_loc_line : l.b_loc_line) == LOC_NONE

is_line_in_area(l::OverlayLabel, index::Integer) =
    (index == 0 ? l.a_loc_line : l.b_loc_line) == LOC_INTERIOR

is_hole(l::OverlayLabel, index::Integer) = index == 0 ? l.a_is_hole : l.b_is_hole

is_collapse(l::OverlayLabel, index::Integer) = dimension(l, index) == DIM_COLLAPSE

is_interior_collapse(l::OverlayLabel) =
    (l.a_dim == DIM_COLLAPSE && l.a_loc_line == LOC_INTERIOR) ||
    (l.b_dim == DIM_COLLAPSE && l.b_loc_line == LOC_INTERIOR)

is_collapse_and_not_part_interior(l::OverlayLabel) =
    (l.a_dim == DIM_COLLAPSE && l.b_dim == DIM_NOT_PART && l.b_loc_line == LOC_INTERIOR) ||
    (l.b_dim == DIM_COLLAPSE && l.a_dim == DIM_NOT_PART && l.a_loc_line == LOC_INTERIOR)

is_line_interior(l::OverlayLabel, index::Integer) =
    (index == 0 ? l.a_loc_line : l.b_loc_line) == LOC_INTERIOR

get_line_location(l::OverlayLabel, index::Integer) = index == 0 ? l.a_loc_line : l.b_loc_line

has_sides(l::OverlayLabel, index::Integer) =
    index == 0 ? (l.a_loc_left != LOC_NONE || l.a_loc_right != LOC_NONE) :
                 (l.b_loc_left != LOC_NONE || l.b_loc_right != LOC_NONE)

# ## Orientation-sensitive location accessor (the shared-label L/R swap)

"""
    get_location(label, index, position, is_forward) -> Int8

The location of `position` (`POS_LEFT` / `POS_RIGHT` / `POS_ON`) of input `index`,
for a containing half-edge whose orientation is `is_forward`. When the edge is the
reverse half-edge (`is_forward == false`) the Left/Right stored sides are swapped,
so a single label serves both members of a symmetric pair (port of JTS
`OverlayLabel.getLocation(index, position, isForward)`).
"""
function get_location(l::OverlayLabel, index::Integer, position::Integer, is_forward::Bool)
    if index == 0
        position == POS_LEFT  && return is_forward ? l.a_loc_left  : l.a_loc_right
        position == POS_RIGHT && return is_forward ? l.a_loc_right : l.a_loc_left
        position == POS_ON    && return l.a_loc_line
    else
        position == POS_LEFT  && return is_forward ? l.b_loc_left  : l.b_loc_right
        position == POS_RIGHT && return is_forward ? l.b_loc_right : l.b_loc_left
        position == POS_ON    && return l.b_loc_line
    end
    return LOC_NONE
end

# The linear (ON) location of a source (port of the 2-arg `getLocation`).
get_location(l::OverlayLabel, index::Integer) = index == 0 ? l.a_loc_line : l.b_loc_line

# For a boundary edge the side location; otherwise the line location — the
# quantity `markInResultArea` tests (port of `getLocationBoundaryOrLine`).
get_location_boundary_or_line(l::OverlayLabel, index::Integer, position::Integer, is_forward::Bool) =
    is_boundary(l, index) ? get_location(l, index, position, is_forward) :
                            get_line_location(l, index)
