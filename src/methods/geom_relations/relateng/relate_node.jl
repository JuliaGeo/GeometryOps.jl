# # RelateNG node-edge topology
#
# Ports of JTS `RelateEdge.java` and `RelateNode.java`, in this order (JTS
# file boundaries preserved as clearly marked sections):
#
# 1. `RelateEdge` — one edge of the node "wheel": a direction around the
#    node carrying, per input geometry, the dimension and the left/on/right
#    locations of that geometry relative to the edge.
# 2. `RelateNode` — the wheel itself: the edges around a node in CCW order,
#    with insertion-or-merge and area-label propagation.
#
# Method order within each section parallels the Java file, so this file
# diffs against its Java counterparts.
#
# The Java classes are mutually recursive (`RelateEdge` holds its parent
# `RelateNode` to reach the node coordinate in `compareToEdge`); here the
# edge stores the symbolic `NodeKey` directly (design D2), and the manifold
# and `exact` flag the angle comparison needs are stored on the `RelateNode`
# (consistent with `RelateGeometry`/`AdjacentEdgeLocator`) and passed into
# `compare_to_edge`.

# Port of JTS Position constants (org.locationtech.jts.geom.Position):
# the indices for the location of a point relative to a directed edge.
const POS_ON = Int8(0)
const POS_LEFT = Int8(1)
const POS_RIGHT = Int8(2)

#==========================================================================
## RelateEdge (port of JTS RelateEdge.java)
==========================================================================#

# Port of RelateEdge.IS_FORWARD / IS_REVERSE.
const IS_FORWARD = true
const IS_REVERSE = false

#=
The dimension of an input geometry which is not known (Java
RelateEdge.DIM_UNKNOWN = -1). Numerically equal to `DIM_FALSE`, but kept as
its own name because the module-level `DIM_UNKNOWN` (from the IMPredicate
port in topology_predicate.jl) is JTS's *other* DIM_UNKNOWN constant, which
equals Dimension.DONTCARE = -3.
=#
const DIM_UNKNOWN_EDGE = Int8(-1)

# Indicates that the location is currently unknown (Java LOC_UNKNOWN =
# Location.NONE); `LOC_NONE` (de9im.jl) is used directly below.

"""
    RelateEdge{P}

An edge of a [`RelateNode`](@ref)'s wheel: the direction `node → dir_pt`,
labeled per input geometry with the dimension of the geometry element the
edge came from and the geometry's location on the left of, right of, and on
the edge. Unknown dimensions are `DIM_UNKNOWN_EDGE`; unknown locations are
`LOC_NONE`.

Port of JTS `RelateEdge`; the Java class stores its parent `RelateNode` to
reach the node coordinate, here the symbolic [`NodeKey`](@ref) is stored
directly (design D2).
"""
mutable struct RelateEdge{P}
    const node::NodeKey{P}
    const dir_pt::P

    a_dim::Int8
    a_loc_left::Int8
    a_loc_right::Int8
    a_loc_line::Int8

    b_dim::Int8
    b_loc_left::Int8
    b_loc_right::Int8
    b_loc_line::Int8
end

#=
Port of the static RelateEdge.create(node, dirPt, isA, dim, isForward)
factory: an area edge for `DIM_A`, a line edge otherwise.
=#
function relate_edge(node::NodeKey, dir_pt, is_a::Bool, dim::Integer, is_forward::Bool)
    if dim == DIM_A
        #-- create an area edge
        return RelateEdge(node, dir_pt, is_a, is_forward)
    end
    #-- create line edge
    return RelateEdge(node, dir_pt, is_a)
end

# Port of the static RelateEdge.findKnownEdgeIndex(edges, isA): index of the
# first edge with a known dimension for the geometry (1-based; 0 if none —
# Java returns -1).
function find_known_edge_index(edges::AbstractVector{<:RelateEdge}, is_a::Bool)
    for (i, e) in pairs(edges)
        is_known(e, is_a) && return i
    end
    return 0
end

# Port of the static RelateEdge.setAreaInterior(edges, isA) (renamed to
# distinguish it from the single-edge `set_area_interior!`).
function set_all_area_interior!(edges::AbstractVector{<:RelateEdge}, is_a::Bool)
    for e in edges
        set_area_interior!(e, is_a)
    end
    return nothing
end

# All-unknown base edge (the Java field initializers).
_relate_edge_unknown(node::NodeKey{P}, pt) where {P} = RelateEdge{P}(node, pt,
    DIM_UNKNOWN_EDGE, LOC_NONE, LOC_NONE, LOC_NONE,
    DIM_UNKNOWN_EDGE, LOC_NONE, LOC_NONE, LOC_NONE)

# Port of RelateEdge(node, pt, isA, isForward): an area edge.
function RelateEdge(node::NodeKey, pt, is_a::Bool, is_forward::Bool)
    e = _relate_edge_unknown(node, pt)
    set_locations_area!(e, is_a, is_forward)
    return e
end

# Port of RelateEdge(node, pt, isA): a line edge.
function RelateEdge(node::NodeKey, pt, is_a::Bool)
    e = _relate_edge_unknown(node, pt)
    set_locations_line!(e, is_a)
    return e
end

# Port of RelateEdge(node, pt, isA, locLeft, locRight, locLine): an area
# edge with explicit locations.
function RelateEdge(node::NodeKey, pt, is_a::Bool,
        loc_left::Integer, loc_right::Integer, loc_line::Integer)
    e = _relate_edge_unknown(node, pt)
    set_locations!(e, is_a, loc_left, loc_right, loc_line)
    return e
end

# Port of RelateEdge.setLocations (private; forces dim 2, as in Java).
function set_locations!(e::RelateEdge, is_a::Bool, loc_left::Integer,
        loc_right::Integer, loc_line::Integer)
    if is_a
        e.a_dim = 2
        e.a_loc_left = loc_left
        e.a_loc_right = loc_right
        e.a_loc_line = loc_line
    else
        e.b_dim = 2
        e.b_loc_left = loc_left
        e.b_loc_right = loc_right
        e.b_loc_line = loc_line
    end
    return nothing
end

# Port of RelateEdge.setLocationsLine (private).
function set_locations_line!(e::RelateEdge, is_a::Bool)
    if is_a
        e.a_dim = 1
        e.a_loc_left = LOC_EXTERIOR
        e.a_loc_right = LOC_EXTERIOR
        e.a_loc_line = LOC_INTERIOR
    else
        e.b_dim = 1
        e.b_loc_left = LOC_EXTERIOR
        e.b_loc_right = LOC_EXTERIOR
        e.b_loc_line = LOC_INTERIOR
    end
    return nothing
end

# Port of RelateEdge.setLocationsArea (private): a forward edge (the exiting
# edge of a CW-oriented corner) has the area interior on the right; a
# reverse (entering) edge has it on the left.
function set_locations_area!(e::RelateEdge, is_a::Bool, is_forward::Bool)
    loc_left = is_forward ? LOC_EXTERIOR : LOC_INTERIOR
    loc_right = is_forward ? LOC_INTERIOR : LOC_EXTERIOR
    if is_a
        e.a_dim = 2
        e.a_loc_left = loc_left
        e.a_loc_right = loc_right
        e.a_loc_line = LOC_BOUNDARY
    else
        e.b_dim = 2
        e.b_loc_left = loc_left
        e.b_loc_right = loc_right
        e.b_loc_line = LOC_BOUNDARY
    end
    return nothing
end

#=
Port of RelateEdge.compareToEdge(edgeDirPt): CCW angle comparison of this
edge's direction against `edge_dir_pt` around the node. The Java calls
`PolygonNodeTopology.compareAngle(node.getCoordinate(), dirPt, edgeDirPt)`;
here the apex is the symbolic node key and the comparison goes through the
kernel (`rk_compare_edge_dir`), which is why the manifold and `exact` flag
are threaded in (stored on the `RelateNode` by the callers).
=#
compare_to_edge(m::Manifold, e::RelateEdge, edge_dir_pt; exact) =
    rk_compare_edge_dir(m, e.node, e.dir_pt, edge_dir_pt; exact)

#=
Port of RelateEdge.merge(isA, dirPt, dim, isForward): merge the labeling of
a coincident (collinear, same-direction) edge of geometry `is_a` into this
edge. If the geometry is so far unknown on this edge, its dimension and
locations are simply installed; otherwise the dimension/on-location merge
(`merge_dim_edge_loc!`: area overrides line) and the side merges
(`merge_side_location!`: INTERIOR takes precedence) apply.
`dir_pt` is unused, as in Java (kept for signature parity).
=#
function merge_edge!(e::RelateEdge, is_a::Bool, dir_pt, dim::Integer, is_forward::Bool)
    loc_edge = LOC_INTERIOR
    loc_left = LOC_EXTERIOR
    loc_right = LOC_EXTERIOR
    if dim == DIM_A
        loc_edge = LOC_BOUNDARY
        loc_left = is_forward ? LOC_EXTERIOR : LOC_INTERIOR
        loc_right = is_forward ? LOC_INTERIOR : LOC_EXTERIOR
    end

    if !is_known(e, is_a)
        set_dimension!(e, is_a, dim)
        set_on!(e, is_a, loc_edge)
        set_left!(e, is_a, loc_left)
        set_right!(e, is_a, loc_right)
        return nothing
    end

    # Assert: node-dirpt is collinear with node-pt
    merge_dim_edge_loc!(e, is_a, loc_edge)
    merge_side_location!(e, is_a, POS_LEFT, loc_left)
    merge_side_location!(e, is_a, POS_RIGHT, loc_right)
    return nothing
end

#=
Port of RelateEdge.mergeDimEdgeLoc (private). Area edges override Line
edges. Merging edges of same dimension is a no-op for the dimension and on
location. But merging an area edge into a line edge sets the dimension to A
and the location to BOUNDARY.
=#
function merge_dim_edge_loc!(e::RelateEdge, is_a::Bool, loc_edge::Integer)
    #TODO: this logic needs work - ie handling A edges marked as Interior
    dim = loc_edge == LOC_BOUNDARY ? DIM_A : DIM_L
    if dim == DIM_A && edge_dim(e, is_a) == DIM_L
        set_dimension!(e, is_a, dim)
        set_on!(e, is_a, LOC_BOUNDARY)
    end
    return nothing
end

# Port of RelateEdge.mergeSideLocation (private): INTERIOR takes precedence
# over EXTERIOR.
function merge_side_location!(e::RelateEdge, is_a::Bool, pos::Integer, loc::Integer)
    curr_loc = edge_location(e, is_a, pos)
    if curr_loc != LOC_INTERIOR
        set_location!(e, is_a, pos, loc)
    end
    return nothing
end

# Port of RelateEdge.setDimension (private).
function set_dimension!(e::RelateEdge, is_a::Bool, dim::Integer)
    if is_a
        e.a_dim = dim
    else
        e.b_dim = dim
    end
    return nothing
end

# Port of RelateEdge.setLocation.
function set_location!(e::RelateEdge, is_a::Bool, pos::Integer, loc::Integer)
    if pos == POS_LEFT
        set_left!(e, is_a, loc)
    elseif pos == POS_RIGHT
        set_right!(e, is_a, loc)
    elseif pos == POS_ON
        set_on!(e, is_a, loc)
    end
    return nothing
end

# Port of RelateEdge.setAllLocations.
function set_all_locations!(e::RelateEdge, is_a::Bool, loc::Integer)
    set_left!(e, is_a, loc)
    set_right!(e, is_a, loc)
    set_on!(e, is_a, loc)
    return nothing
end

# Port of RelateEdge.setUnknownLocations: fill only the still-unknown
# positions with `loc`.
function set_unknown_locations!(e::RelateEdge, is_a::Bool, loc::Integer)
    if !is_known(e, is_a, POS_LEFT)
        set_location!(e, is_a, POS_LEFT, loc)
    end
    if !is_known(e, is_a, POS_RIGHT)
        set_location!(e, is_a, POS_RIGHT, loc)
    end
    if !is_known(e, is_a, POS_ON)
        set_location!(e, is_a, POS_ON, loc)
    end
    return nothing
end

# Port of RelateEdge.setLeft (private).
function set_left!(e::RelateEdge, is_a::Bool, loc::Integer)
    if is_a
        e.a_loc_left = loc
    else
        e.b_loc_left = loc
    end
    return nothing
end

# Port of RelateEdge.setRight (private).
function set_right!(e::RelateEdge, is_a::Bool, loc::Integer)
    if is_a
        e.a_loc_right = loc
    else
        e.b_loc_right = loc
    end
    return nothing
end

# Port of RelateEdge.setOn (private).
function set_on!(e::RelateEdge, is_a::Bool, loc::Integer)
    if is_a
        e.a_loc_line = loc
    else
        e.b_loc_line = loc
    end
    return nothing
end

# Port of RelateEdge.location(isA, position). (Java asserts unreachable for
# a bad position; here an ArgumentError.)
function edge_location(e::RelateEdge, is_a::Bool, position::Integer)
    if is_a
        position == POS_LEFT && return e.a_loc_left
        position == POS_RIGHT && return e.a_loc_right
        position == POS_ON && return e.a_loc_line
    else
        position == POS_LEFT && return e.b_loc_left
        position == POS_RIGHT && return e.b_loc_right
        position == POS_ON && return e.b_loc_line
    end
    throw(ArgumentError("invalid position: $position"))
end

# Port of RelateEdge.dimension (private).
edge_dim(e::RelateEdge, is_a::Bool) = is_a ? e.a_dim : e.b_dim

# Port of RelateEdge.isKnown(isA) (private): whether the geometry's
# dimension on this edge is known.
is_known(e::RelateEdge, is_a::Bool) =
    is_a ? e.a_dim != DIM_UNKNOWN_EDGE : e.b_dim != DIM_UNKNOWN_EDGE

# Port of RelateEdge.isKnown(isA, pos) (private): whether the location at
# `pos` is known.
is_known(e::RelateEdge, is_a::Bool, pos::Integer) = edge_location(e, is_a, pos) != LOC_NONE

# Port of RelateEdge.isInterior.
is_interior(e::RelateEdge, is_a::Bool, position::Integer) =
    edge_location(e, is_a, position) == LOC_INTERIOR

# Port of RelateEdge.setDimLocations.
function set_dim_locations!(e::RelateEdge, is_a::Bool, dim::Integer, loc::Integer)
    if is_a
        e.a_dim = dim
        e.a_loc_left = loc
        e.a_loc_right = loc
        e.a_loc_line = loc
    else
        e.b_dim = dim
        e.b_loc_left = loc
        e.b_loc_right = loc
        e.b_loc_line = loc
    end
    return nothing
end

# Port of RelateEdge.setAreaInterior(isA): all locations become INTERIOR
# (the dimension is untouched — a line edge inside an area keeps dim L).
function set_area_interior!(e::RelateEdge, is_a::Bool)
    if is_a
        e.a_loc_left = LOC_INTERIOR
        e.a_loc_right = LOC_INTERIOR
        e.a_loc_line = LOC_INTERIOR
    else
        e.b_loc_left = LOC_INTERIOR
        e.b_loc_right = LOC_INTERIOR
        e.b_loc_line = LOC_INTERIOR
    end
    return nothing
end

# Port of RelateEdge.toString (+ labelString/locationString), as a debugging
# aid. The node is symbolic, so crossing nodes print their NodeKey segment
# pair instead of a coordinate (`_point_string`, node_sections.jl).
function Base.show(io::IO, e::RelateEdge)
    print(io, _point_string(e.node), " - ", _point_string(e.dir_pt), " - ", _label_string(e))
end

_label_string(e::RelateEdge) =
    string("A:", _location_string(e, true), "/B:", _location_string(e, false))

_location_string(e::RelateEdge, is_a::Bool) = string(
    _loc_symbol(edge_location(e, is_a, POS_LEFT)),
    _loc_symbol(edge_location(e, is_a, POS_ON)),
    _loc_symbol(edge_location(e, is_a, POS_RIGHT)))

# Port of Location.toLocationSymbol.
function _loc_symbol(loc::Integer)
    loc == LOC_EXTERIOR && return 'e'
    loc == LOC_BOUNDARY && return 'b'
    loc == LOC_INTERIOR && return 'i'
    loc == LOC_NONE && return '-'
    throw(ArgumentError("Unknown location value: $loc"))
end

#==========================================================================
## RelateNode (port of JTS RelateNode.java)
==========================================================================#

"""
    RelateNode(m::Manifold, node::NodeKey; exact)

The topology at a node between the edges of two input geometries: a list of
the [`RelateEdge`](@ref)s around the node in CCW order, ordered by their CCW
angle with the positive X-axis.

Port of JTS `RelateNode`; the Java class is keyed by the node `Coordinate`,
here by the symbolic [`NodeKey`](@ref) (design D2). The manifold and the
`exact` flag (absent in Java) are stored for the edge-angle comparisons in
[`add_edge!`](@ref add_edge!(::RelateNode, ::Bool, ::Any, ::Integer, ::Bool)).
"""
struct RelateNode{M <: Manifold, E, P}
    m::M
    exact::E
    node::NodeKey{P}
    edges::Vector{RelateEdge{P}}
end

RelateNode(m::Manifold, node::NodeKey{P}; exact) where {P} =
    RelateNode(m, exact, node, RelateEdge{P}[])

# Port of RelateNode.getCoordinate. The Java method returns the node
# Coordinate; here the node is its symbolic NodeKey (design D2).
get_coordinate(n::RelateNode) = n.node

# Port of RelateNode.getEdges.
get_edges(n::RelateNode) = n.edges

# Port of RelateNode.addEdges(List<NodeSection>).
function add_edges!(n::RelateNode, nss::AbstractVector{<:NodeSection})
    for ns in nss
        add_edges!(n, ns)
    end
    return nothing
end

# Port of RelateNode.addEdges(NodeSection).
function add_edges!(n::RelateNode, ns::NodeSection)
    dim = section_dim(ns)
    isa_g = is_a(ns)
    if dim == DIM_L
        add_line_edge!(n, isa_g, get_vertex(ns, 0))
        add_line_edge!(n, isa_g, get_vertex(ns, 1))
    elseif dim == DIM_A
        #-- assumes node edges have CW orientation (as per JTS norm)
        #-- entering edge - interior on L
        e0 = add_area_edge!(n, isa_g, get_vertex(ns, 0), false)
        #-- exiting edge - interior on R
        e1 = add_area_edge!(n, isa_g, get_vertex(ns, 1), true)
        # Zero-length edges are skipped (Java addEdge returns null; they
        # only arise from invalid rings with repeated points).
        (e0 === nothing || e1 === nothing) && return nothing

        index0 = findfirst(e -> e === e0, n.edges)
        index1 = findfirst(e -> e === e1, n.edges)
        update_edges_in_area!(n, isa_g, index0, index1)
        update_if_area_prev!(n, isa_g, index0)
        update_if_area_next!(n, isa_g, index1)
    end
    return nothing
end

# Port of RelateNode.updateEdgesInArea (private): mark every edge strictly
# between the entering and exiting edge (in CCW order) as area interior.
function update_edges_in_area!(n::RelateNode, is_a::Bool, index_from::Integer, index_to::Integer)
    index = next_index(n.edges, index_from)
    while index != index_to
        edge = n.edges[index]
        set_area_interior!(edge, is_a)
        index = next_index(n.edges, index)
    end
    return nothing
end

# Port of RelateNode.updateIfAreaPrev (private): if the CW-previous edge has
# the area interior on its left, this edge lies inside the area.
function update_if_area_prev!(n::RelateNode, is_a::Bool, index::Integer)
    index_prev = prev_index(n.edges, index)
    edge_prev = n.edges[index_prev]
    if is_interior(edge_prev, is_a, POS_LEFT)
        edge = n.edges[index]
        set_area_interior!(edge, is_a)
    end
    return nothing
end

# Port of RelateNode.updateIfAreaNext (private): if the CCW-next edge has
# the area interior on its right, this edge lies inside the area.
function update_if_area_next!(n::RelateNode, is_a::Bool, index::Integer)
    index_next = next_index(n.edges, index)
    edge_next = n.edges[index_next]
    if is_interior(edge_next, is_a, POS_RIGHT)
        edge = n.edges[index]
        set_area_interior!(edge, is_a)
    end
    return nothing
end

# Port of RelateNode.addLineEdge (private).
add_line_edge!(n::RelateNode, is_a::Bool, dir_pt) =
    add_edge!(n, is_a, dir_pt, DIM_L, false)

# Port of RelateNode.addAreaEdge (private).
add_area_edge!(n::RelateNode, is_a::Bool, dir_pt, is_forward::Bool) =
    add_edge!(n, is_a, dir_pt, DIM_A, is_forward)

"""
    add_edge!(n::RelateNode, is_a, dir_pt, dim, is_forward)

Adds or merges an edge to the node, keeping the wheel sorted by CCW angle
with the positive X-axis. `dim` is the dimension of the geometry element
containing the edge, `is_forward` the direction of the edge. Returns the
created or merged edge for this point, or `nothing` for a malformed
(`nothing` or zero-length) input edge.

Port of RelateNode.addEdge.
"""
function add_edge!(n::RelateNode, is_a::Bool, dir_pt, dim::Integer, is_forward::Bool)
    #-- check for well-formed edge - skip null or zero-len input
    dir_pt === nothing && return nothing
    # Java: nodePt.equals2D(dirPt). A proper-crossing node lies strictly in
    # the interior of its defining segments, so a (vertex) direction point
    # can never coincide with it; only vertex nodes need the check.
    (!n.node.is_crossing && _equals2(n.node.pt, dir_pt)) && return nothing

    insert_index = 0
    for (i, e) in pairs(n.edges)
        comp = compare_to_edge(n.m, e, dir_pt; exact = n.exact)
        if comp == 0
            merge_edge!(e, is_a, dir_pt, dim, is_forward)
            return e
        end
        if comp > 0
            #-- found further edge, so insert a new edge at this position
            insert_index = i
            break
        end
    end
    #-- add a new edge
    e = relate_edge(n.node, dir_pt, is_a, dim, is_forward)
    if insert_index == 0
        #-- add edge at end of list
        push!(n.edges, e)
    else
        #-- add edge before higher edge found
        insert!(n.edges, insert_index, e)
    end
    return e
end

"""
    finish!(n::RelateNode, is_area_interior_a::Bool, is_area_interior_b::Bool)

Computes the final topology for the edges around this node. Although nodes
lie on the boundary of areas or the interior of lines, in a mixed GC they
may also lie in the interior of an area. This changes the locations of the
sides and line to Interior.

Port of RelateNode.finish.
"""
function finish!(n::RelateNode, is_area_interior_a::Bool, is_area_interior_b::Bool)
    finish_node!(n, true, is_area_interior_a)
    finish_node!(n, false, is_area_interior_b)
    return nothing
end

# Port of RelateNode.finishNode (private).
function finish_node!(n::RelateNode, is_a::Bool, is_area_interior::Bool)
    if is_area_interior
        set_all_area_interior!(n.edges, is_a)
    else
        start_index = find_known_edge_index(n.edges, is_a)
        #-- only interacting nodes are finished, so this should never happen
        #Assert: start_index > 0, "Node does not have AB interaction"
        propagate_side_locations!(n, is_a, start_index)
    end
    return nothing
end

# Port of RelateNode.propagateSideLocations (private): walk the wheel CCW
# from the first known edge, filling unknown locations with the latest known
# LEFT location (the location of the angular sector CCW of each edge).
function propagate_side_locations!(n::RelateNode, is_a::Bool, start_index::Integer)
    curr_loc = edge_location(n.edges[start_index], is_a, POS_LEFT)
    #-- edges are stored in CCW order
    index = next_index(n.edges, start_index)
    while index != start_index
        e = n.edges[index]
        set_unknown_locations!(e, is_a, curr_loc)
        curr_loc = edge_location(e, is_a, POS_LEFT)
        index = next_index(n.edges, index)
    end
    return nothing
end

# Ports of the static RelateNode.prevIndex / nextIndex (1-based, circular).
prev_index(edges::AbstractVector, index::Integer) =
    index > 1 ? index - 1 : length(edges)

next_index(edges::AbstractVector, i::Integer) =
    i >= length(edges) ? 1 : i + 1

# Port of RelateNode.toString, as a debugging aid.
function Base.show(io::IO, n::RelateNode)
    print(io, "Node[", _point_string(n.node), "]:")
    for e in n.edges
        print(io, "\n", e)
    end
end

# Port of RelateNode.hasExteriorEdge(isA): whether any edge has the geometry
# in its exterior on either side.
function has_exterior_edge(n::RelateNode, is_a::Bool)
    for e in n.edges
        if LOC_EXTERIOR == edge_location(e, is_a, POS_LEFT) ||
                LOC_EXTERIOR == edge_location(e, is_a, POS_RIGHT)
            return true
        end
    end
    return false
end
