# # RelateNG node sections
#
# Ports of JTS `NodeSection.java` and `NodeSections.java`, in this order
# (JTS file boundaries preserved as clearly marked sections):
#
# 1. `NodeSection`   — a geometry component's contribution to a node
#    (struct created in Task 11 for `AdjacentEdgeLocator`; full API Task 15)
# 2. `NodeSections`  — the collector of all sections at one node
#
# Method order within each section parallels the Java file, so this file
# diffs against its Java counterparts.

#==========================================================================
## NodeSection (port of JTS NodeSection.java)
==========================================================================#

"""
    NodeSection{P, G}

Represents a computed node along with the incident edges on either side of
it (if they exist). This captures the information about a node in a geometry
component required to determine the component's contribution to the node
topology. A node in an area geometry always has edges on both sides of the
node. A node in a linear geometry may have one or other incident edge
missing, if the node occurs at an endpoint of the line.

The edges of an area node are assumed to be provided with CW-shell
orientation (as per JTS norm). This must be enforced by the caller.

Port of JTS `NodeSection`, with one symbolic twist (design D2): the Java
class stores the node as a `Coordinate`; here `node` is a [`NodeKey`](@ref),
so proper-crossing nodes never need a constructed coordinate. `v0`/`v1` are
coordinate tuples of type `P` (or `nothing` for a missing incident edge at a
line endpoint); `polygonal` is the parent polygonal geometry of an area
section, or `nothing` if the section is not on a polygon boundary.

The field order matches the Java constructor argument order
`(isA, dimension, id, ringId, poly, isNodeAtVertex, v0, nodePt, v1)`.
(The struct is declared before the comparator helpers below because their
signatures reference it; the Java file declares `EdgeAngleComparator` and
`isAreaArea` first.)
"""
struct NodeSection{P, G}
    is_a::Bool
    dim::Int8
    id::Int32
    ring_id::Int32
    polygonal::G
    is_node_at_vertex::Bool
    v0::Union{P, Nothing}
    node::NodeKey{P}
    v1::Union{P, Nothing}
end

"""
    edge_angle_compare(m::Manifold, ns1::NodeSection, ns2::NodeSection; exact)

Compares sections by the angle the entering edge (`get_vertex(ns, 0)`) makes
with the positive X axis at the node, angles increasing CCW.

Port of `NodeSection.EdgeAngleComparator` (a static nested `Comparator`
class in Java, `compareAngle(ns1.nodePt, ns1.getVertex(0), ns2.getVertex(0))`);
here a comparator function over `rk_compare_edge_dir` with the
symbolic node of `ns1` as apex, taking the manifold and `exact` flag the
kernel comparison needs. Use as a sort predicate via
`lt = (a, b) -> edge_angle_compare(m, a, b; exact) < 0`.

At a crossing node (a `NodeKey` with `is_crossing`) the sections' `v0` are
normally among the four endpoints of the node's defining segments (the
sections are built from those segments themselves), where the comparison
is derived from the original endpoints. Sections merged onto the node by
the D3 coincidence pass (`TopologyComputer`) may carry foreign directions;
`rk_compare_edge_dir` then compares around the exact rational apex.
"""
edge_angle_compare(m::Manifold, ns1::NodeSection, ns2::NodeSection; exact) =
    rk_compare_edge_dir(m, ns1.node, get_vertex(ns1, 0), get_vertex(ns2, 0); exact)

# Port of NodeSection.isAreaArea: whether both sections are area sections.
is_area_area(a::NodeSection, b::NodeSection) =
    section_dim(a) == DIM_A && section_dim(b) == DIM_A

# Port of NodeSection.getVertex(i): the incident edge vertex before (0) or
# after (1) the node, or `nothing` if that edge does not exist.
get_vertex(ns::NodeSection, i::Integer) = i == 0 ? ns.v0 : ns.v1

# Port of NodeSection.nodePt. The Java method returns the node Coordinate;
# here the node is its symbolic NodeKey (design D2).
node_pt(ns::NodeSection) = ns.node

# Port of NodeSection.dimension.
section_dim(ns::NodeSection) = ns.dim

# Port of NodeSection.id.
section_id(ns::NodeSection) = ns.id

# Port of NodeSection.ringId.
ring_id(ns::NodeSection) = ns.ring_id

"""
    get_polygonal(ns::NodeSection)

Gets the polygon this section is part of.
Will be `nothing` if section is not on a polygon boundary.

Port of NodeSection.getPolygonal.
"""
get_polygonal(ns::NodeSection) = ns.polygonal

# Port of NodeSection.isShell.
is_shell(ns::NodeSection) = ns.ring_id == 0

# Port of NodeSection.isArea.
is_area(ns::NodeSection) = ns.dim == DIM_A

# Port of NodeSection.isA.
is_a(ns::NodeSection) = ns.is_a

# Port of NodeSection.isSameGeometry: both sections from the same input
# geometry (A or B).
is_same_geometry(ns::NodeSection, other::NodeSection) = is_a(ns) == is_a(other)

# Port of NodeSection.isSamePolygon: both sections from the same polygon
# element of the same input geometry. (Element ids are only unique within
# one input geometry.)
is_same_polygon(ns::NodeSection, other::NodeSection) =
    is_a(ns) == is_a(other) && section_id(ns) == section_id(other)

# Port of NodeSection.isNodeAtVertex.
is_node_at_vertex(ns::NodeSection) = ns.is_node_at_vertex

# Port of NodeSection.isProper: a node is "proper" for a section if it lies
# in the interior of an edge (i.e. NOT at a vertex of the component).
is_proper(ns::NodeSection) = !ns.is_node_at_vertex

# Port of the static NodeSection.isProper(a, b): both sections proper.
is_proper(a::NodeSection, b::NodeSection) = is_proper(a) && is_proper(b)

# Port of NodeSection.toString (+ edgeRep), as a debugging aid. The node is
# symbolic, so crossing nodes print their NodeKey segment pair instead of a
# coordinate.
function Base.show(io::IO, ns::NodeSection)
    at_vertex_ind = ns.is_node_at_vertex ? "-V-" : "---"
    poly_id = ns.id >= 0 ? "[$(ns.id):$(ns.ring_id)]" : ""
    print(io, geom_name(ns.is_a), ns.dim, poly_id, ": ",
        _edge_rep(ns.v0, ns.node), " ", at_vertex_ind, " ", _edge_rep(ns.node, ns.v1))
end

_edge_rep(p0, p1) =
    (p0 === nothing || p1 === nothing) ? "null" : string(_pt_rep(p0), " - ", _pt_rep(p1))
_pt_rep(p) = "($(GI.x(p)) $(GI.y(p)))"
_pt_rep(k::NodeKey) = k.is_crossing ?
    string("X[", _pt_rep(k.pt), " - ", _pt_rep(k.a1), " × ", _pt_rep(k.b0), " - ", _pt_rep(k.b1), "]") :
    _pt_rep(k.pt)

"""
    compare_to(ns::NodeSection, other::NodeSection)

Compare node sections by parent geometry, dimension, element id and ring id,
and edge vertices. Sections are assumed to be at the same node point.
Returns a negative/zero/positive `Int` (Java `compareTo` contract).

Port of NodeSection.compareTo.
"""
function compare_to(ns::NodeSection, o::NodeSection)
    # Assert: ns.node == o.node

    #-- sort A before B
    if ns.is_a != o.is_a
        ns.is_a && return -1
        return 1
    end
    #-- sort on dimensions
    comp_dim = _compare_int(ns.dim, o.dim)
    comp_dim != 0 && return comp_dim

    #-- sort on id and ring id
    comp_id = _compare_int(ns.id, o.id)
    comp_id != 0 && return comp_id

    comp_ring_id = _compare_int(ns.ring_id, o.ring_id)
    comp_ring_id != 0 && return comp_ring_id

    #-- sort on edge coordinates
    comp_v0 = _compare_with_null(ns.v0, o.v0)
    comp_v0 != 0 && return comp_v0

    return _compare_with_null(ns.v1, o.v1)
end

# Java Integer.compare.
_compare_int(a, b) = a < b ? -1 : (a > b ? 1 : 0)

# Port of NodeSection.compareWithNull: `nothing` (Java null) sorts below
# any coordinate; coordinates compare via Coordinate.compareTo.
_compare_with_null(::Nothing, ::Nothing) = 0
_compare_with_null(::Nothing, v1) = -1   # null is lower than non-null
_compare_with_null(v0, ::Nothing) = 1
_compare_with_null(v0, v1) = _compare_pt(v0, v1)

# JTS Coordinate.compareTo: lexicographic on (x, y). (Lived in the Task-11
# AdjacentEdgeLocator slice in point_locator.jl until Task 15.)
function _compare_pt(p, q)
    GI.x(p) < GI.x(q) && return -1
    GI.x(p) > GI.x(q) && return 1
    GI.y(p) < GI.y(q) && return -1
    GI.y(p) > GI.y(q) && return 1
    return 0
end

#==========================================================================
## NodeSections (port of JTS NodeSections.java)
==========================================================================#

"""
    NodeSections(node::NodeKey)

Collects the [`NodeSection`](@ref)s of all geometry components incident on
one node, and assembles them into the node's edge topology
([`create_node`](@ref)).

Port of JTS `NodeSections`; the Java class is keyed by the node
`Coordinate`, here by the symbolic [`NodeKey`](@ref) (design D2).
"""
mutable struct NodeSections{P}
    const node::NodeKey{P}
    const sections::Vector{NodeSection}
end

NodeSections(node::NodeKey) = NodeSections(node, NodeSection[])

# Port of NodeSections.getCoordinate. The Java method returns the node
# Coordinate; here the node is its symbolic NodeKey (design D2).
get_coordinate(nss::NodeSections) = nss.node

# Port of NodeSections.addNodeSection.
function add_node_section!(nss::NodeSections, e::NodeSection)
    push!(nss.sections, e)
    return nothing
end

# Port of NodeSections.hasInteractionAB: whether both input geometries
# contribute a section at this node.
function has_interaction_ab(nss::NodeSections)
    found_a = false
    found_b = false
    for ns in nss.sections
        if is_a(ns)
            found_a = true
        else
            found_b = true
        end
        found_a && found_b && return true
    end
    return false
end

"""
    get_polygonal(nss::NodeSections, is_a::Bool)

The parent polygonal geometry of the first section of input geometry
`is_a` that has one, or `nothing`.

Port of NodeSections.getPolygonal(boolean isA).
"""
function get_polygonal(nss::NodeSections, is_a_target::Bool)
    for ns in nss.sections
        if is_a(ns) == is_a_target
            poly = get_polygonal(ns)
            poly !== nothing && return poly
        end
    end
    return nothing
end

"""
    create_node(m::Manifold, nss::NodeSections; exact)

Creates the node topology: prepares the sections, builds a
[`RelateNode`](@ref) at the node and feeds it the sections via
`add_edges!`. Per-polygon section groups are first rewritten into
maximal-ring structure by [`polygon_node_convert`](@ref) (the
`PolygonNodeConverter.convert` port). Returns the assembled node.

Port of NodeSections.createNode. The manifold/`exact` parameters (absent in
Java, where `createNode()` is nullary) are threaded through for the angle
comparisons in the converter and the node's edge wheel. (`RelateNode` is
defined in relate_node.jl, included after this file; the reference resolves
at call time.)
"""
function create_node(m::Manifold, nss::NodeSections; exact)
    prepare_sections!(nss)

    node = RelateNode(m, nss.node; exact)
    i = 1
    while i <= length(nss.sections)
        ns = nss.sections[i]
        #-- if there multiple polygon sections incident at node convert them to maximal-ring structure
        if is_area(ns) && _has_multiple_polygon_sections(nss.sections, i)
            poly_sections = _collect_polygon_sections(nss.sections, i)
            ns_convert = polygon_node_convert(m, poly_sections; exact)
            add_edges!(node, ns_convert)
            i += length(poly_sections)
        else
            #-- the most common case is a line or a single polygon ring section
            add_edges!(node, ns)
            i += 1
        end
    end
    return node
end

"""
    prepare_sections!(nss::NodeSections)

Sorts the sections (by [`compare_to`](@ref), the Java natural ordering)
so that:

- lines are before areas
- edges from the same polygon are contiguous

Port of NodeSections.prepareSections.
"""
function prepare_sections!(nss::NodeSections)
    sort!(nss.sections; lt = (a, b) -> compare_to(a, b) < 0)
    #TODO: remove duplicate sections
    return nothing
end

# Port of NodeSections.hasMultiplePolygonSections (1-based index).
function _has_multiple_polygon_sections(sections::Vector{<:NodeSection}, i::Integer)
    #-- if last section can only be one
    i >= length(sections) && return false
    #-- check if there are at least two sections for same polygon
    ns = sections[i]
    ns_next = sections[i + 1]
    return is_same_polygon(ns, ns_next)
end

# Port of NodeSections.collectPolygonSections (1-based index).
function _collect_polygon_sections(sections::Vector{<:NodeSection}, i::Integer)
    poly_sections = NodeSection[]
    #-- note ids are only unique to a geometry
    poly_section = sections[i]
    while i <= length(sections) && is_same_polygon(poly_section, sections[i])
        push!(poly_sections, sections[i])
        i += 1
    end
    return poly_sections
end
