# # RelateNG node sections
#
# Port of JTS `NodeSection.java` (data shape only, for now). Created in
# Task 11 because `AdjacentEdgeLocator` builds `NodeSection`s to test polygon
# edge adjacency; the full `NodeSection` API (comparator, `EdgeAngleComparator`)
# and the `NodeSections` collector land in Task 15.

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

# Port of NodeSection.getVertex(i): the incident edge vertex before (0) or
# after (1) the node, or `nothing` if that edge does not exist.
get_vertex(ns::NodeSection, i::Integer) = i == 0 ? ns.v0 : ns.v1
