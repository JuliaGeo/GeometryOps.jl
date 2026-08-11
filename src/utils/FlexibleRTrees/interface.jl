# # SpatialTreeInterface

using ..SpatialTreeInterface
import ..SpatialTreeInterface: isspatialtree, isleaf, nchild, getchild,
    child_indices_extents, depth_first_search

"""
    RTreeNode{T, E}

A cursor into one node of an [`RTree`](@ref): the tree, the node's level
(0-based; the children of a level-`l` node live in `levels[l + 1]`), its
position within that level, and its extent.  All SpatialTreeInterface
methods traverse the tree through these cursors.  The children of one node
occupy one contiguous run of the next level's extent vector, which each
per-child method resolves once per node and then indexes into.  At the leaf
level, `child_indices_extents` maps leaf slots through `tree.indices`, so
queries return indices into the original collection despite the packed
reordering.
"""
struct RTreeNode{T <: RTree, E <: Extents.Extent}
    tree::T
    level::Int      # 0-based; children of a level-l node live in levels[l + 1]
    index::Int      # position within its level
    extent::E
end

Extents.extent(node::RTreeNode) = node.extent

isspatialtree(::Type{<:RTree}) = true
isspatialtree(::Type{<:RTreeNode}) = true

@inline _child_extents(node::RTreeNode) = node.tree.levels[node.level + 1]

@inline function _child_range(node::RTreeNode, child_extents)
    start_idx = (node.index - 1) * node.tree.nodecapacity + 1
    stop_idx = min(start_idx + node.tree.nodecapacity - 1, length(child_extents))
    return start_idx:stop_idx
end

isleaf(node::RTreeNode) = node.level == length(node.tree.levels) - 1

nchild(node::RTreeNode) = length(_child_range(node, _child_extents(node)))

function getchild(node::RTreeNode, i::Int)
    child_index = (node.index - 1) * node.tree.nodecapacity + i
    return RTreeNode(node.tree, node.level + 1, child_index, _child_extents(node)[child_index])
end

function getchild(node::RTreeNode)
    extents = _child_extents(node)
    tree, childlevel = node.tree, node.level + 1
    range = _child_range(node, extents)
    return (RTreeNode(tree, childlevel, ci, @inbounds extents[ci]) for ci in range)
end

function child_indices_extents(node::RTreeNode)
    extents = _child_extents(node)
    indices = node.tree.indices
    range = _child_range(node, extents)
    return ((@inbounds(indices[i]), @inbounds(extents[i])) for i in range)
end

# The tree itself acts as the (implicit) root node.
_rootnode(tree::RTree) = RTreeNode(tree, 0, 1, tree.extent)

isleaf(tree::RTree) = length(tree.levels) == 1
nchild(tree::RTree) = length(tree.levels[1])
getchild(tree::RTree) = getchild(_rootnode(tree))
getchild(tree::RTree, i) = getchild(_rootnode(tree), i)
child_indices_extents(tree::RTree) = child_indices_extents(_rootnode(tree))

# ## Queries

"""
    query(tree::RTree, extent_or_geom)

Indices (into the collection the tree was built from) of every leaf whose
extent intersects the given extent — or the extent of the given geometry —
in ascending order.
"""
query(tree::RTree, ext::Extents.Extent) =
    sort!(depth_first_search(Base.Fix1(Extents.intersects, ext), tree))
function query(tree::RTree, geom)
    ext = GI.extent(geom)
    isnothing(ext) && throw(ArgumentError("no extent found on $(typeof(geom))"))
    return query(tree, ext)
end
