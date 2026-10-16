# # SpatialTreeInterface

using ..SpatialTreeInterface
import ..SpatialTreeInterface: spatialtree, isspatialtree, isleaf, nchild, getchild,
    child_indices_extents, split_weight, subtree_indices

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

# The leaf slots under a node form one contiguous window.  A node at level `l`
# spans `nodecapacity^(nlev - l)` slots, since the leaf entries sit at level
# `nlev`; the last node of a level is clipped to the leaf count.
@inline function _leaf_range(node::RTreeNode)
    span = node.tree.nodecapacity^(length(node.tree.levels) - node.level)
    lo = (node.index - 1) * span + 1
    hi = min(node.index * span, length(node.tree.levels[end]))
    return lo:hi
end

split_weight(node::RTreeNode) = length(_leaf_range(node))
subtree_indices(node::RTreeNode) = view(node.tree.indices, _leaf_range(node))

# The tree itself acts as the (implicit) root node.
_rootnode(tree::RTree) = RTreeNode(tree, 0, 1, tree.extent)

isleaf(tree::RTree) = length(tree.levels) == 1
nchild(tree::RTree) = length(tree.levels[1])
getchild(tree::RTree) = getchild(_rootnode(tree))
getchild(tree::RTree, i) = getchild(_rootnode(tree), i)
child_indices_extents(tree::RTree) = child_indices_extents(_rootnode(tree))
split_weight(tree::RTree) = length(tree.indices)
subtree_indices(tree::RTree) = tree.indices

# ## Queries
# `SpatialTreeInterface.query(tree, extent_or_geom_or_predicate)` works on an
# `RTree` as it does on every other tree; there is no `RTree`-specific query.

"""
    spatialtree([manifold::Manifold], geometries)

Build a default packed STR [`RTree`](@ref) from an iterable of geometries, over
their extents on `manifold` (`Planar()` if not given).  This is GeometryOps'
method for `SpatialTreeInterface.spatialtree`.

`missing` and `nothing` entries are skipped, as are empty geometries.  Queries
return indices into `geometries` regardless, so the skipped entries do not
shift the answers.  If no valid geometries remain, return `nothing`.
"""
spatialtree(geometries) = spatialtree(Planar(), geometries)
function spatialtree(manifold::Manifold, geometries)
    items = geometries isa AbstractVector ? geometries : collect(geometries)
    positions = Int[]
    for (i, geom) in pairs(items)
        if (isnothing(geom) || ismissing(geom) || GI.isempty(geom))
            continue
        end
        push!(positions, i)
    end
    isempty(positions) && return nothing
    return RTree(manifold, STR(), items; indices = positions)
end
