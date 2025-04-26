# # Interface implementations
# Below are some basic implementations of the interface,
# for STRTree and a "no-tree" implementation that is a flat list of extents.

using SortTileRecursiveTree: STRtree, STRNode, STRLeafNode

# ## SortTileRecursiveTree

isspatialtree(::Type{<: STRtree}) = true
node_extent(tree::STRtree) = node_extent(tree.rootnode)
nchild(tree::STRtree) = nchild(tree.rootnode)
getchild(tree::STRtree) = getchild(tree.rootnode)
getchild(tree::STRtree, i) = getchild(tree.rootnode, i)
isleaf(tree::STRtree) = isleaf(tree.rootnode)
child_indices_extents(tree::STRtree) = child_indices_extents(tree.rootnode)

isspatialtree(::Type{<: STRNode}) = true
node_extent(node::STRNode) = node.extent
nchild(node::STRNode) = length(node.children)
getchild(node::STRNode) = node.children
getchild(node::STRNode, i) = node.children[i]
isleaf(node::STRNode) = false # STRNodes are not leaves by definition

isspatialtree(::Type{<: STRLeafNode}) = true
node_extent(node::STRLeafNode) = node.extent
isleaf(node::STRLeafNode) = true
child_indices_extents(node::STRLeafNode) = zip(node.indices, node.extents)

# ## FlatNoTree
"""
    FlatNoTree(iterable_of_geoms_or_extents)

Represents a flat collection with no tree structure, i.e., a brute force search.
This is cost free, so particularly useful when you don't want to build a tree!
"""
struct FlatNoTree{T}
    geometries::T
end

isspatialtree(::Type{<: FlatNoTree}) = true
isleaf(tree::FlatNoTree) = true
node_extent(tree::FlatNoTree) = mapreduce(GI.extent, Extents.union, tree.geometries)

# NOTE: use pairs instead of enumerate here, so that we can support 
# iterators or collections that define custom `pairs` methods.
# This includes things like filtered extent lists, for example,
# so we can perform extent thinning with no allocations.
function child_indices_extents(tree::FlatNoTree{T}) where T
    # This test only applies at compile time and should be optimized away in any case.
    # And we can use multiple dispatch to override anyway, but it should be cost free I think.
    if applicable(Base.keys, T) 
        return ((i, node_extent(obj)) for (i, obj) in pairs(tree.geometries))
    else
        return ((i, node_extent(obj)) for (i, obj) in enumerate(tree.geometries))
    end
end