# # SortTileRecursiveTree
# `STRtree` and its nodes implement the interface directly: the tree forwards
# to its root node, inner nodes expose their children, and leaf nodes zip
# their stored indices with their extents.
module SpatialTreeInterfaceSortTileRecursiveTreeExt

import SpatialTreeInterface: isspatialtree, isleaf, getchild, nchild, child_indices_extents, node_extent
using SortTileRecursiveTree: STRtree, STRNode, STRLeafNode

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

end # module
