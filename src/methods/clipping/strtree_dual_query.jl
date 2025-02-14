"""
    STRDualQuery

A module for performing dual-tree traversals on STRtrees to find potentially overlapping geometry pairs.

The main entry point is `maybe_overlapping_geoms_and_query_lists_in_order`.
"""
module STRDualQuery

using SortTileRecursiveTree
using AbstractTrees
using SortTileRecursiveTree: STRtree, STRNode, STRLeafNode
using GeoInterface.Extents

# first define AbstractTrees interface for STRtree
include("strtree_abstracttrees.jl")

"""
    maybe_overlapping_geoms_and_query_lists_in_order(tree_a::STRtree, tree_b::STRtree)

Performs an efficient dual-tree traversal to find potentially overlapping geometry pairs.
Returns a vector of pairs, where each pair contains an index from tree_a and a sorted vector 
of indices from tree_b that might overlap.

The result looks like this:
```
[
    a1 => [b1, b2, b3],
    a2 => [b4, b5],
    a3 => [b6, b7, b8, b9],
    ...
]
```
in which the overlap map is sorted by the tree_a indices, and within each group, the tree_b indices are sorted.
"""
function maybe_overlapping_geoms_and_query_lists_in_order(tree_a::STRtree, tree_b::STRtree)
    # Use DefaultDict to automatically create empty vectors for new keys
    overlap_map = Dict{Int, Vector{Int}}()
    
    # Start the recursive traversal from the root nodes
    _dual_tree_traverse(tree_a.rootnode, tree_b.rootnode, overlap_map)
    
    # Convert to the required output format and sort
    result = [(k, sort!(v)) for (k, v) in pairs(overlap_map)]
    sort!(result, by=first)  # Sort by tree_a indices
    
    return result
end

"""
    _dual_tree_traverse(node_a::Union{STRNode,STRLeafNode}, node_b::Union{STRNode,STRLeafNode}, 
                       overlap_map::DefaultDict{Int,Vector{Int}})

Recursive helper function that performs the dual-tree traversal.
"""
function _dual_tree_traverse(node_a::Union{STRNode,STRLeafNode}, 
                           node_b::Union{STRNode,STRLeafNode}, 
                           overlap_map::Dict{Int,Vector{Int}})
    
    # Early exit if bounding boxes don't overlap
    if !Extents.intersects(nodevalue(node_a), nodevalue(node_b))
        return
    end
    
    # Case 1: Both nodes are leaves
    if node_a isa STRLeafNode && node_b isa STRLeafNode
        for idx_a in node_a.indices
            dict_vec = get!(() -> Int[], overlap_map, idx_a)
            append!(dict_vec, node_b.indices)
        end
        return
    end
    
    # Case 2: node_a is a leaf, node_b is internal
    if node_a isa STRLeafNode
        for child_b in children(node_b)
            _dual_tree_traverse(node_a, child_b, overlap_map)
        end
        return
    end
    
    # Case 3: node_b is a leaf, node_a is internal
    if node_b isa STRLeafNode
        for child_a in children(node_a)
            _dual_tree_traverse(child_a, node_b, overlap_map)
        end
        return
    end
    
    # Case 4: Both nodes are internal
    for child_a in children(node_a)
        for child_b in children(node_b)
            _dual_tree_traverse(child_a, child_b, overlap_map)
        end
    end
end

export maybe_overlapping_geoms_and_query_lists_in_order

end # module