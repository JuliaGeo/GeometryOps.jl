"""
    STRDualQuery

A module for performing dual-tree traversals on STRtrees to find potentially overlapping geometry pairs.

The main entry point is `maybe_overlapping_geoms_and_query_lists_in_order`.
"""
module STRDualQuery

using SortTileRecursiveTree

using SortTileRecursiveTree: STRtree, STRNode, STRLeafNode
using GeoInterface.Extents

import GeoInterface as GI

"helper function to get the extent of any STR node, since leaf nodes don't store global extent."

"""
    maybe_overlapping_geoms_and_query_lists_in_order(tree_a::STRtree, tree_b::STRtree, edges_a::Vector{<: GI.Line}, edges_b::Vector{<: GI.Line})

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
function maybe_overlapping_geoms_and_query_lists_in_order(tree_a::STRtree, tree_b::STRtree, edges_a::Vector{<: GI.Line}, edges_b::Vector{<: GI.Line})
    # Use DefaultDict to automatically create empty vectors for new keys
    overlap_map = Dict{Int, Vector{Int}}()
    
    # Start the recursive traversal from the root nodes
    _dual_tree_traverse!(overlap_map, tree_a.rootnode, tree_b.rootnode, edges_a, edges_b)
    
    # Convert to the required output format and sort
    result = [(k, sort!(v)) for (k, v) in pairs(overlap_map)]
    sort!(result, by=first)  # Sort by tree_a indices
    
    return result
end

"""
    _dual_tree_traverse!(
        overlap_map::Dict{Int,Vector{Int}}, 
        node_a::Union{STRNode,STRLeafNode}, node_b::Union{STRNode,STRLeafNode}, 
        edges_a::Vector{<: GI.Line}, edges_b::Vector{<: GI.Line}
    )

Recursive helper function that performs the dual-tree traversal and stores results in `overlap_map`.
"""
function _dual_tree_traverse!(
    overlap_map::Dict{Int,Vector{Int}}, 
    node_a::Union{STRNode,STRLeafNode}, 
    node_b::Union{STRNode,STRLeafNode}, 
    edges_a::Vector{<: GI.Line}, 
    edges_b::Vector{<: GI.Line}
)
    
    # Early exit if bounding boxes don't overlap
    if !Extents.intersects(GI.extent(node_a), GI.extent(node_b))
        return
    end
    
    # Case 1: Both nodes are leaves
    if node_a isa STRLeafNode && node_b isa STRLeafNode
        for idx_a in node_a.indices
            dict_vec = get!(() -> Int[], overlap_map, idx_a)
            for idx_b in node_b.indices
                # Final extent rejection, this is cheaper than the allocation for `push!`
                if Extents.intersects(GI.extent(edges_a[idx_a]), GI.extent(edges_b[idx_b]))
                    push!(dict_vec, idx_b)
                end
            end
        end
        return
    end
    
    # Case 2: node_a is a leaf, node_b is internal
    if node_a isa STRLeafNode
        for child_b in node_b.children
            _dual_tree_traverse!(overlap_map, node_a, child_b, edges_a, edges_b)
        end
        return
    end
    
    # Case 3: node_b is a leaf, node_a is internal
    if node_b isa STRLeafNode
        for child_a in node_a.children
            _dual_tree_traverse!(overlap_map, child_a, node_b, edges_a, edges_b)
        end
        return
    end
    
    # Case 4: Both nodes are internal
    for child_a in node_a.children
        for child_b in node_b.children
            _dual_tree_traverse!(overlap_map, child_a, child_b, edges_a, edges_b)
        end
    end
end

export maybe_overlapping_geoms_and_query_lists_in_order

end # module

# The code to use this is:
# elseif accelerator isa DoubleSTRtree
#     # If both of the polygons are quite large,
#     # then we do a dual-tree traversal of the STRtrees
#     # and find all potentially overlapping edges.
#     # This is kind of like an adjacency list of a graph or a sparse matrix that we're constructing.

#     # First, we materialize an edge list for each polygon as a vector of `GI.Line` objects.
#     ext_a = GI.extent(poly_a)
#     ext_b = GI.extent(poly_b)
#     edges_a, indices_a = to_edgelist(ext_b, poly_a, T) # ::Vector{GI.Line} # with precalculated extents
#     edges_b, indices_b = to_edgelist(ext_a, poly_b, T) # ::Vector{GI.Line} # with precalculated extents

#     # Now, construct STRtrees from the edge lists.
#     # TODO: we can optimize the strtrees by passing in only edges that reside within the area of extent overlap
#     # between poly_a and poly_b.
#     # this would greatly help for e.g. coverage transactions
#     # BUT this needs fixes in SortTileRecursiveTree.jl to allow you to pass in a vector of indices
#     # as well as extents / geometries.
#     tree_a = STRtree(edges_a)
#     tree_b = STRtree(edges_b)

#     # Now we perform a dual-tree traversal to find 
#     # all potentially overlapping pairs of edges.
#     # The `STRDualQuery` module is defined in this folder,
#     # in `strtree_dual_query.jl`.
#     index_list = STRDualQuery.maybe_overlapping_geoms_and_query_lists_in_order(
#         tree_a, tree_b, edges_a, edges_b
#     )

#     # Track the last index in `poly_a` that we've processed.
#     last_a_idx_orig = indices_a[1]

#     # Iterate over the indices in `poly_a` that may potentially overlap with any edge in `poly_b`.
#     for (a_idx, b_idxs) in index_list
#         a_idx_orig = indices_a[a_idx]
#         # If we've skipped any indices in `poly_a`,
#         # then we need to process those indices.
#         # But we know that they will never intersect, so we can skip calling `f_intersect`.
#         if a_idx_orig > (last_a_idx_orig + 1)
#             for i in (last_a_idx_orig + 1:a_idx_orig - 1)
#                 # TODO explain this check
#                 # a1, a2 = edges_a[i].geom
#                 # a1 == a2 && continue
#                 point_a = GI.getpoint(poly_a, i)
#                 f_on_each_a(point_a, i)
#                 f_after_each_a(point_a, i)
#             end
#         end
#         last_a_idx_orig = a_idx_orig
        
#         # Get the edge in `poly_a` that we're currently processing.
#         aedge = edges_a[a_idx]
#         aedge_extent = GI.extent(aedge)
#         a1t, a2t = aedge.geom

#         # Call `f_a` for the start point of the edge.
#         f_on_each_a(a1t, a_idx_orig)

#         for b_idx in b_idxs
#             b_idx_orig = indices_b[b_idx]
#             bedge = edges_b[b_idx]
#             b1t, b2t = bedge.geom
#             b1t == b2t && continue

#             if GI.Extents.intersects(aedge_extent, GI.extent(bedge))
#                 LoopStateMachine.@controlflow f_on_each_maybe_intersect(
#                     ((a1t, a2t), a_idx_orig), ((b1t, b2t), b_idx_orig)
#                 )
#             end
#         end

#         # Call `f_after_each_a` for the start point of the edge - postprocess once
#         # we're done with tracing.
#         f_after_each_a(a1t, a_idx_orig)

#         last_a_idx_orig = a_idx_orig

#     end

#     # println("DOUBLE TREE TRAVERSAL")
# end

# using Test
# using GeometryOps
# using GeometryOps.STRDualQuery
# using SortTileRecursiveTree
# using StaticArrays
# using GI.Extents

# @testset "STRDualQuery" begin
#     @testset "Basic overlapping rectangles" begin
#         # Create two sets of rectangles represented by their corner points
#         # Tree A:                Tree B:
#         # [0,0]--[1,1]          [0.5,0.5]--[1.5,1.5]    (overlaps with A1)
#         # [2,2]--[3,3]          [2.5,2.5]--[3.5,3.5]    (overlaps with A2)
        
#         edges_a = [
#             ((0.0, 0.0), (1.0, 1.0)),  # A1
#             ((2.0, 2.0), (3.0, 3.0))   # A2
#         ]
#         edges_b = [
#             ((0.5, 0.5), (1.5, 1.5)),  # B1
#             ((2.5, 2.5), (3.5, 3.5))   # B2
#         ]

#         # Convert edges to STRtree format
#         tree_a = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_a])
#         tree_b = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_b])

#         result = maybe_overlapping_geoms_and_query_lists_in_order(tree_a, tree_b)
        
#         # Check results
#         @test length(result) == 2
#         @test result[1][1] == 1  # First edge from tree_a
#         @test result[1][2] == [1]  # Overlaps with first edge from tree_b
#         @test result[2][1] == 2  # Second edge from tree_a
#         @test result[2][2] == [2]  # Overlaps with second edge from tree_b
#     end

#     @testset "Non-overlapping geometries" begin
#         edges_a = [((0.0, 0.0), (1.0, 1.0))]
#         edges_b = [((10.0, 10.0), (11.0, 11.0))]

#         tree_a = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_a])
#         tree_b = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_b])

#         result = maybe_overlapping_geoms_and_query_lists_in_order(tree_a, tree_b)
#         @test isempty(result)
#     end

#     @testset "Multiple overlaps" begin
#         # Create a scenario where one edge overlaps with multiple others
#         edges_a = [((0.0, 0.0), (2.0, 2.0))]  # One long diagonal line
#         edges_b = [
#             ((0.5, 0.5), (1.0, 1.0)),  # B1 overlaps
#             ((1.0, 1.0), (1.5, 1.5)),  # B2 overlaps
#             ((1.5, 1.5), (2.0, 2.0)),  # B3 overlaps
#             ((3.0, 3.0), (4.0, 4.0))   # B4 doesn't overlap
#         ]

#         tree_a = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_a])
#         tree_b = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_b])

#         result = maybe_overlapping_geoms_and_query_lists_in_order(tree_a, tree_b)
        
#         @test length(result) == 1
#         @test result[1][1] == 1
#         @test result[1][2] == [1, 2, 3]  # Should find first three edges from tree_b
#     end

#     @testset "Empty trees" begin
#         empty_tree = STRtree(GI.Line{2, Float64}[])
#         edges_a = [((0.0, 0.0), (1.0, 1.0))]
#         non_empty_tree = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_a])

#         # Test empty tree with non-empty tree
#         result1 = maybe_overlapping_geoms_and_query_lists_in_order(empty_tree, non_empty_tree)
#         @test isempty(result1)

#         # Test non-empty tree with empty tree
#         result2 = maybe_overlapping_geoms_and_query_lists_in_order(non_empty_tree, empty_tree)
#         @test isempty(result2)

#         # Test empty tree with empty tree
#         result3 = maybe_overlapping_geoms_and_query_lists_in_order(empty_tree, empty_tree)
#         @test isempty(result3)
#     end
# end