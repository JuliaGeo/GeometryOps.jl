"""
    dual_depth_first_search(f, predicate, tree1, tree2)

Executes a dual depth-first search over two trees, descending into the children of 
nodes `i` and `j` when `predicate(node_extent(i), node_extent(j))` is true, 
and pruning that branch when `predicate(node_extent(i), node_extent(j))` is false.

Finally, calls `f(i1, i2)` for each leaf-level index `i1::Int` in `tree1` and `i2::Int` in `tree2` 
that satisfies `predicate(extent(i1), extent(i2))`.

Here, `f(i1::Int, i2::Int)` may be any function that takes two integers as arguments.
It may optionally return an [`Action`](@ref LoopStateMachine.Action) to alter the control
flow of the `Action(:full_return, true)` to return `Action(:full_return, true)` from this 
function and break out of the recursion.

This is generic to anything that implements the SpatialTreeInterface, particularly the methods
[`isleaf`](@ref), [`getchild`](@ref), and [`child_indices_extents`](@ref).

## Examples

```julia
using NaturalEarth, 
```
"""
function dual_depth_first_search(f::F, predicate::P, node1::N1, node2::N2) where {F, P, N1, N2}
    if isleaf(node1) && isleaf(node2)
        # both nodes are leaves, so we can just iterate over the indices and extents
        for (i1, extent1) in child_indices_extents(node1)
            for (i2, extent2) in child_indices_extents(node2)
                if predicate(extent1, extent2)
                    @controlflow f(i1, i2)
                end
            end
        end
    elseif isleaf(node1) # node2 is not a leaf, node1 is - recurse further into node2
        for child in getchild(node2)
            if predicate(node_extent(node1), node_extent(child))
                @controlflow dual_depth_first_search(f, predicate, node1, child)
            end
        end
    elseif isleaf(node2) # node1 is not a leaf, node2 is - recurse further into node1
        for child in getchild(node1)
            if predicate(node_extent(child), node_extent(node2))
                @controlflow dual_depth_first_search(f, predicate, child, node2)
            end
        end
    else # neither node is a leaf, recurse into both children
        for child1 in getchild(node1)
            for child2 in getchild(node2)
                if predicate(node_extent(child1), node_extent(child2))
                    @controlflow dual_depth_first_search(f, predicate, child1, child2)
                end
            end
        end
    end
end
