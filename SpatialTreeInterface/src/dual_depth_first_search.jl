"""
    dual_depth_first_search(f, predicate, tree1, tree2)

Executes a dual depth-first search over two trees, descending into the children of
nodes `i` and `j` when `predicate(node_extent(i), node_extent(j))` is true,
and pruning that branch when `predicate(node_extent(i), node_extent(j))` is false.

Finally, calls `f(i1, i2)` for each leaf-level index `i1::Int` in `tree1` and `i2::Int` in `tree2`
that satisfies `predicate(extent(i1), extent(i2))`.

Here, `f(i1::Int, i2::Int)` may be any function that takes two integers as arguments.

It may optionally return an [`Action`](@ref LoopStateMachine.Action) to alter the control
flow of the `Action(:full_return, true)`. Return `Action(:full_return, true)` from this
function and break out of the recursion.

This is generic to anything that implements the SpatialTreeInterface, particularly the methods
[`isleaf`](@ref), [`getchild`](@ref), [`node_extent`](@ref) and [`child_indices_extents`](@ref).

Each visited node's extent is computed once and carried into the recursion.
Trees that derive their extents rather than storing them should also define
[`node_extent_is_expensive`](@ref).
"""
function dual_depth_first_search(f::F, predicate::P, node1::N1, node2::N2) where {F, P, N1, N2}
    return dual_depth_first_search(f, predicate, node1, node_extent(node1), node2, node_extent(node2))
end

# Extents of `node`'s children, derived once, or `nothing` to mean "call
# `node_extent` in the loop".  `node_extent_is_expensive` is a function of `N`
# alone, so the branch folds and only one arm is compiled.
@inline function _child_extents(node::N) where {N}
    node_extent_is_expensive(N) || return nothing
    return [node_extent(child) for child in getchild(node)]
end

@inline _child_extent(::Nothing, child, i) = node_extent(child)
@inline _child_extent(extents, child, i) = extents[i]

# `extent1` and `extent2` are a precondition, not a hint: they must equal
# `node_extent(node1)` and `node_extent(node2)`.
function dual_depth_first_search(
    f::F, predicate::P, node1::N1, extent1::E1, node2::N2, extent2::E2
) where {F, P, N1, E1, N2, E2}
    leaf1 = isleaf(node1)
    leaf2 = isleaf(node2)
    if leaf1 && leaf2
        # bound once each - `cie_2` would otherwise be rebuilt per cell of node1
        cie_1 = child_indices_extents(node1)
        cie_2 = child_indices_extents(node2)
        for (i1, cell_extent1) in cie_1
            for (i2, cell_extent2) in cie_2
                if predicate(cell_extent1, cell_extent2)
                    @controlflow f(i1, i2)
                end
            end
        end
    elseif leaf1 # node2 is not a leaf, node1 is - recurse further into node2
        for child in getchild(node2)
            child_extent = node_extent(child)
            if predicate(extent1, child_extent)
                @controlflow dual_depth_first_search(f, predicate, node1, extent1, child, child_extent)
            end
        end
    elseif leaf2 # node1 is not a leaf, node2 is - recurse further into node1
        for child in getchild(node1)
            child_extent = node_extent(child)
            if predicate(child_extent, extent2)
                @controlflow dual_depth_first_search(f, predicate, child, child_extent, node2, extent2)
            end
        end
    else # neither node is a leaf, recurse into both children
        extents2 = _child_extents(node2)
        for child1 in getchild(node1)
            child_extent1 = node_extent(child1)
            i2 = 0
            for child2 in getchild(node2)
                i2 += 1
                child_extent2 = _child_extent(extents2, child2, i2)
                if predicate(child_extent1, child_extent2)
                    @controlflow dual_depth_first_search(
                        f, predicate, child1, child_extent1, child2, child_extent2
                    )
                end
            end
        end
    end
end
