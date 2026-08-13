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
[`isleaf`](@ref), [`getchild`](@ref), [`node_extent`](@ref) and [`child_indices_extents`](@ref).

Each visited node's extent is computed once and carried into the recursion, so
implementations are free to make `node_extent` as precise as they like without
paying for it again at every level of a descent.  Trees that *derive* their
extents rather than storing them should also define
[`node_extent_is_expensive`](@ref) - see there for what that buys.

## Examples

```julia
using NaturalEarth,
```
"""
function dual_depth_first_search(f::F, predicate::P, node1::N1, node2::N2) where {F, P, N1, N2}
    return dual_depth_first_search(f, predicate, node1, node_extent(node1), node2, node_extent(node2))
end

# Children paired with their extents, for a loop that makes a single pass.  Lazy,
# so this costs exactly one `node_extent` per child and allocates nothing.
@inline _children_with_extents(node) = ((child, node_extent(child)) for child in getchild(node))

# Children paired with their extents, for the *inner* loop of the both-internal
# branch below, which is re-entered once per child of the opposing node.
#
# The lazy form re-derives every child's extent on each pass, which is what the
# loops here used to do inline.  For a tree that stores an extent per node that
# is a field load and re-deriving is free, so we keep the generator and allocate
# nothing.  For a tree that computes its extents, `node_extent_is_expensive`
# opts into materializing them once per node pair instead - trading a small
# vector for a factor-of-fanout reduction in `node_extent` calls.
#
# The branch is on a trait that is constant per type, so exactly one of these is
# compiled into any given traversal.
@inline _reusable_children_with_extents(node) =
    _reusable_children_with_extents(Val(node_extent_is_expensive(node)), node)
@inline _reusable_children_with_extents(::Val{false}, node) = _children_with_extents(node)
@inline _reusable_children_with_extents(::Val{true}, node) = collect(_children_with_extents(node))

# The extent arguments are a precondition, not a hint: `extent1` must equal
# `node_extent(node1)` and `extent2` must equal `node_extent(node2)`.  Callers
# that already hold a node's extent can use this form to avoid recomputing it.
function dual_depth_first_search(
    f::F, predicate::P, node1::N1, extent1::E1, node2::N2, extent2::E2
) where {F, P, N1, E1, N2, E2}
    leaf1 = isleaf(node1)
    leaf2 = isleaf(node2)
    if leaf1 && leaf2
        # both nodes are leaves, so we can just iterate over the indices and extents.
        # `child_indices_extents` is bound once here rather than being re-evaluated
        # per outer iteration - implementations are encouraged to materialize a
        # vector there, and that would otherwise be rebuilt for every cell in node1.
        inner = child_indices_extents(node2)
        for (i1, cell_extent1) in child_indices_extents(node1)
            for (i2, cell_extent2) in inner
                if predicate(cell_extent1, cell_extent2)
                    @controlflow f(i1, i2)
                end
            end
        end
    elseif leaf1 # node2 is not a leaf, node1 is - recurse further into node2
        for (child, child_extent) in _children_with_extents(node2)
            if predicate(extent1, child_extent)
                @controlflow dual_depth_first_search(f, predicate, node1, extent1, child, child_extent)
            end
        end
    elseif leaf2 # node1 is not a leaf, node2 is - recurse further into node1
        for (child, child_extent) in _children_with_extents(node1)
            if predicate(child_extent, extent2)
                @controlflow dual_depth_first_search(f, predicate, child, child_extent, node2, extent2)
            end
        end
    else # neither node is a leaf, recurse into both children
        children2 = _reusable_children_with_extents(node2)
        for (child1, child_extent1) in _children_with_extents(node1)
            for (child2, child_extent2) in children2
                if predicate(child_extent1, child_extent2)
                    @controlflow dual_depth_first_search(
                        f, predicate, child1, child_extent1, child2, child_extent2
                    )
                end
            end
        end
    end
end
