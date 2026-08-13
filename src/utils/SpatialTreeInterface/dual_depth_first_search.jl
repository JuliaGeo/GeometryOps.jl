"""
    node_size(node, extent)

How much of the search space a node covers, used by
[`dual_depth_first_search_balanced`](@ref) to pick which of two internal nodes to
descend.  Any measure monotone in the size of the node will do - only the ratio
between the two matters.

Defaults to the area (or volume) of an `Extents.Extent`; trees whose
`node_extent` returns anything else must overload it.  It takes the node as well
as the extent so that a tree which knows how many cells sit beneath a node can
measure that instead, e.g. `node_size(n::MyCursor, extent) = prod(ncells(n))`.
"""
function node_size end

node_size(node, extent::Extents.Extent) = prod(map(b -> b[2] - b[1], values(extent)))
node_size(node, extent) = throw(ArgumentError("""
    `node_size` is not defined for extents of type $(typeof(extent)).

    `dual_depth_first_search_balanced` needs to compare how much space two nodes
    cover.  The default only knows how to measure an `Extents.Extent`; define

        SpatialTreeInterface.node_size(node::$(typeof(node)), extent) = ...

    returning any real number that is monotone in the size of the node.
    """))

# How to descend when both nodes of a pair are internal.  These are internal:
# the two public entry points below select between them.
struct Lockstep end

struct Balanced{T}
    threshold::T
end

# Which side of an internal/internal pair to split.  `Lockstep` returns a
# constant, so the two single-sided branches below are eliminated at compile
# time and the default traversal is exactly what it was.
@inline _descent_side(::Lockstep, node1, extent1, node2, extent2) = :both
@inline function _descent_side(descent::Balanced, node1, extent1, node2, extent2)
    size1 = node_size(node1, extent1)
    size2 = node_size(node2, extent2)
    size1 > descent.threshold * size2 && return :first
    size2 > descent.threshold * size1 && return :second
    return :both
end

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

Each visited node's extent is computed once and carried into the recursion.
Trees that derive their extents rather than storing them should also define
[`node_extent_is_expensive`](@ref).

When both nodes of a pair are internal, this descends both at once.  For trees
that are structurally mismatched, see
[`dual_depth_first_search_balanced`](@ref).

## Examples

```julia
using NaturalEarth,
```
"""
function dual_depth_first_search(f::F, predicate::P, node1::N1, node2::N2) where {F, P, N1, N2}
    return dual_depth_first_search(f, predicate, node1, node_extent(node1), node2, node_extent(node2))
end

"""
    dual_depth_first_search_balanced(f, predicate, tree1, tree2; threshold = 4)

As [`dual_depth_first_search`](@ref), but when both nodes of a pair are internal
and one covers more than `threshold` times as much space as the other, descend
only the larger and carry the smaller through unchanged.  This closes the gap
between trees of different depth and branching factor, which lockstep descent
cannot - it halves both sides at every step, so one side reaches its leaves while
the other is still many levels too coarse.

`threshold` is a ratio of [`node_size`](@ref).  `1` always splits the larger
node; a very large threshold recovers `dual_depth_first_search`.

Prunes on the same condition as `dual_depth_first_search` - a failed `predicate`
at a node pair - so on any tree where a parent's extent covers its children's,
both reach the same leaf pairs.  It is not a general improvement: on well-matched
trees it visits more node pairs to get there.  Measure before switching.
"""
function dual_depth_first_search_balanced(
    f::F, predicate::P, node1::N1, node2::N2; threshold = 4
) where {F, P, N1, N2}
    return dual_depth_first_search_balanced(
        f, predicate, node1, node_extent(node1), node2, node_extent(node2); threshold
    )
end

function dual_depth_first_search_balanced(
    f::F, predicate::P, node1::N1, extent1::E1, node2::N2, extent2::E2; threshold = 4
) where {F, P, N1, E1, N2, E2}
    threshold >= 1 || throw(ArgumentError("`threshold` must be at least 1, got $threshold"))
    return _dual_depth_first_search(
        f, predicate, node1, extent1, node2, extent2, Balanced(threshold)
    )
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
    return _dual_depth_first_search(f, predicate, node1, extent1, node2, extent2, Lockstep())
end

function _dual_depth_first_search(
    f::F, predicate::P, node1::N1, extent1::E1, node2::N2, extent2::E2, descent::D
) where {F, P, N1, E1, N2, E2, D}
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
                @controlflow _dual_depth_first_search(
                    f, predicate, node1, extent1, child, child_extent, descent
                )
            end
        end
    elseif leaf2 # node1 is not a leaf, node2 is - recurse further into node1
        for child in getchild(node1)
            child_extent = node_extent(child)
            if predicate(child_extent, extent2)
                @controlflow _dual_depth_first_search(
                    f, predicate, child, child_extent, node2, extent2, descent
                )
            end
        end
    else # neither node is a leaf
        side = _descent_side(descent, node1, extent1, node2, extent2)
        if side === :both # recurse into both children
            extents2 = _child_extents(node2)
            for child1 in getchild(node1)
                child_extent1 = node_extent(child1)
                i2 = 0
                for child2 in getchild(node2)
                    i2 += 1
                    child_extent2 = _child_extent(extents2, child2, i2)
                    if predicate(child_extent1, child_extent2)
                        @controlflow _dual_depth_first_search(
                            f, predicate, child1, child_extent1, child2, child_extent2, descent
                        )
                    end
                end
            end
        elseif side === :first # node1 is much larger - split it, carry node2 through
            for child1 in getchild(node1)
                child_extent1 = node_extent(child1)
                if predicate(child_extent1, extent2)
                    @controlflow _dual_depth_first_search(
                        f, predicate, child1, child_extent1, node2, extent2, descent
                    )
                end
            end
        else # node2 is much larger - split it, carry node1 through
            for child2 in getchild(node2)
                child_extent2 = node_extent(child2)
                if predicate(extent1, child_extent2)
                    @controlflow _dual_depth_first_search(
                        f, predicate, node1, extent1, child2, child_extent2, descent
                    )
                end
            end
        end
    end
end
