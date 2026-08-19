"""
    dual_depth_first_search(f, predicate, tree1, tree2)

Executes a dual depth-first search over two trees, descending into the children of
nodes `i` and `j` when `predicate(node_extent(i), node_extent(j))` is true,
and pruning that branch when `predicate(node_extent(i), node_extent(j))` is false.

Finally, calls `f(i1, i2)` for each leaf-level index `i1::Int` in `tree1` and `i2::Int` in `tree2`
that satisfies `predicate(extent(i1), extent(i2))`.

Here, `f(i1::Int, i2::Int)` may be any function that takes two integers as arguments.

It may optionally return an [`Action`](@ref GeometryOps.LoopStateMachine.Action) to alter the control
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

# The descent's scratch stack of derived child extents, or `nothing` to mean "call
# `node_extent` in the loop".  `node_extent_is_expensive` is a function of `N`
# alone, so the branch folds and only one arm is compiled.
@inline _extent_stack(stack, node::N, ::E) where {N, E} = node_extent_is_expensive(N) ? _as_stack(stack, E) : nothing
@inline _as_stack(stack, ::Type) = stack
@inline _as_stack(::Nothing, ::Type{E}) where {E} = E[]

# Hand an ancestor's stack on through levels that do not need one themselves.
@inline _carry(::Nothing, stack) = stack
@inline _carry(stack, _) = stack

@inline _fill_child_extents!(::Nothing, node) = nothing
@inline function _fill_child_extents!(stack, node)
    for child in getchild(node)
        push!(stack, node_extent(child))
    end
    return nothing
end

@inline _child_extent(::Nothing, base, child, i) = node_extent(child)
@inline _child_extent(stack, base, child, i) = stack[base + i]

# `extent1` and `extent2` are a precondition, not a hint: they must equal
# `node_extent(node1)` and `node_extent(node2)`.
function dual_depth_first_search(
    f::F, predicate::P, node1::N1, extent1::E1, node2::N2, extent2::E2
) where {F, P, N1, E1, N2, E2}
    return dual_depth_first_search(f, predicate, node1, extent1, node2, extent2, nothing)
end

function dual_depth_first_search(
    f::F, predicate::P, node1::N1, extent1::E1, node2::N2, extent2::E2, stack::S
) where {F, P, N1, E1, N2, E2, S}
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
                @controlflow dual_depth_first_search(f, predicate, node1, extent1, child, child_extent, stack)
            end
        end
    elseif leaf2 # node1 is not a leaf, node2 is - recurse further into node1
        for child in getchild(node1)
            child_extent = node_extent(child)
            if predicate(child_extent, extent2)
                @controlflow dual_depth_first_search(f, predicate, child, child_extent, node2, extent2, stack)
            end
        end
    else # neither node is a leaf, recurse into both children
        # node2's child extents go on the shared stack and come off again on the
        # way out, so the descent allocates no buffer per visited node pair
        stack2 = _extent_stack(stack, node2, extent2)
        base = stack2 === nothing ? 0 : length(stack2)
        _fill_child_extents!(stack2, node2)
        child_stack = _carry(stack2, stack)
        for child1 in getchild(node1)
            child_extent1 = node_extent(child1)
            i2 = 0
            for child2 in getchild(node2)
                i2 += 1
                child_extent2 = _child_extent(stack2, base, child2, i2)
                if predicate(child_extent1, child_extent2)
                    @controlflow dual_depth_first_search(
                        f, predicate, child1, child_extent1, child2, child_extent2, child_stack
                    )
                end
            end
        end
        stack2 === nothing || resize!(stack2, base)
    end
end
