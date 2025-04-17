
"""
    depth_first_search(f, predicate, tree)

Call `f(i)` for each index `i` in the tree that satisfies `predicate(extent(i))`.

This is generic to anything that implements the SpatialTreeInterface, particularly the methods
[`isleaf`](@ref), [`getchild`](@ref), and [`child_extents`](@ref).
"""
function depth_first_search(f::F, predicate::P, node::N) where {F, P, N}
    if isleaf(node)
        for (i, leaf_geometry_extent) in child_indices_extents(node)
            if predicate(leaf_geometry_extent)
                @controlflow f(i)
            end
        end
    else
        for child in getchild(node)
            if predicate(node_extent(child))
                @controlflow depth_first_search(f, predicate, child)
            end
        end
    end
end
function depth_first_search(predicate, node)
    a = Int[]
    depth_first_search(Base.Fix1(push!, a), predicate, node)
    return a
end
