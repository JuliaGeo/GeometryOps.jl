# # Interface implementations
# Below is a "no-tree" implementation of the interface that is a flat list of extents.
# The implementation for `SortTileRecursiveTree.STRtree` lives in the package
# extension `SpatialTreeInterfaceSortTileRecursiveTreeExt`, and loads when
# SortTileRecursiveTree does.

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

# A flat collection with a known length is its own leaf: its weight is that
# length, and its indices are its own.  Sizeless iterators fall back to the
# defaults (weight 1, and a leaf walk through `child_indices_extents`).
function split_weight(tree::FlatNoTree)
    Base.IteratorSize(tree.geometries) isa Union{Base.HasLength, Base.HasShape} || return 1
    return length(tree.geometries)
end
subtree_indices(tree::FlatNoTree{<:AbstractVector}) = eachindex(tree.geometries)
