import GeoInterface as GI
import Extents

import GeometryOps.LoopStateMachine: @controlflow


struct NaturalLevel{E <: Extents.Extent}
    # level::Int      # level of node in tree
    # node_index::Int # index of the node in the level
    # extent::E       # extent of the node - precomputed and cached
    extents::Vector{E} # child extents
end

struct NaturalIndex{E <: Extents.Extent}
    nodecapacity::Int # "spread", number of children per node
    extent::E
    levels::Vector{NaturalLevel{E}}
end

GI.extent(idx::NaturalIndex) = idx.extent
Extents.extent(idx::NaturalIndex) = idx.extent

function Base.show(io::IO, ::MIME"text/plain", idx::NaturalIndex)
    println(io, "NaturalIndex with $(length(idx.levels)) levels and $(idx.nodecapacity) children per node")
    println(io, "extent: $(idx.extent)")
end

function Base.show(io::IO, idx::NaturalIndex)
    println(io, "NaturalIndex($(length(idx.levels)) levels, $(idx.extent))")
end

function NaturalIndex(geoms; nodecapacity = 32)
    e1 = GI.extent(first(geoms))
    E = typeof(e1)
    return NaturalIndex{E}(geoms; nodecapacity = nodecapacity)
end

function NaturalIndex{E}(geoms; nodecapacity = 32) where E <: Extents.Extent

    last_level_extents = GI.extent.(geoms)
    ngeoms = length(last_level_extents)
    last_level = NaturalLevel(last_level_extents)

    nlevels = _number_of_levels(nodecapacity, ngeoms)

    levels = Vector{NaturalLevel{E}}(undef, nlevels)
    levels[end] = last_level

    for level_index in (nlevels-1):(-1):1
        prev_level = levels[level_index+1] # this is always instantiated
        nrects = _number_of_keys(nodecapacity, nlevels - (level_index), ngeoms)
        # @show level_index nrects
        extents = [
            begin
                start = (rect_index - 1) * nodecapacity + 1
                stop = min(start + nodecapacity - 1, length(prev_level.extents))
                reduce(Extents.union, view(prev_level.extents, start:stop))
            end
            for rect_index in 1:nrects
        ]
        levels[level_index] = NaturalLevel(extents)
    end

    return NaturalIndex(nodecapacity, reduce(Extents.union, levels[1].extents), levels)

end

function _number_of_keys(nodecapacity::Int, level::Int, ngeoms::Int)
    return ceil(Int, ngeoms / (nodecapacity ^ (level)))
end

"""
    _number_of_levels(nodecapacity::Int, ngeoms::Int)

Calculate the number of levels in a natural tree for a given number of geometries and node capacity.

## How this works

The number of keys in a level is given by `ngeoms / nodecapacity ^ level`.

The number of levels is the smallest integer such that the number of keys in the last level is 1.
So it goes - if that makes sense.
"""
function _number_of_levels(nodecapacity::Int, ngeoms::Int)
    level = 1
    while _number_of_keys(nodecapacity, level, ngeoms) > 1
        level += 1
    end
    return level
end


# This is like a pointer to a node in the tree.
struct NaturalTreeNode{E <: Extents.Extent}
    parent_index::NaturalIndex{E}
    level::Int
    index::Int
    extent::E
end

Extents.extent(node::NaturalTreeNode) = node.extent

"""
    query(f, index::NaturalIndex, pred)

Query the index for all extents that satisfy the predicate `pred`.

`pred` must be a 1-argument function that returns a Boolean.
`pred` may also be an Extent or a geometry.

Whenever a leaf node is encountered, which 
"""
function query end

"""
    sanitize_pred(pred)

Convert a predicate to a function that returns a Boolean.

If `pred` is an Extent, convert it to a function that returns a Boolean by intersecting with the extent.
If `pred` is a geometry, convert it to an extent first, then wrap in Extents.intersects.

Otherwise, return the predicate unchanged.


Users and developers may overload this function to provide custom behaviour when something is passed in.
"""
function sanitize_query_predicate(pred::P) where P
    sanitize_query_predicate(GI.trait(pred), pred)
end

sanitize_query_predicate(::Nothing, pred::P) where P = pred
sanitize_query_predicate(::GI.AbstractTrait, pred) = sanitize_pred(GI.extent(pred))

sanitize_query_predicate(pred::Extents.Extent) = Base.Fix1(Extents.intersects, pred)

# What does SpatialTreeInterface require of trees?
# - Parents completely cover their children
# - `GI.extent(node)` returns `Extent` 
#   - can mean that `Extents.extent(node)` returns the extent of the node
# - `nchild(node)` returns the number of children of the node
# - `getchild(node)` returns an iterator over all children of the node
# - `getchild(node, i)` returns the i-th child of the node
# - `isleaf(node)` returns a boolean indicating whether the node is a leaf
# - `child_indices_extents(node)` returns an iterator over the indices and extents of the children of the node

function nchild(node::NaturalTreeNode)
    start_idx = (node.index - 1) * node.parent_index.nodecapacity + 1
    stop_idx = min(start_idx + node.parent_index.nodecapacity - 1, length(node.parent_index.levels[node.level+1].extents))
    return stop_idx - start_idx + 1
end

function getchild(node::NaturalTreeNode, i::Int)
    child_index = (node.index - 1) * node.parent_index.nodecapacity + i
    return NaturalTreeNode(
        node.parent_index, 
        node.level + 1, # increment level by 1
        child_index, # index of this particular child
        node.parent_index.levels[node.level+1].extents[child_index] # the extent of this child
    )
end

# Get all children of a node
function getchild(node::NaturalTreeNode)
    return (getchild(node, i) for i in 1:nchild(node))
end

isleaf(node::NaturalTreeNode) = node.level == length(node.parent_index.levels) - 1

function child_indices_extents(node::NaturalTreeNode)
    start_idx = (node.index - 1) * node.parent_index.nodecapacity + 1
    stop_idx = min(start_idx + node.parent_index.nodecapacity - 1, length(node.parent_index.levels[node.level+1].extents))
    return ((i, node.parent_index.levels[node.level+1].extents[i]) for i in start_idx:stop_idx)
end

# implementation for "root node" / top level tree

isleaf(node::NaturalIndex) = length(node.levels) == 1

nchild(node::NaturalIndex) = length(node.levels[1].extents)

getchild(node::NaturalIndex) = getchild(NaturalTreeNode(node, 0, 1, node.extent))
getchild(node::NaturalIndex, i) = getchild(NaturalTreeNode(node, 0, 1, node.extent), i)

child_indices_extents(node::NaturalIndex) = (i_ext for i_ext in enumerate(node.levels[1].extents))


"""
    do_query(f, predicate, tree)

Call `f(i)` for each index `i` in the tree that satisfies `predicate(extent(i))`.

This is generic to anything that implements the SpatialTreeInterface, particularly the methods
[`isleaf`](@ref), [`getchild`](@ref), and [`child_extents`](@ref).
"""
function do_query(f::F, predicate::P, node::N) where {F, P, N}
    if isleaf(node)
        for (i, leaf_geometry_extent) in child_indices_extents(node)
            if predicate(leaf_geometry_extent)
                @controlflow f(i)
            end
        end
    else
        for child in getchild(node)
            if predicate(GI.extent(child))
                do_query(f, predicate, child)
            end
        end
    end
end

function do_query(predicate, node)
    a = Int[]
    do_query(Base.Fix1(push!, a), predicate, node)
    return a
end

# implement spatial tree interface for SortTileRecursiveTree.jl

nchild(tree::STRtree) = nchild(tree.rootnode)
getchild(tree::STRtree) = getchild(tree.rootnode)
getchild(tree::STRtree, i) = getchild(tree.rootnode, i)
isleaf(tree::STRtree) = isleaf(tree.rootnode)
child_indices_extents(tree::STRtree) = child_indices_extents(tree.rootnode)


nchild(node::STRNode) = length(node.children)
getchild(node::STRNode) = node.children
getchild(node::STRNode, i) = node.children[i]
isleaf(node::STRNode) = false # STRNodes are not leaves by definition

isleaf(node::STRLeafNode) = true
child_indices_extents(node::STRLeafNode) = zip(node.indices, node.extents)

# now a `do_query` function call "just works"!


"""
    do_dual_query(f, predicate, node1, node2)

Call `f(i1, i2)` for each index `i1` in `node1` and `i2` in `node2` that satisfies `predicate(extent(i1), extent(i2))`.

This is generic to anything that implements the SpatialTreeInterface, particularly the methods
[`isleaf`](@ref), [`getchild`](@ref), and [`child_extents`](@ref).
"""
function do_dual_query(f::F, predicate::P, node1::N1, node2::N2) where {F, P, N1, N2}
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
            if predicate(GI.extent(node1), GI.extent(child))
                do_dual_query(f, predicate, node1, child)
            end
        end
    elseif isleaf(node2) # node1 is not a leaf, node2 is - recurse further into node1
        for child in getchild(node1)
            if predicate(GI.extent(child), GI.extent(node2))
                do_dual_query(f, predicate, child, node2)
            end
        end
    else # neither node is a leaf, recurse into both children
        for child1 in getchild(node1)
            for child2 in getchild(node2)
                if predicate(GI.extent(child1), GI.extent(child2))
                    do_dual_query(f, predicate, child1, child2)
                end
            end
        end
    end
end