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
function sanitize_pred(pred::P) where P
    sanitize_pred(GI.trait(pred), pred)
end

sanitize_pred(::Nothing, pred::P) where P = pred
sanitize_pred(::GI.AbstractTrait, pred) = sanitize_pred(GI.extent(pred))

sanitize_pred(pred::Extents.Extent) = Base.Fix1(Extents.intersects, pred)


function query(f, index::NaturalIndex, pred)
    # pred = sanitize_pred(ipred)
    # At the top level, iterate over all toplevel nodes and query the extent
    first_level = index.levels[1]
    for (node_idx, node_extent) in enumerate(first_level.extents)
        if pred(node_extent)
            query(f, NaturalTreeNode(index, 1, node_idx, node_extent), pred)
        end
    end
end

function query(f, node::NaturalTreeNode, pred)
    # pred = sanitize_pred(ipred)
    # If we're at the bottom level, iterate over all geometries in the node
    # while performing extent rejection
    if node.level == length(node.parent_index.levels) - 1
        
        extents = node.parent_index.levels[end].extents
        # At the leaf level, indices are equivalent to the geometry index
        start_idx = (node.index - 1) * node.parent_index.nodecapacity + 1
        stop_idx = min(start_idx + node.parent_index.nodecapacity - 1, length(extents))
        
        for i in start_idx:stop_idx
            if pred(extents[i])
                # The @controlflow macro allows `f` to return an `Action`
                # from LoopStateMachine.jl, which allows us to break out of the loop, 
                # or return a value directly if we wish.
                @controlflow f(i) # provide `f` the index of the geometry in question
            end
        end
    else # not a leaf node, recurse lower down the tree

        extents = node.parent_index.levels[node.level+1].extents

        start_idx = (node.index - 1) * node.parent_index.nodecapacity + 1
        stop_idx = min(start_idx + node.parent_index.nodecapacity - 1, length(extents))
        
        for i in start_idx:stop_idx
            if pred(extents[i])
                query(f, NaturalTreeNode(node.parent_index, node.level + 1, i, extents[i]), pred)
            end
        end
    end

end