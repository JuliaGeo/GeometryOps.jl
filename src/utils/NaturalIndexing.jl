module NaturalIndexing

import GeoInterface as GI
import Extents

using ..SpatialTreeInterface

import ..GeometryOps as GO # TODO: only needed for NaturallyIndexedRing, remove when that is removed.

export NaturalIndex, NaturallyIndexedRing, prepare_naturally

"""
    NaturalLevel{E <: Extents.Extent}

A level in the natural tree.  Stored in a vector in [`NaturalIndex`](@ref).

- `extents` is a vector of extents of the children of the node
"""
struct NaturalLevel{E <: Extents.Extent}
    extents::Vector{E} # child extents
end

Base.show(io::IO, level::NaturalLevel) = print(io, "NaturalLevel($(length(level.extents)) extents)")
Base.show(io::IO, ::MIME"text/plain", level::NaturalLevel) = Base.show(io, level)

"""
    NaturalIndex{E <: Extents.Extent}

A natural tree index.  Stored in a vector in [`NaturalIndex`](@ref).

- `nodecapacity` is the "spread", number of children per node
- `extent` is the extent of the tree
- `levels` is a vector of [`NaturalLevel`](@ref)s
"""
struct NaturalIndex{E <: Extents.Extent}
    nodecapacity::Int # "spread", number of children per node
    extent::E
    levels::Vector{NaturalLevel{E}}
end

Extents.extent(idx::NaturalIndex) = idx.extent

function Base.show(io::IO, ::MIME"text/plain", idx::NaturalIndex)
    println(io, "NaturalIndex with $(length(idx.levels)) levels and $(idx.nodecapacity) children per node")
    println(io, "extent: $(idx.extent)")
end
function Base.show(io::IO, idx::NaturalIndex)
    println(io, "NaturalIndex($(length(idx.levels)) levels, $(idx.extent))")
end

function NaturalIndex(geoms; nodecapacity = 32)
    # Get the extent type initially (coord order, coord type, etc.)
    # so that the construction is type stable.
    e1 = GI.extent(first(geoms))
    E = typeof(e1)
    return NaturalIndex{E}(geoms; nodecapacity = nodecapacity)
end
function NaturalIndex(last_level_extents::Vector{E}; nodecapacity = 32) where E <: Extents.Extent
    # If we are passed a vector of extents - inflate immediately!
    return NaturalIndex{E}(last_level_extents; nodecapacity = nodecapacity)
end

function NaturalIndex{E}(geoms; nodecapacity = 32) where E <: Extents.Extent
    # If passed a vector of geometries, and we know the type of the extent,
    # then simply retrieve the extents so they can serve as the "last-level" 
    # extents.
    # Finally, call the lowest level method that performs inflation.
    last_level_extents = GI.extent.(geoms)
    return NaturalIndex{E}(last_level_extents; nodecapacity = nodecapacity)
end
# This is the main constructor that performs inflation.
function NaturalIndex{E}(last_level_extents::Vector{E}; nodecapacity = 32) where E <: Extents.Extent
    ngeoms = length(last_level_extents)
    last_level = NaturalLevel(last_level_extents)

    nlevels = _number_of_levels(nodecapacity, ngeoms)

    levels = Vector{NaturalLevel{E}}(undef, nlevels)
    levels[end] = last_level
    # Iterate backwards, from bottom to top level,
    # and build up the level extent vectors.
    for level_index in (nlevels-1):(-1):1
        prev_level = levels[level_index+1] # this is always instantiated, since we are iterating backwards
        nrects = _number_of_keys(nodecapacity, nlevels - (level_index), ngeoms)
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
"""
    NaturalIndexNode{E <: Extents.Extent}

A reference to a node in the natural tree.  Kind of like a tree cursor.

- `parent_index` is a pointer to the parent index
- `level` is the level of the node in the tree
- `index` is the index of the node in the level
- `extent` is the extent of the node
"""
struct NaturalIndexNode{E <: Extents.Extent}
    parent_index::NaturalIndex{E}
    level::Int
    index::Int
    extent::E
end

Extents.extent(node::NaturalIndexNode) = node.extent

# What does SpatialTreeInterface require of trees?
# - Parents completely cover their children
# - `GI.extent(node)` returns `Extent` 
#   - can mean that `Extents.extent(node)` returns the extent of the node
# - `nchild(node)` returns the number of children of the node
# - `getchild(node)` returns an iterator over all children of the node
# - `getchild(node, i)` returns the i-th child of the node
# - `isleaf(node)` returns a boolean indicating whether the node is a leaf
# - `child_indices_extents(node)` returns an iterator over the indices and extents of the children of the node

SpatialTreeInterface.isspatialtree(::Type{<: NaturalIndex}) = true
SpatialTreeInterface.isspatialtree(::Type{<: NaturalIndexNode}) = true

function SpatialTreeInterface.nchild(node::NaturalIndexNode)
    start_idx = (node.index - 1) * node.parent_index.nodecapacity + 1
    stop_idx = min(start_idx + node.parent_index.nodecapacity - 1, length(node.parent_index.levels[node.level+1].extents))
    return stop_idx - start_idx + 1
end

function SpatialTreeInterface.getchild(node::NaturalIndexNode, i::Int)
    child_index = (node.index - 1) * node.parent_index.nodecapacity + i
    return NaturalIndexNode(
        node.parent_index, 
        node.level + 1, # increment level by 1
        child_index, # index of this particular child
        node.parent_index.levels[node.level+1].extents[child_index] # the extent of this child
    )
end

# Get all children of a node
function SpatialTreeInterface.getchild(node::NaturalIndexNode)
    return (SpatialTreeInterface.getchild(node, i) for i in 1:SpatialTreeInterface.nchild(node))
end

SpatialTreeInterface.isleaf(node::NaturalIndexNode) = node.level == length(node.parent_index.levels) - 1

function SpatialTreeInterface.child_indices_extents(node::NaturalIndexNode)
    start_idx = (node.index - 1) * node.parent_index.nodecapacity + 1
    stop_idx = min(start_idx + node.parent_index.nodecapacity - 1, length(node.parent_index.levels[node.level+1].extents))
    return ((i, node.parent_index.levels[node.level+1].extents[i]) for i in start_idx:stop_idx)
end

# implementation for "root node" / top level tree

SpatialTreeInterface.isleaf(node::NaturalIndex) = length(node.levels) == 1

SpatialTreeInterface.nchild(node::NaturalIndex) = length(node.levels[1].extents)

SpatialTreeInterface.getchild(node::NaturalIndex) = SpatialTreeInterface.getchild(NaturalIndexNode(node, 0, 1, node.extent))
SpatialTreeInterface.getchild(node::NaturalIndex, i) = SpatialTreeInterface.getchild(NaturalIndexNode(node, 0, 1, node.extent), i)

SpatialTreeInterface.child_indices_extents(node::NaturalIndex) = (i_ext for i_ext in enumerate(node.levels[1].extents))

"""
    NaturallyIndexedRing(points; nodecapacity = 32)

A linear ring that contains a natural index.

!!! warning
    This will be removed in favour of prepared geometry - the idea here
    is just to test what interface works best to store things in.
"""
struct NaturallyIndexedRing
    points::Vector{Tuple{Float64, Float64}}
    index::NaturalIndex{Extents.Extent{(:X, :Y), NTuple{2, NTuple{2, Float64}}}}
end

function NaturallyIndexedRing(points::Vector{Tuple{Float64, Float64}}; nodecapacity = 32)
    index = NaturalIndex(GO.edge_extents(GI.LinearRing(points)); nodecapacity)
    return NaturallyIndexedRing(points, index)
end
NaturallyIndexedRing(ring::NaturallyIndexedRing) = ring

function GI.convert(::Type{NaturallyIndexedRing}, ::GI.LinearRingTrait, geom)
    points = GO.tuples(geom).geom
    return NaturallyIndexedRing(points)
end

Base.show(io::IO, ::MIME"text/plain", ring::NaturallyIndexedRing) = Base.show(io, ring)
Base.show(io::IO, ring::NaturallyIndexedRing) = print(io, "NaturallyIndexedRing($(length(ring.points)) points) with index $(sprint(show, ring.index))")

GI.ncoord(::GI.LinearRingTrait, ring::NaturallyIndexedRing) = 2
GI.is3d(::GI.LinearRingTrait, ring::NaturallyIndexedRing) = false
GI.ismeasured(::GI.LinearRingTrait, ring::NaturallyIndexedRing) = false

GI.ngeom(::GI.LinearRingTrait, ring::NaturallyIndexedRing) = length(ring.points)
GI.getgeom(::GI.LinearRingTrait, ring::NaturallyIndexedRing) = ring.points
GI.getgeom(::GI.LinearRingTrait, ring::NaturallyIndexedRing, i::Int) = ring.points[i]

Extents.extent(ring::NaturallyIndexedRing) = ring.index.extent

GI.isgeometry(::Type{<: NaturallyIndexedRing}) = true
GI.geomtrait(::NaturallyIndexedRing) = GI.LinearRingTrait()

function prepare_naturally(geom)
    return GO.apply(GI.PolygonTrait(), geom) do poly
        return GI.Polygon([GI.convert(NaturallyIndexedRing, GI.LinearRingTrait(), ring) for ring in GI.getring(poly)])
    end
end

end # module NaturalIndexing