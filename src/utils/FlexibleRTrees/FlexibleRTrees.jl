# # FlexibleRTrees

#=
A packed (bulk-loaded, static) R-tree over `Extents.Extent`s of any
dimensionality, with a pluggable bulk-load algorithm — sort-tile-recursive
([`STR`](@ref)), Hilbert-packed ([`HPR`](@ref)), or none
([`Unsorted`](@ref)) — behind one tree type.

Storage is flat: `RTree{A, E}` holds a vector of per-level extent vectors
plus a leaf permutation, and is a concrete type at any size or depth.  A
bulk-load algorithm chooses only the *leaf order*, via [`loadorder`](@ref);
packing always unions consecutive runs of `nodecapacity` extents, bottom-up.
Upper levels therefore group runs of the leaf order rather than re-tiling
each level: Hilbert order is spatially local at every scale so `HPR` packs
tightly, while `STR`'s upper levels are slightly looser than a re-tiled
pointer tree's.

The tree implements SpatialTreeInterface, so `depth_first_search` /
`dual_depth_first_search` (and everything built on them) work unchanged.
Leaf queries yield indices into the *original* input collection.

Parts of the construction logic are adapted from
[SortTileRecursiveTree.jl](https://github.com/maxfreu/SortTileRecursiveTree.jl) (MIT).

```julia
tree = RTree(HPR(), extents)                    # or STR(), Unsorted()
hits = query(tree, Extents.Extent(X = (0, 1), Y = (0, 1)))
```
=#

module FlexibleRTrees

import GeoInterface as GI
import Extents
using StaticArrays: MVector

using ..SpatialTreeInterface
import ..SpatialTreeInterface: isspatialtree, isleaf, nchild, getchild,
    child_indices_extents, node_extent, depth_first_search

export RTree, BulkLoadAlgorithm, STR, HPR, Unsorted, query

# ## Bulk-load algorithms

"""
    BulkLoadAlgorithm

Supertype for the algorithms that decide the *leaf order* of an [`RTree`](@ref).
Packing is always "union consecutive runs of `nodecapacity`, bottom-up"; the
algorithm only chooses the order, via a [`loadorder`](@ref) method.
"""
abstract type BulkLoadAlgorithm end

"""
    STR()

Sort-tile-recursive ordering (Leutenegger et al., 1997), generalized to any
dimensionality: sort by center along the first dimension, cut into slabs,
recurse within each slab on the remaining dimensions.
"""
struct STR <: BulkLoadAlgorithm end

"""
    HPR()

Hilbert-packed ordering, as in JTS's `HPRtree`: sort by the Hilbert-curve
index of each extent's center.  Hilbert order is spatially local at every
scale, which suits this tree's consecutive-run packing particularly well.
"""
struct HPR <: BulkLoadAlgorithm end

"""
    Unsorted()

Keep the input order (no sort).  Equivalent to natural indexing — good when
the input is already spatially coherent (e.g. the edges of a ring), and the
baseline the sorting algorithms have to beat.
"""
struct Unsorted <: BulkLoadAlgorithm end

# ## The tree

"""
    RTree(algorithm::BulkLoadAlgorithm, data; nodecapacity = 16)

A packed R-tree over the extents of `data` (anything `GI.extent` accepts —
geometries, or `Extents.Extent`s themselves), of any dimensionality, bulk
loaded in the order chosen by `algorithm`.

The tree is flat and fully concrete: `levels[1]` is the coarsest level and
`levels[end]` holds the leaf extents in packed order, with `indices` mapping
each leaf slot back to its position in `data`.  Queries through
SpatialTreeInterface therefore return indices into the original collection.
"""
struct RTree{A <: BulkLoadAlgorithm, E <: Extents.Extent}
    algorithm::A
    nodecapacity::Int
    extent::E
    levels::Vector{Vector{E}}   # levels[1] = coarsest, levels[end] = leaf extents (packed order)
    indices::Vector{Int}        # leaf slot -> index into the original collection
end

function RTree(algorithm::A, data; nodecapacity::Int = 16) where A <: BulkLoadAlgorithm
    nodecapacity >= 2 || throw(ArgumentError("`nodecapacity` must be at least 2, got $nodecapacity"))
    isnothing(iterate(data)) && throw(ArgumentError("cannot build an `RTree` from an empty collection"))
    E = typeof(GI.extent(first(data)))
    extents = E[GI.extent(x) for x in data]
    perm = loadorder(algorithm, extents, nodecapacity)
    leaves = extents[perm]
    levels = _pack_levels(leaves, nodecapacity)
    total = reduce(Extents.union, levels[1])
    return RTree{A, E}(algorithm, nodecapacity, total, levels, perm)
end

Extents.extent(tree::RTree) = tree.extent

function Base.show(io::IO, tree::RTree{A}) where A
    print(io, "RTree{", nameof(A), "}(", length(tree.indices), " leaves, ",
        length(tree.levels), " levels, capacity ", tree.nodecapacity, ")")
end
Base.show(io::IO, ::MIME"text/plain", tree::RTree) = Base.show(io, tree)

# ## Leaf ordering

"""
    loadorder(algorithm, extents::Vector{<:Extents.Extent}, nodecapacity)::Vector{Int}

The permutation in which `algorithm` packs `extents` into leaves.  Implement
this for a new `BulkLoadAlgorithm` subtype to plug in another ordering.
"""
loadorder(::Unsorted, extents, nodecapacity) = collect(1:length(extents))
loadorder(::HPR, extents, nodecapacity) = sortperm(_hilbert_keys(extents))
function loadorder(::STR, extents::Vector{E}, nodecapacity) where E
    centers = [_center(e) for e in extents]
    perm = collect(1:length(extents))
    _str_tile!(perm, centers, 1, length(extents), 1, _ndims(E), nodecapacity)
    return perm
end

#=
One recursion level of N-dimensional sort-tile-recursive: sort the range by
the current dimension's center, cut it into `S ≈ P^(1/remaining)` slabs of
whole leaf pages, and tile each slab along the remaining dimensions.  After
the last dimension the consecutive `nodecapacity`-runs are the leaf tiles.
=#
function _str_tile!(perm, centers, lo, hi, dim, ndims, nodecapacity)
    len = hi - lo + 1
    len <= nodecapacity && return # a single leaf: internal order doesn't matter
    sort!(view(perm, lo:hi); by = i -> @inbounds(centers[i][dim]))
    dim == ndims && return
    P = cld(len, nodecapacity)                        # leaf pages in this range
    S = ceil(Int, P^(1 / (ndims - dim + 1)))          # slabs along this dimension
    slab = cld(P, S) * nodecapacity                   # items per slab (whole pages)
    i = lo
    while i <= hi
        _str_tile!(perm, centers, i, min(i + slab - 1, hi), dim + 1, ndims, nodecapacity)
        i += slab
    end
    return
end

# ## Packing

# Union consecutive runs of `nodecapacity` extents, bottom-up, until a level
# fits in one (implicit) root node.  Returns levels coarsest-first.
function _pack_levels(leaves::Vector{E}, nodecapacity::Int) where E
    levels = [leaves]
    current = leaves
    while length(current) > nodecapacity
        nparents = cld(length(current), nodecapacity)
        parents = Vector{E}(undef, nparents)
        for p in 1:nparents
            lo = (p - 1) * nodecapacity + 1
            hi = min(p * nodecapacity, length(current))
            acc = current[lo]
            for j in (lo + 1):hi
                acc = Extents.union(acc, @inbounds current[j])
            end
            parents[p] = acc
        end
        push!(levels, parents)
        current = parents
    end
    return reverse!(levels)
end

# ## Extent helpers

_ndims(::Type{Extents.Extent{K, V}}) where {K, V} = length(K)
_center(ext::Extents.Extent) = map(b -> (b[1] + b[2]) / 2, values(ext))

include("hilbert.jl")

# ## SpatialTreeInterface

"""
    RTreeNode{A, E}

A cursor into one node of an [`RTree`](@ref): the tree, the node's level
(0-based; the children of a level-`l` node live in `levels[l + 1]`), its
position within that level, and its extent.  All SpatialTreeInterface
methods traverse the tree through these cursors.  The children of one node
occupy one contiguous run of the next level's extent vector, which each
per-child method resolves once per node and then indexes into.  At the leaf
level, `child_indices_extents` maps leaf slots through `tree.indices`, so
queries return indices into the original collection despite the packed
reordering.
"""
struct RTreeNode{A <: BulkLoadAlgorithm, E <: Extents.Extent}
    tree::RTree{A, E}
    level::Int      # 0-based; children of a level-l node live in levels[l + 1]
    index::Int      # position within its level
    extent::E
end

Extents.extent(node::RTreeNode) = node.extent

isspatialtree(::Type{<:RTree}) = true
isspatialtree(::Type{<:RTreeNode}) = true

@inline _child_extents(node::RTreeNode) = node.tree.levels[node.level + 1]

@inline function _child_range(node::RTreeNode, child_extents)
    start_idx = (node.index - 1) * node.tree.nodecapacity + 1
    stop_idx = min(start_idx + node.tree.nodecapacity - 1, length(child_extents))
    return start_idx:stop_idx
end

isleaf(node::RTreeNode) = node.level == length(node.tree.levels) - 1

nchild(node::RTreeNode) = length(_child_range(node, _child_extents(node)))

function getchild(node::RTreeNode, i::Int)
    child_index = (node.index - 1) * node.tree.nodecapacity + i
    return RTreeNode(node.tree, node.level + 1, child_index, _child_extents(node)[child_index])
end

function getchild(node::RTreeNode)
    extents = _child_extents(node)
    tree, childlevel = node.tree, node.level + 1
    range = _child_range(node, extents)
    return (RTreeNode(tree, childlevel, ci, @inbounds extents[ci]) for ci in range)
end

function child_indices_extents(node::RTreeNode)
    extents = _child_extents(node)
    indices = node.tree.indices
    range = _child_range(node, extents)
    return ((@inbounds(indices[i]), @inbounds(extents[i])) for i in range)
end

# The tree itself acts as the (implicit) root node.
_rootnode(tree::RTree) = RTreeNode(tree, 0, 1, tree.extent)

isleaf(tree::RTree) = length(tree.levels) == 1
nchild(tree::RTree) = length(tree.levels[1])
getchild(tree::RTree) = getchild(_rootnode(tree))
getchild(tree::RTree, i) = getchild(_rootnode(tree), i)
child_indices_extents(tree::RTree) = child_indices_extents(_rootnode(tree))

# ## Queries

"""
    query(tree::RTree, extent_or_geom)

Indices (into the collection the tree was built from) of every leaf whose
extent intersects the given extent — or the extent of the given geometry —
in ascending order.
"""
query(tree::RTree, ext::Extents.Extent) =
    sort!(depth_first_search(Base.Fix1(Extents.intersects, ext), tree))
function query(tree::RTree, geom)
    ext = GI.extent(geom)
    isnothing(ext) && throw(ArgumentError("no extent found on $(typeof(geom))"))
    return query(tree, ext)
end

end # module FlexibleRTrees
