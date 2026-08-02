# # Types

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
    RTree(algorithm::BulkLoadAlgorithm, data; nodecapacity = 16, extents = nothing, indices = nothing)

A packed R-tree over the extents of `data` (anything `GI.extent` accepts —
geometries, or `Extents.Extent`s themselves), of any dimensionality, bulk
loaded in the order chosen by `algorithm`.

Pass a vector as `extents` (one per element of `data`, in order) to index
`data` by precomputed extents instead of `GI.extent` — for payload elements
that carry no extent of their own, or extents computed in another coordinate
space.  The tree takes ownership of the vector (`Unsorted` aliases it as the
leaf level rather than copying).

Pass a vector as `indices` (one per element of `data`, in order) to have
queries report those indices instead of positions in `data` — for indexing a
filtered subset of a larger collection, say one whose `missing` entries were
dropped, while still addressing the collection it came from.

The tree is flat and fully concrete: `levels[1]` is the coarsest level and
`levels[end]` holds the leaf extents in packed order, with `indices` mapping
each leaf slot back to its position in `data` (or, when `indices` was given,
to the caller's index for that element).  Queries through
SpatialTreeInterface therefore return indices into `data`, which the tree
keeps as `tree.data` so hits map straight back to elements wherever the
tree travels.
"""
struct RTree{A <: BulkLoadAlgorithm, E <: Extents.Extent, D <: AbstractVector, I <: AbstractVector{Int}}
    algorithm::A
    nodecapacity::Int
    extent::E
    levels::Vector{Vector{E}}   # levels[1] = coarsest, levels[end] = leaf extents (packed order)
    indices::I                  # leaf slot -> index into `data` (`Base.OneTo` when unpermuted)
    data::D                     # the indexed collection
end

function RTree(algorithm::A, data; nodecapacity::Int = 16,
        extents::Union{Nothing, Vector{<:Extents.Extent}} = nothing,
        indices::Union{Nothing, AbstractVector{Int}} = nothing) where A <: BulkLoadAlgorithm
    nodecapacity >= 2 || throw(ArgumentError("`nodecapacity` must be at least 2, got $nodecapacity"))
    items = data isa AbstractVector ? data : collect(data)
    isempty(items) && throw(ArgumentError("cannot build an `RTree` from an empty collection"))
    exts = if extents === nothing
        E = typeof(GI.extent(first(items)))
        E[GI.extent(x) for x in items]
    else
        length(extents) == length(items) || throw(ArgumentError(
            "`extents` must have one entry per element of `data`, got $(length(extents)) for $(length(items))"))
        extents
    end
    isnothing(indices) || length(indices) == length(items) || throw(ArgumentError(
        "`indices` must have one entry per element of `data`, got $(length(indices)) for $(length(items))"))
    perm = loadorder(algorithm, exts, nodecapacity)
    leaves = perm isa Base.OneTo ? exts : exts[perm]
    levels = _pack_levels(leaves, nodecapacity)
    total = reduce(Extents.union, levels[1])
    leafindices = isnothing(indices) ? perm : indices[perm]
    return RTree(algorithm, nodecapacity, total, levels, leafindices, items)
end

"""
    RTree(m::Manifold, algorithm::BulkLoadAlgorithm, data; nodecapacity = 16, indices = nothing)

Build the tree over each element's extent *on the manifold `m`*, via
`Extents.extent(m, x)`.  On `Spherical()` the leaves are the 3D Cartesian
boxes of the elements as regions on the unit sphere, covering the arc bulge
and enclosed poles that vertex extents miss.
"""
function RTree(m::Manifold, algorithm::BulkLoadAlgorithm, data; nodecapacity::Int = 16,
        indices::Union{Nothing, AbstractVector{Int}} = nothing)
    items = data isa AbstractVector ? data : collect(data)
    isempty(items) && throw(ArgumentError("cannot build an `RTree` from an empty collection"))
    return RTree(algorithm, items; nodecapacity, indices,
        extents = [Extents.extent(m, x) for x in items])
end

Extents.extent(tree::RTree) = tree.extent

function Base.show(io::IO, tree::RTree{A}) where A
    print(io, "RTree{", nameof(A), "}(", length(tree.indices), " leaves, ",
        length(tree.levels), " levels, capacity ", tree.nodecapacity, ")")
end
Base.show(io::IO, ::MIME"text/plain", tree::RTree) = Base.show(io, tree)
