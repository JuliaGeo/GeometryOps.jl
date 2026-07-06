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

include("types.jl")         # `BulkLoadAlgorithm`s and the `RTree` type
include("bulk_loading.jl")  # `loadorder` methods and bottom-up packing
include("hilbert.jl")       # Hilbert keys for `HPR`'s sort
include("interface.jl")     # SpatialTreeInterface implementation and `query`

end # module FlexibleRTrees
