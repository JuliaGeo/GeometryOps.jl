module SpatialTreeInterface

import ..LoopStateMachine: @controlflow

import Extents
import GeoInterface as GI
import AbstractTrees

# public isspatialtree, isleaf, getchild, nchild, child_indices_extents, node_extent
export query
export FlatNoTree

# The spatial tree interface and its implementations are defined here.
include("interface.jl")
include("implementations.jl")

# Here we have some algorithms that use the spatial tree interface.
# The first file holds a single depth-first search, i.e., a single-tree query.
include("depth_first_search.jl")

# The second file holds a dual depth-first search, i.e., a dual-tree query.
# This iterates over two trees simultaneously, and is substantially more efficient
# than two separate single-tree queries since it can prune branches in tandem as it
# descends into the trees.
include("dual_depth_first_search.jl")


"""
    query(tree, predicate)

Return a sorted list of indices of the tree that satisfy the predicate.
"""
function query(tree, predicate)
    a = Int[]
    depth_first_search(Base.Fix1(push!, a), sanitize_predicate(predicate), tree)
    return sort!(a)
end


"""
    sanitize_predicate(pred)

Convert a predicate to a function that returns a Boolean.

If `pred` is an Extent, convert it to a function that returns a Boolean by intersecting with the extent.
If `pred` is a geometry, convert it to an extent first, then wrap in Extents.intersects.

Otherwise, return the predicate unchanged.


Users and developers may overload this function to provide custom behaviour when something is passed in.
"""
sanitize_predicate(pred) = sanitize_predicate(GI.trait(pred), pred)
sanitize_predicate(::Nothing, pred) = pred
sanitize_predicate(::GI.AbstractTrait, pred) = sanitize_predicate(GI.extent(pred))
sanitize_predicate(pred::Extents.Extent) = Base.Fix1(Extents.intersects, pred)


end # module SpatialTreeInterface