"""
    SpatialTreeInterface

A minimal interface for spatial tree types (R-trees, STR-trees, natural
indices, flat lists, ...), plus the two traversals most spatial algorithms
need: a single-tree [`depth_first_search`](@ref) and a
[`dual_depth_first_search`](@ref) that walks two trees in tandem.

A type implements the interface by defining [`isspatialtree`](@ref),
[`isleaf`](@ref), [`getchild`](@ref), [`nchild`](@ref) and
[`child_indices_extents`](@ref); see the README for the full contract.

The [`LoopStateMachine`](@ref) submodule provides `Action` and `@controlflow`,
which let a callback passed to a traversal `continue`, `break` or return out
of the loop it is called from.
"""
module SpatialTreeInterface

include("LoopStateMachine.jl")
using .LoopStateMachine: @controlflow, Action

import Extents
import GeoInterface as GI

export query
export FlatNoTree
export spatialtree

@static if VERSION >= v"1.11"
    eval(Meta.parse("""
    public isspatialtree, isleaf, getchild, nchild, child_indices_extents,
        node_extent, node_extent_is_expensive,
        depth_first_search, dual_depth_first_search, sanitize_predicate,
        LoopStateMachine
    """))
end

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
    spatialtree([manifold], geometries)

Build a spatial tree over an iterable of geometries.

SpatialTreeInterface defines no methods for this function; it exists so that
packages which provide trees (GeometryOps.jl, for example) can all extend one
shared entry point, and code written against the interface can build a tree
without naming a concrete tree type.
"""
function spatialtree end


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
function sanitize_predicate(::GI.AbstractTrait, pred)
    ext = GI.extent(pred)
    isnothing(ext) && throw(ArgumentError("no extent found on $(typeof(pred)); cannot query a tree with it"))
    return sanitize_predicate(ext)
end
sanitize_predicate(pred::Extents.Extent) = Base.Fix1(Extents.intersects, pred)


end # module SpatialTreeInterface
