# # Interface
# Interface definition for spatial tree types.
# There is no abstract supertype here since it's impossible to enforce,
# but we do have a few methods that are common to all spatial tree types.

"""
    isspatialtree(tree)::Bool

Return true if the object is a spatial tree, false otherwise.

## Implementation notes

For type stability, if your spatial tree type is `MyTree`, you should define
`isspatialtree(::Type{MyTree}) = true`, and `isspatialtree(::MyTree)` will forward
to that method automatically.
"""
isspatialtree(::T) where T = isspatialtree(T)
isspatialtree(::Type{<: Any}) = false


"""
    getchild(node)
    getchild(node, i)

Accessor function to get the children of a node.

If invoked as `getchild(node)`, return an iterator over all the children of a node.  
This may be lazy, like a `Base.Generator`, or it may be materialized.

If invoked as `getchild(node, i)`, return the `i`-th child of a node.
"""
function getchild end 

getchild(node) = (getchild(node, i) for i in 1:nchild(node))

"""
    getchild(node, i)

Return the `i`-th child of a node.
"""
getchild(node, i) = error("getchild(node, i) is not implemented for node type $(typeof(node)), it must be implemented!")

"""
    nchild(node)

Return the number of children of a node.
"""
nchild(node) = error("nchild is not implemented for node type $(typeof(node)), it must be implemented!")

"""
    isleaf(node)

Return true if the node is a leaf node, i.e., there are no "children" below it.
[`getchild`](@ref) should still work on leaf nodes, though, returning an iterator over the extents stored in the node - and similarly for `getnodes.`
"""
isleaf(node) = error("isleaf is not implemented for node type $(typeof(node))")

"""
    child_indices_extents(node)

Return an iterator over the indices and extents of the children of a node.

Each value of the iterator should take the form `(i, extent)`.

This can only be invoked on leaf nodes!
"""
function child_indices_extents(node)
    children = getchild(node)
    if applicable(Base.keys, typeof(children)) 
        return ((i, node_extent(obj)) for (i, obj) in pairs(children))
    else
        return ((i, node_extent(obj)) for (i, obj) in enumerate(children))
    end
end

"""
    node_extent(node)

Return the extent like object of the node.  
Falls back to `GI.extent` by default, which falls back
to `Extents.extent`.  

Generally, defining `Extents.extent(node)` is sufficient here, and you
won't need to define this

The reason we don't use that directly is to give users of this interface
a way to define bounding boxes that are not extents, like spherical caps
and other such things.
"""
node_extent(node) = GI.extent(node)

"""
    node_extent_is_expensive(node)::Bool

Return true if [`node_extent`](@ref) computes the node's extent instead of
reading one the node already stores.  Defaults to false.

When true, [`dual_depth_first_search`](@ref) caches a node's child extents rather
than re-deriving them once per opposing child, at the cost of a small vector per
visited node.

## Implementation notes

Define this on the type - `node_extent_is_expensive(::Type{MyNode}) = true` - so
that it is known at compile time; `node_extent_is_expensive(::MyNode)` forwards
there automatically.
"""
node_extent_is_expensive(::T) where {T} = node_extent_is_expensive(T)
node_extent_is_expensive(::Type{<:Any}) = false
