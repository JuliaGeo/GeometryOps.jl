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

getchild(node) = AbstractTrees.children(node)

"""
    getchild(node, i)

Return the `i`-th child of a node.
"""
getchild(node, i) = getchild(node)[i]

"""
    nchild(node)

Return the number of children of a node.
"""
nchild(node) = length(getchild(node))

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
    return zip(1:nchild(node), getchild(node))
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
