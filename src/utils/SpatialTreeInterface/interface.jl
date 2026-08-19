# # Interface
# Interface definition for spatial tree types.
# There is no abstract supertype here since it's impossible to enforce,
# but we do have a few methods that are common to all spatial tree types.
#
# ```@meta
# CollapsedDocStrings = true
# ```
#
# The methods that make up the interface are:
#
# ```@docs; canonical=false
# isspatialtree
# isleaf
# getchild
# nchild
# child_indices_extents
# node_extent
# ```

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
than re-deriving them once per opposing child.  The cache is one scratch stack
per traversal, so opting in does not allocate per visited node.

A node's children must all share an extent type, since they share one stack; see
[`children_extent_type`](@ref) for trees where that type changes with depth.

## Implementation notes

Define this on the type - `node_extent_is_expensive(::Type{MyNode}) = true` - so
that it is known at compile time; `node_extent_is_expensive(::MyNode)` forwards
there automatically.
"""
node_extent_is_expensive(::T) where {T} = node_extent_is_expensive(T)
node_extent_is_expensive(::Type{<:Any}) = false

"""
    children_extent_type(node)::Union{Type, Nothing}

Return the type [`node_extent`](@ref) gives for `node`'s *children*, or `nothing`
- the default - to say it is the same type as `node`'s own extent.

Only [`dual_depth_first_search`](@ref) on a tree that opts into
[`node_extent_is_expensive`](@ref) consults this: it caches a node's child
extents in a typed scratch stack, and needs the element type before it touches
a child.  Siblings must agree, but a tree whose extent type changes with depth
can say so here and still be traversed.

## Implementation notes

Define it on the node - `children_extent_type(node::MyNode) = MyChildExtent` - or
on its type if that reads better, `children_extent_type(::Type{MyNode}) = MyChildExtent`;
the node method falls back to the type one.  Either way the answer is a constant
per node type, so it folds away.
"""
children_extent_type(node) = children_extent_type(typeof(node))
children_extent_type(::Type{<:Any}) = nothing
