# SpatialTreeInterface.jl

A simple interface for spatial tree types.

## What is a spatial tree?

- 2 dimensional extents
- Parent nodes encompass all leaf nodes
- Leaf nodes contain references to the geometries they represent as indices (or so we assume here)

## Why is this useful?

- It allows us to write algorithms that can work with any spatial tree type, without having to know the details of the tree type.
    - for example, dual tree traversal / queries
- It allows us to flexibly and easily swap out and use different tree types, depending on the problem at hand.

This is also a zero cost interface if implemented correctly!  Verified implementations exist for "flat" trees like the "Natural Index" from `tg`, and "hierarchical" trees like the `STRtree` from `SortTileRecursiveTree.jl`.

## Interface

- `isspatialtree(tree)::Bool`
- `isleaf(node)::Bool` - is the node a leaf node?  In this context, a leaf node is a node that does not have other nodes as its children, but stores a list of indices and extents (even if implicit).
- `getchild(node)` - get the children of a node.  This may be materialized if necessary or available, but can also be lazy (like a generator).
- `getchild(node, i)` - get the `i`-th child of a node.
- `nchild(node)::Int` - the number of children of a node.
- `child_indices_extents(node)` - an iterator over the indices and extents of the children of a **leaf** node.

These are the only methods that are required to be implemented.  

Optionally, one may define:
- `node_extent(node)` - get the extent of a node.  This falls back to `GI.extent` but can potentially be overridden if you want to return a different but extent-like object.

They enable the generic query functions described below:

## Query functions

- `do_query(f, predicate, node)` - call `f(i)` for each index `i` in `node` that satisfies `predicate(extent(i))`.
- `do_dual_query(f, predicate, tree1, tree2)` - call `f(i1, i2)` for each index `i1` in `tree1` and `i2` in `tree2` that satisfies `predicate(extent(i1), extent(i2))`.

These are both completely non-allocating, and will only call `f` for indices that satisfy the predicate.
You can of course build a standard query interface on top of `do_query` if you want - that's simply:
```julia
a = Int[]
do_query(Base.Fix1(push!, a), predicate, node)
```
where `predicate` might be `Base.Fix1(Extents.intersects, extent_to_query)`.

