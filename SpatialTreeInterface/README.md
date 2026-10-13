# SpatialTreeInterface.jl


A simple interface for spatial tree types, and the traversals that run over them.

## What is a spatial tree?

- Nodes carry extents (bounding boxes); a parent's extent encompasses all of its children.
- Leaf nodes refer to the geometries they hold by index into the original collection.

## Why is this useful?

- Algorithms can be written against the interface and work with any spatial tree type, without knowing its details — for example, dual-tree traversal.
- Tree types can be swapped freely depending on the problem at hand.

This is a zero-cost interface if implemented correctly.  Verified implementations exist for "flat" trees like the "Natural Index" from `tg`, and "hierarchical" trees like the `STRtree` from `SortTileRecursiveTree.jl` (which this package supports out of the box, through a package extension that loads when SortTileRecursiveTree is loaded).

## Interface

- `isspatialtree(tree)::Bool`
- `isleaf(node)::Bool` - is the node a leaf node?  In this context, a leaf node is a node that does not have other nodes as its children, but stores a list of indices and extents (even if implicit).
- `getchild(node)` - get the children of a node.  This may be materialized if necessary or available, but can also be lazy (like a generator).
- `getchild(node, i)` - get the `i`-th child of a node.
- `nchild(node)::Int` - the number of children of a node.
- `child_indices_extents(node)` - an iterator over the indices and extents of the children of a **leaf** node.

These are the only methods that are required to be implemented.

Optionally, one may define:
- `node_extent(node)` - get the extent of a node.  This falls back to `GeoInterface.extent` but can potentially be overridden if you want to return a different but extent-like object.
- `node_extent_is_expensive(node)::Bool` - whether `node_extent` computes the extent rather than reading a stored one.  Defaults to `false`; set it on trees that derive extents on access, and `dual_depth_first_search` will cache child extents instead of re-deriving them.

They enable the generic query functions described below.

## Traversals

- `depth_first_search(f, predicate, tree)` - call `f(i)` for each index `i` in `tree` that satisfies `predicate(extent(i))`.
- `dual_depth_first_search(f, predicate, tree1, tree2)` - call `f(i1, i2)` for each index `i1` in `tree1` and `i2` in `tree2` that satisfies `predicate(extent(i1), extent(i2))`.

These are both completely non-allocating, and will only call `f` for indices that satisfy the predicate.  (Trees that opt into `node_extent_is_expensive` trade one small vector per visited internal node for far fewer `node_extent` calls.)

A standard query is built on top of `depth_first_search`:

```julia
using SpatialTreeInterface, SortTileRecursiveTree, Extents

tree = STRtree([Extents.Extent(X = (x, x + 1), Y = (y, y + 1)) for x in 1:10, y in 1:10])
query(tree, Extents.Extent(X = (1, 2), Y = (1, 2)))   # indices of extents that intersect
query(tree, (1.5, 1.5))                               # any GeoInterface geometry works too
query(tree, ext -> ext.X[1] > 5)                      # or any predicate on extents
```

`FlatNoTree(collection)` wraps any collection as a single leaf, for brute-force search without building a tree.

`spatialtree([manifold], geometries)` is an empty generic function that tree-providing packages extend, so that code written against the interface has one entry point for building a default tree.

## LoopStateMachine

The `SpatialTreeInterface.LoopStateMachine` submodule lets the callback `f` steer the traversal it is called from.  Return an `Action` from `f`:

- `Action(:continue)` / `Action(:break)` - the matching loop keyword, in the loop that called `f`.
- `Action(:return, x)` - make the function running the loop return `x`.
- `Action(:full_return, x)` - make the function running the loop return the `Action` itself, which propagates out of every enclosing `@controlflow` - i.e. out of a whole recursive traversal.  Use it to stop a tree query as soon as the first hit is found.

```julia
using SpatialTreeInterface.LoopStateMachine: Action

first_hit = SpatialTreeInterface.depth_first_search(Base.Fix1(Extents.intersects, ext), tree) do i
    return Action(:full_return, i)
end
first_hit.x   # the first matching index
```

The `@controlflow f(...)` macro is what processes these actions; use it in your own loops to give callers the same control.
