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
