# GeometryOpsCore

This is a "core" package for [GeometryOps.jl](https://github.com/JuliaGeo/GeometryOps.jl), that defines some basic primitive functions and types for GeometryOps.

Generally, you would depend on this to use either the GeometryOps types (like `Linear`, `Spherical`, etc) or the primitive functions like `apply`, `applyreduce`, `flatten`, etc. 
All of these are also accessible from GeometryOps, so it's preferable that you use GeometryOps directly.

Tests are in the main GeometryOps tests, we don't have separate tests for GeometryOpsCore since it's in a monorepo structure.