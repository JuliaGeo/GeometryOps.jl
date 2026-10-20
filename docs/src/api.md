```@meta
CurrentModule = GeometryOps
```

# Full GeometryOps API documentation

!!! warning
    This page is still very much WIP!

Documentation for [GeometryOps](https://github.com/JuliaGeo/GeometryOps.jl)'s full API (only for reference!).

```@index
```

## [`apply` and associated functions](@id Primitive-functions)
```@docs
apply
applyreduce
reproject
transform
```

## General geometry methods

### OGC methods
```@docs
GeometryOps.contains
coveredby
covers
crosses
disjoint
intersects
overlaps
touches
within
```

### Other general methods
```@docs
equals
centroid
distance
signed_distance
area
signed_area
angles
embed_extent
```

## Barycentric coordinates

```@docs
barycentric_coordinates
barycentric_coordinates!
barycentric_interpolate
```

## Other methods

```@autodocs
Modules = [GeometryOps]
# `reproject` is documented under "`apply` and associated functions" above.
Filter = t -> t !== GeometryOps.reproject
```

## Core types

```@autodocs
Modules = [GeometryOpsCore]
```

## Submodules

These are internal submodules of GeometryOps.  They are not part of the public API, so
anything here may change or disappear in a patch release, but they are documented for
reference (and because GeometryOps' own source code links to them).

### UnitSpherical

`GeometryOps.UnitSpherical` adds the orientation predicates, arc intersection and arc
extents on top of the UnitSphericalGeometryCore sublibrary, whose API it re-exports.

```@autodocs
Modules = [GeometryOps.UnitSpherical]
```

### UnitSphericalGeometryCore

The unit-spherical point and cap types, geographic transformations, `slerp` and bounding
caps live in the UnitSphericalGeometryCore sublibrary of this repository, which
GeometryOps depends on.

```@autodocs
Modules = [UnitSphericalGeometryCore]
```

### RobustCrossProduct

```@autodocs
Modules = [UnitSphericalGeometryCore.RobustCrossProduct]
```

### SpatialTreeInterface

The spatial tree interface and its traversals live in the SpatialTreeInterface
sublibrary of this repository, which GeometryOps depends on and re-exposes as
`GeometryOps.SpatialTreeInterface`.

```@autodocs
Modules = [SpatialTreeInterface]
```

### FlexibleRTrees

```@autodocs
Modules = [GeometryOps.FlexibleRTrees]
```

### NaturalIndexing

```@autodocs
Modules = [GeometryOps.NaturalIndexing]
```

### LoopStateMachine

`LoopStateMachine` is a submodule of SpatialTreeInterface, available as
`GeometryOps.LoopStateMachine`.

```@autodocs
Modules = [SpatialTreeInterface.LoopStateMachine]
```
