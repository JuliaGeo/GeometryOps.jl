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

```@autodocs
Modules = [GeometryOps.UnitSpherical]
```

### RobustCrossProduct

```@autodocs
Modules = [GeometryOps.UnitSpherical.RobustCrossProduct]
```

### SpatialTreeInterface

```@autodocs
Modules = [GeometryOps.SpatialTreeInterface]
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

```@autodocs
Modules = [GeometryOps.LoopStateMachine]
```
