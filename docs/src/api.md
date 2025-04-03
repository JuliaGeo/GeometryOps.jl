```@meta
CurrentModule = GeometryOps
```

# Full GeometryOps API documentation

!!! warning
    This page is still very much WIP!

Documentation for [GeometryOps](https://github.com/JuliaGeo/GeometryOps.jl)'s full API (only for reference!).

## `apply` and associated functions
```@docs
apply
applyreduce
reproject
transform
```

## Manifolds

```@docs
Manifold
Planar
Spherical
Geodesic
AutoManifold
```

## Algorithms

```@docs
Algorithm
ManifoldIndependentAlgorithm
SingleManifoldAlgorithm
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

## All other methods


```@index
```

```@autodocs
Modules = [GeometryOps, GeometryOpsCore]
```
