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
```

## Core types
```@autodocs
Modules = [GeometryOpsCore]
```
