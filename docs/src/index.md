```@meta
CurrentModule = GeometryOps
```

# GeometryOps.jl

GeometryOps.jl is a package for geometric calculations on (primarily 2D) geometries.

The driving idea behind this package is to unify all the disparate packages for geometric calculations in Julia, and make them [GeoInterface.jl](https://github.com/JuliaGeo/GeoInterface.jl)-compatible. We seem to be focusing primarily on 2/2.5D geometries for now.

Most of the usecases are driven by GIS and similar Earth data workflows, so this might be a bit specialized towards that, but methods should always be general to any coordinate space.

We welcome contributions, either as pull requests or discussion on issues!

## Main concepts

### The `apply` paradigm

!!! note
    See the [Primitive Functions](@ref) page for more information on this.