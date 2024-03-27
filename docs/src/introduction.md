# Introduction


GeometryOps.jl is a package for geometric calculations on (primarily 2D) geometries.

The driving idea behind this package is to unify all the disparate packages for geometric calculations in Julia, and make them [GeoInterface.jl](https://github.com/JuliaGeo/GeoInterface.jl)-compatible. We seem to be focusing primarily on 2/2.5D geometries for now.

Most of the usecases are driven by GIS and similar Earth data workflows, so this might be a bit specialized towards that, but methods should always be general to any coordinate space.

We welcome contributions, either as pull requests or discussion on issues!

## Main concepts

### The `apply` paradigm

!!! note
    See the [Primitive Functions](@ref Primitive-functions) page for more information on this.

The `apply` function allows you to decompose a given collection of geometries down to a certain level, and then operate on it. 

Functionally, it's similar to `map` in the way you apply it to geometries.

### What's this `GeoInterface.Wrapper` thing?

Write a comment about GeoInterface.Wrapper and why it helps in type stability to guarantee a particular return type.

