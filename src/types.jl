#=
# Types

This file defines some fundamental types used in GeometryOps.

!!! warning
    Unlike in other Julia packages, only some types are defined in this file, not all. 
    This is because we define types in the files where they are used, to make it easier to understand the code.

=#
export GEOS
#=

## `GEOS`

`GEOS` is a struct which instructs the method it's passed to as an algorithm
to use the appropriate GEOS function via `LibGEOS.jl` for the operation.

It's generally a lot slower than the native Julia implementations, but it's
useful for two reasons:
1. Functionality which doesn't exist in GeometryOps can be accessed through the GeometryOps API, but use GEOS in the backend until someone implements a native Julia version.
2. It's a good way to test the correctness of the native implementations.

=#
"""
    abstract type CLibraryPlanarAlgorithm <: GeometryOpsCore.SingleManifoldAlgorithm{Planar} end

This is a type which extends `GeometryOpsCore.SingleManifoldAlgorithm{Planar}`,
and is used as an abstract supertype for some C library based algorithms.

The type requires that algorithm structs be arranged as:
```
struct MyAlgorithm <: CLibraryPlanarAlgorithm
    manifold::Planar
    params::NamedTuple
end
```

Then you get a nice constructor for free, as well as the 
`get(alg, key, value)` and `get(alg, key) do ...` syntax.
Plus the [`enforce`](@ref) method, which will check that given keyword arguments
are present.
"""
abstract type CLibraryPlanarAlgorithm <: GeometryOpsCore.SingleManifoldAlgorithm{Planar} end


function (::Type{T})(; params...) where {T <: CLibraryPlanarAlgorithm}
    nt = NamedTuple(params)
    return T(Planar(), nt)
end
(T::Type{<: CLibraryPlanarAlgorithm})(params::NamedTuple) = T(Planar(), params)


# These are definitions for convenience, so we don't have to type out 
# `alg.params` every time.

Base.get(alg::CLibraryPlanarAlgorithm, key, value) = Base.get(alg.params, key, value)
Base.get(f::Function, alg::CLibraryPlanarAlgorithm, key) = Base.get(f, alg.params, key)

"""
    enforce(alg::CLibraryPlanarAlgorithm, kw::Symbol, f)

Enforce the presence of a keyword argument in a `GEOS` algorithm, and return `alg.params[kw]`.

Throws an error if the key is not present, and mentions `f` in the error message (since there isn't 
a good way to get the name of the function that called this method).

This applies to all `CLibraryPlanarAlgorithm` types, like [`GEOS`](@ref) and [`TG`](@ref).
"""
function enforce(alg::CLibraryPlanarAlgorithm, kw::Symbol, f)
    if haskey(alg.params, kw)
        return alg.params[kw]
    else
        error("$(f) requires a `$(kw)` keyword argument to the `GEOS` algorithm, which was not provided.")
    end
end



"""
    GEOS(; params...)

A struct which instructs the method it's passed to as an algorithm
to use the appropriate GEOS function via `LibGEOS.jl` for the operation.

Dispatch is generally carried out using the names of the keyword arguments.
For example, `segmentize` will only accept a `GEOS` struct with only a 
`max_distance` keyword, and no other.

It's generally somewhat slower than the native Julia implementations, since
it must convert to the LibGEOS implementation and back - so be warned!

## Extended help

This uses the [LibGEOS.jl](https://github.com/JuliaGeometry/LibGEOS.jl) package,
which is a Julia wrapper around the C library GEOS (https://trac.osgeo.org/geos).
"""
struct GEOS <: CLibraryPlanarAlgorithm # SingleManifoldAlgorithm{Planar}
    manifold::Planar
    params::NamedTuple
end

"""
    TG(; params...)

A struct which instructs the method it's passed to as an algorithm
to use the appropriate TG function via `TGGeometry.jl` for the operation.

It's generally a lot faster than the native Julia implementations, but only
supports planar manifolds / operations.  Also, it only supports geometric predicates,
specifically the ones which the underlying `tg` library supports.  These are:

[`equals`](@ref), [`intersects`](@ref), [`disjoint`](@ref), [`contains`](@ref), 
[`within`](@ref), [`covers`](@ref), [`coveredby`](@ref), and [`touches`](@ref).

## Extended help

This uses the [TGGeometry.jl](https://github.com/JuliaGeo/TGGeometry.jl) package,
which is a Julia wrapper around the `tg` C library (https://github.com/tidwall/tg).
"""
struct TG <: CLibraryPlanarAlgorithm
    manifold::Planar
    params::NamedTuple
end


