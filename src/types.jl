#=
# Types

This file defines some fundamental types used in GeometryOps.

!!! warning
    Unlike in other Julia packages, only some types are defined in this file, not all. 
    This is because we define types in the files where they are used, to make it easier to understand the code.

=#
export TraitTarget, GEOS
#=
## `TraitTarget`

This struct holds a trait parameter or a union of trait parameters.
It's essentially a way to construct unions.
=#

"""
    TraitTarget{T}

This struct holds a trait parameter or a union of trait parameters.

It is primarily used for dispatch into methods which select trait levels, 
like `apply`, or as a parameter to `target`.

## Constructors
```julia
TraitTarget(GI.PointTrait())
TraitTarget(GI.LineStringTrait(), GI.LinearRingTrait()) # and other traits as you may like
TraitTarget(TraitTarget(...))
# There are also type based constructors available, but that's not advised.
TraitTarget(GI.PointTrait)
TraitTarget(Union{GI.LineStringTrait, GI.LinearRingTrait})
# etc.
```

"""
struct TraitTarget{T} end
TraitTarget(::Type{T}) where T = TraitTarget{T}()
TraitTarget(::T) where T<:GI.AbstractTrait = TraitTarget{T}()
TraitTarget(::TraitTarget{T}) where T = TraitTarget{T}()
TraitTarget(::Type{<:TraitTarget{T}}) where T = TraitTarget{T}()
TraitTarget(traits::GI.AbstractTrait...) = TraitTarget{Union{map(typeof, traits)...}}()


Base.in(::Trait, ::TraitTarget{Target}) where {Trait <: GI.AbstractTrait, Target} = Trait <: Target

#=
## `BoolsAsTypes`

In `apply` and `applyreduce`, we pass `threading` and `calc_extent` as types, not simple boolean values.  

This is to help compilation - with a type to hold on to, it's easier for 
the compiler to separate threaded and non-threaded code paths.

Note that if we didn't include the parent abstract type, this would have been really 
type unstable, since the compiler couldn't tell what would be returned!

We had to add the type annotation on the `_booltype(::Bool)` method for this reason as well.
=#
abstract type BoolsAsTypes end
struct _True <: BoolsAsTypes end
struct _False <: BoolsAsTypes end

@inline _booltype(x::Bool)::BoolsAsTypes = x ? _True() : _False()
@inline _booltype(x::BoolsAsTypes)::BoolsAsTypes = x

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
    GEOS(; params...)

A struct which instructs the method it's passed to as an algorithm
to use the appropriate GEOS function via `LibGEOS.jl` for the operation.

Dispatch is generally carried out using the names of the keyword arguments.
For example, `segmentize` will only accept a `GEOS` struct with only a 
`max_distance` keyword, and no other.

It's generally a lot slower than the native Julia implementations, since
it must convert to the LibGEOS implementation and back - so be warned!
"""
struct GEOS{Keys}
    params::NamedTuple{Keys}
    function GEOS{Keys}(arg) where Keys
        return new{Keys}(NamedTuple{Keys}(arg))
    end
end

function GEOS(; params...)
    nt = NamedTuple(params)
    return GEOS{keys(nt)}(nt)
end

