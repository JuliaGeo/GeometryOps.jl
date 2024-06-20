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

TODO: should we switch to `Static.jl`?
=#
abstract type BoolsAsTypes end
struct _True <: BoolsAsTypes end
struct _False <: BoolsAsTypes end

@inline _booltype(x::Bool)::BoolsAsTypes = x ? _True() : _False()
@inline _booltype(x::BoolsAsTypes)::BoolsAsTypes = x


const TuplePoint{N, T} = NTuple{N, T} where {N, T}
const TupleEdge{N, T} = Tuple{TuplePoint{N, T}, TuplePoint{N, T}} where {N, T}

TuplePoint(geom, ::Type{T} = Float64) where T = _TuplePoint(T, GI.trait(geom), geom)
_TuplePoint(::Type{T}, ::GI.PointTrait, geom) where T = tuples(geom, T)
_TuplePoint(::Type{T}, trait::GI.AbstractTrait, _) where T = throw(ArgumentError("Geometry with trait $trait cannot be made into a point."))

TuplePoint_2D(vals, ::Type{T} = Float64) where T = TuplePoint{2, T}((GI.x(vals), GI.y(vals)))

TuplePoint_3D(vals, ::Type{T} = Float64, M = _False()) where T = _TuplePoint_3D(T, _booltype(M), vals)
_TuplePoint_3D(::Type{T}, ::_False, vals) where T = TuplePoint{3, T}((GI.x(vals), GI.y(vals), GI.z(vals)))
_TuplePoint_3D(::Type{T}, ::_True, vals) where T = TuplePoint{3, T}((GI.x(vals), GI.y(vals), GI.m(vals)))

TuplePoint_4D(vals, ::Type{T} = Float64) where T = TuplePoint{4, T}((GI.x(vals), GI.y(vals), GI.z(vals), GI.m(vals)))

#=
## `SVPoint`


=#
struct SVPoint{N, T, Z, M} <: GeometryBasics.StaticArraysCore.StaticVector{N,T}
    vals::NTuple{N,T}
end
Base.getindex(p::SVPoint, i::Int64) = p.vals[i]
# TODO: overload `similar_type``

const SVEdge{N, T, Z, M} = Tuple{SVPoint{N,T,Z,M}, SVPoint{N,T,Z,M}} where {N,T,Z,M}

# General SVPoint constructor when point type/size isn't known
SVPoint(geom, ::Type{T} = Float64) where T = _SVPoint(T, GI.trait(geom), geom)
_SVPoint(::Type{T}, ::GI.PointTrait, geom) where T = svpoints(geom, T)
_SVPoint(::Type{T}, trait::GI.AbstractTrait, _) where T = throw(ArgumentError("Geometry with trait $trait cannot be made into a point."))

# Syntactic sugar for type stability within functions with known point types
const PointType2D{T} = SVPoint{2, T, false, false} where T
const PointType3D{T} = SVPoint{3, T, true, false} where T
const PointType3DM{T} = SVPoint{3, T, false, true} where T
const PointType4D{T} = SVPoint{4, T, true, true} where T

SVPoint_2D(vals, ::Type{T} = Float64) where T <: AbstractFloat = PointType2D{T}(TuplePoint_2D(vals, T))

SVPoint_3D(vals, ::Type{T} = Float64, M = _False()) where {T <: AbstractFloat} = _SVPoint_3D(T, _booltype(M), vals)
_SVPoint_3D(T, M::_False, vals) = PointType3D{T}(TuplePoint_3D(vals, T, M))
_SVPoint_3D(T, M::_True, vals) = PointType3DM{T}(TuplePoint_3D(vals, T, M))

SVPoint_4D(vals, ::Type{T} = Float64) where T <: AbstractFloat = PointType4D{T}(TuplePoint_4D(vals, T))

#=
Get type of points and polygons made through library functionality (e.g. clipping)
TODO: Increase type options as library expands capabilities
=#
const PolyType2D{T} = GI.Polygon{false, false, Vector{GI.LinearRing{false, false, Vector{PointType2D{T}}, Nothing, Nothing}}, Nothing, Nothing} where T

const Edge{N, T} = Union{TupleEdge{N, T}, SVEdge{N, T, Z, M}} where {N, T, Z, M}
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
struct GEOS
    params::NamedTuple
end

function GEOS(; params...)
    nt = NamedTuple(params)
    return GEOS(nt)
end
# These are definitions for convenience, so we don't have to type out 
# `alg.params` every time.
Base.get(alg::GEOS, key, value) = Base.get(alg.params, key, value)
Base.get(f::Function, alg::GEOS, key) = Base.get(f, alg.params, key)

"""
    enforce(alg::GO.GEOS, kw::Symbol, f)

Enforce the presence of a keyword argument in a `GEOS` algorithm, and return `alg.params[kw]`.

Throws an error if the key is not present, and mentions `f` in the error message (since there isn't 
a good way to get the name of the function that called this method).
"""
function enforce(alg::GEOS, kw::Symbol, f)
    if haskey(alg.params, kw)
        return alg.params[kw]
    else
        error("$(f) requires a `$(kw)` keyword argument to the `GEOS` algorithm, which was not provided.")
    end
end
