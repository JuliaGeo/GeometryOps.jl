#=
# Types

This defines core types that the GeometryOps ecosystem uses, 
and that are usable in more than just GeometryOps.

=#

#=
## `Manifold`

A manifold is mathematically defined as a topological space that resembles Euclidean space locally.

In GeometryOps (and geodesy more generally), there are three manifolds we care about:
- [`Linear`](@ref): the 2d plane, a completely Euclidean manifold
- [`Spherical`](@ref): the unit sphere, but one where areas are multiplied by the radius of the Earth.  This is not Euclidean globally, but all map projections attempt to represent the sphere on the Euclidean 2D plane to varying degrees of success.
- [`Geodesic`](@ref): the ellipsoid, the closest we can come to representing the Earth by a simple geometric shape.  Parametrized by `semimajor_axis` and `inv_flattening`.

Generally, we aim to have `Linear` and `Spherical` be operable everywhere, whereas `Geodesic` will only apply in specific circumstances.
Currently, those circumstances are `area` and `segmentize`, but this could be extended with time and https://github.com/JuliaGeo/SphericalGeodesics.jl.
=#

export Linear, Spherical, Geodesic
export TraitTarget
export BoolsAsTypes, _True, _False, _booltype

"""
    abstract type Manifold

A manifold is mathematically defined as a topological space that resembles Euclidean space locally.

We use the manifold definition to define the space in which an operation should be performed, or where a geometry lies.

Currently we have [`Linear`](@ref), [`Spherical`](@ref), and [`Geodesic`](@ref) manifolds.
"""
abstract type Manifold end

"""
    Linear()

A linear manifold means that the space is completely Euclidean,
and planar geometry suffices.
"""
struct Linear <: Manifold
end

"""
    Spherical(; radius)

A spherical manifold means that the geometry is on the 3-sphere (but is represented by 2-D longitude and latitude).  

## Extended help

!!! note
    The traditional definition of spherical coordinates in physics and mathematics, 
    ``r, \\theta, \\phi``, uses the _colatitude_, that measures angular displacement from the `z`-axis.  
    
    Here, we use the geographic definition of longitude and latitude, meaning
    that `lon` is longitude between -180 and 180, and `lat` is latitude between 
    `-90` (south pole) and `90` (north pole).
"""
Base.@kwdef struct Spherical{T} <: Manifold
    radius::T = 6371008.8
end

"""
    Geodesic(; semimajor_axis, inv_flattening)

A geodesic manifold means that the geometry is on a 3-dimensional ellipsoid, parameterized by `semimajor_axis` (``a`` in mathematical parlance)
and `inv_flattening` (``1/f``).

Usually, this is only relevant for area and segmentization calculations.  It becomes more relevant as one grows closer to the poles (or equator).
"""
Base.@kwdef struct Geodesic{T} <: Manifold
    semimajor_axis::T = 6378137,0
    inv_flattening::T = 298.257223563
end

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

"""
    abstract type BoolsAsTypes

"""
abstract type BoolsAsTypes end

"""
    struct _True <: BoolsAsTypes

A struct that means `true`.
"""
struct _True <: BoolsAsTypes end

"""
    struct _False <: BoolsAsTypes

A struct that means `false`.
"""
struct _False <: BoolsAsTypes end

"""
    _booltype(x)

Returns a [`BoolsAsTypes`](@ref) from `x`, whether it's a boolean or a BoolsAsTypes.
"""
function _booltype end

@inline _booltype(x::Bool)::BoolsAsTypes = x ? _True() : _False()
@inline _booltype(x::BoolsAsTypes)::BoolsAsTypes = x
