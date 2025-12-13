#=
# `Algorithm`s

An `Algorithm` is a type that describes the algorithm used to perform some [`Operation`](@ref).

An algorithm may be associated with one or many [`Manifold`](@ref)s.  It may either have the manifold as a field, or have it as a static parameter (e.g. `struct GEOS <: Manifold{Planar}`).


Algorithms are:
* Ways to perform an operation
* For example: LHuilier, Bessel, Ericsson for spherical area
* May be manifold agnostic (like simplification) or restrictive (like GEOS only works on planar, PROJ algorithm for arclength and area only works on geodesic)
* May or may not carry manifolds around, but manifold should always be accessible from manifold(alg) - it's not necessary that fixed manifold args can skip carrying the manifold around, eg in the case of Proj{Geodesic}.


## Single manifold vs manifold independent algorithms

Some algorithms only work on a single manifold (shoelace area only works on `Planar`)
and others work on any manifold (`simplify`, `segmentize`).  They are allowed to dispatch
on the manifold type but store that manifold as a field, and the fundamental algorithm is the 
same.  For example the segmentize algorithm is the same (distance along line) but the implementation
varies slightly depending on the manifold (planar, spherical, geodesic, etc).

Here's a simple example of two algorithms, one only on planar and one manifold independent:

```julia

struct MyExternalArbitraryPackageAlgorithm <: SingleManifoldAlgorithm{Planar}
    kw1::Int
    kw2::String
end # this already has the methods specified

struct MyIndependentAlgorithm{M <: Manifold} <: ManifoldIndependentAlgorithm{M}
    m::M
    kw1::Int
    kw2::String
end

MyIndependentAlgorithm(m::Manifold; kw1 = 1, kw2 = "hello") = MyIndependentAlgorithm(m, kw1, kw2)
```
=#

export Algorithm, AutoAlgorithm, ManifoldIndependentAlgorithm, SingleManifoldAlgorithm, NoAlgorithm

"""
    abstract type Algorithm{M <: Manifold}

The abstract supertype for all GeometryOps algorithms.  
These define how to perform a particular [`Operation`](@ref).

An algorithm may be associated with one or many [`Manifold`](@ref)s.  
It may either have the manifold as a field, or have it as a static parameter 
(e.g. `struct GEOS <: Algorithm{Planar}`).

## Interface

All `Algorithm`s must implement the following methods:

- `rebuild(alg, manifold::Manifold)` Rebuild algorithm `alg` with a new manifold 
  as passed in the second argument.  This may error and throw a [`WrongManifoldException`](@ref)
  if the manifold is not compatible with that algorithm.
- `manifold(alg::Algorithm)` Return the manifold associated with the algorithm.
- `best_manifold(alg::Algorithm, input)`: Return the best manifold for that algorithm, in the absence of
  any other context.  WARNING: this may change in future and is not stable!

The actual implementation is left to the implementation of that particular [`Operation`](@ref).

## Notable subtypes

- [`AutoAlgorithm`](@ref): Tells the [`Operation`](@ref) receiving 
  it to automatically select the best algorithm for its input data.
- [`ManifoldIndependentAlgorithm`](@ref): An abstract supertype for an algorithm that works on any manifold.
  The manifold must be stored in the algorithm for a `ManifoldIndependentAlgorithm`, and accessed via `manifold(alg)`.
- [`SingleManifoldAlgorithm`](@ref): An abstract supertype for an algorithm that only works on a 
  single manifold, specified in its type parameter.  `SingleManifoldAlgorithm{Planar}` is a special case
  that does not have to store its manifold, since that doesn't contain any information.  All other 
  `SingleManifoldAlgorithm`s must store their manifold, since they do contain information.
- [`NoAlgorithm`](@ref): A type that indicates no algorithm is to be used, essentially the equivalent
  of `nothing`.
"""
abstract type Algorithm{M <: Manifold} end

"""
    manifold(alg::Algorithm)::Manifold

Return the manifold associated with the algorithm.  

May be any subtype of [`Manifold`](@ref).
"""
function manifold end

# The below definition is a special case, since [`Planar`](@ref) has no contents, being a 
# singleton struct.
# If that changes in the future, then this method must be deleted.
manifold(::Algorithm{<: Planar}) = Planar()

"""
    best_manifold(alg::Algorithm, input)::Manifold
    
Return the best [`Manifold`](@ref) for the algorithm `alg` based on the given `input`.

May be any subtype of [`Manifold`](@ref).
"""
function best_manifold end

# ## Implementation of basic algorithm types

# ### `AutoAlgorithm`

"""
    AutoAlgorithm{T, M <: Manifold}(manifold::M, x::T)

Indicates that the [`Operation`](@ref) should automatically select the best algorithm for
its input data, based on the passed in manifold (may be an [`AutoManifold`](@ref)) and data 
`x`.

The actual implementation is left to the implementation of that particular [`Operation`](@ref).
"""
struct AutoAlgorithm{T, M <: Manifold} <: Algorithm{M} 
    manifold::M
    x::T
end

AutoAlgorithm(m::Manifold; kwargs...) = AutoAlgorithm(m, kwargs)
AutoAlgorithm(; kwargs...) = AutoAlgorithm(AutoManifold(), kwargs)

manifold(a::AutoAlgorithm) = a.manifold
rebuild(a::AutoAlgorithm, m::Manifold) = AutoAlgorithm(m, a.x)


# ### `ManifoldIndependentAlgorithm`

"""
    abstract type ManifoldIndependentAlgorithm{M <: Manifold} <: Algorithm{M}

The abstract supertype for a manifold-independent algorithm, i.e., one which may work on any manifold.

The manifold is stored in the algorithm for a `ManifoldIndependentAlgorithm`, and accessed via `manifold(alg)`.
"""
abstract type ManifoldIndependentAlgorithm{M <: Manifold} <: Algorithm{M} end


# ### `SingleManifoldAlgorithm`

"""
    abstract type SingleManifoldAlgorithm{M <: Manifold} <: Algorithm{M}

The abstract supertype for a single-manifold algorithm, i.e., one which is known to only work 
on a single manifold.

The manifold may be accessed via `manifold(alg)`.
"""
abstract type SingleManifoldAlgorithm{M <: Manifold} <: Algorithm{M} end

function (Alg::Type{<: SingleManifoldAlgorithm{M}})(m::M; kwargs...) where {M}
    # successful - the algorithm is designed for this manifold
    # in this case, just return `Alg(; kwargs...)`
    return Alg(; kwargs...)
end

function (Alg::Type{<: SingleManifoldAlgorithm{M}})(m::Manifold; kwargs...) where {M}
    # this catches the case where the algorithm doesn't match the manifold
    # throw a WrongManifoldException and be done with it
    throw(WrongManifoldException{typeof(m), M, Alg}())
end


# ### `NoAlgorithm`

"""
    NoAlgorithm(manifold)

A type that indicates no algorithm is to be used, essentially the equivalent
of `nothing`.

Stores a manifold within itself.
"""
struct NoAlgorithm{M <: Manifold} <: Algorithm{M} 
    m::M
end

NoAlgorithm() = NoAlgorithm(Planar()) # TODO: add a NoManifold or AutoManifold type?

manifold(a::NoAlgorithm) = a.m
# Maybe AutoManifold
# and then we have DD.format like materialization  