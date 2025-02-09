#=
# `Algorithm`s

An `Algorithm` is a type that describes the algorithm used to perform some [`Operation`](@ref).

An algorithm may be associated with one or many [`Manifold`](@ref)s.  It may either have the manifold as a field, or have it as a static parameter (e.g. `struct GEOS <: Manifold{Planar}`).


Algorithms are:
* Ways to perform an operation
* For example: LHuilier, Bessel, Ericsson for spherical area
* May be manifold agnostic (like simplification) or restrictive (like GEOS only works on planar, PROJ algorithm for arclength and area only works on geodesic)
* May or may not carry manifolds around, but manifold should always be accessible from manifold(alg) - it's not necessary that fixed manifold args can skip carrying the manifold around, eg in the case of Proj{Geodesic}.

=#

abstract type Algorithm{M <: Manifold} end

abstract type ManifoldIndependentAlgorithm{M <: Manifold} <: Algorithm{M} end

abstract type SingleManifoldAlgorithm{M <: Manifold} <: Algorithm{M} end

struct NoAlgorithm{M <: Manifold} <: Algorithm{M} 
    m::M
end

NoAlgorithm() = NoAlgorithm(Planar()) # TODO: add a NoManifold or AutoManifold type?
# Maybe AutoManifold
# and then we have DD.format like materialization  

function (Alg::Type{<: SingleManifoldAlgorithm{M}})(m::M; kwargs...) where {M}
    # successful - the algorithm is designed for this manifold
    # in this case, just return `Alg(; kwargs...)`
    return Alg(; kwargs...)
end

function (Alg::Type{<: ManifoldIndependentAlgorithm{M}})(m::Manifold; kwargs...) where {M}
    # this catches the case where the algorithm doesn't match the manifold
    # throw a WrongManifoldException and be done with it
    throw(WrongManifoldException{typeof(m), M, Alg}())
end

# for example

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