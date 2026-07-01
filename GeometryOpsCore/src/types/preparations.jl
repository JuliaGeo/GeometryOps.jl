#=
# Prepared geometry

Interface types for prepared geometry: a GeoInterface-compatible wrapper that carries a
tuple of immutable, eagerly built acceleration structures ("preparations"), retrieved by
*capability* trait — algorithms only care that a geometry *has* an edge index, not how it
was built. Design: `docs/plans/2026-07-01-prepared-geometry-design.md` (repo root).

The consumption interface lives here in GeometryOpsCore so third-party packages can
provide preparations without depending on GeometryOps. The `prepare` builder and the
concrete preparations live in GeometryOps (extents are manifold-aware).
=#

export AbstractPreparation, AbstractPreparationTrait,
    SpatialIndexLike, SpatialEdgeIndexLike, PointInAreaLike,
    preptrait, PreparationSpec, Prepared, getprep, prepare

"""
    AbstractPreparation

Supertype of built, immutable acceleration structures attached to exactly one geometry via
[`Prepared`](@ref). Concrete subtypes implement [`preptrait`](@ref).
"""
abstract type AbstractPreparation end

"""
    AbstractPreparationTrait

Capability traits for preparations: what a preparation can answer, not what it is.
Retrieval via `get(prepared, trait)` returns the first preparation with that capability.
"""
abstract type AbstractPreparationTrait end

"Capability: an extent tree over child elements (rings of a polygon, polygons of a multi)."
struct SpatialIndexLike <: AbstractPreparationTrait end
"Capability: an extent tree over the segments/edges of a curve or geometry."
struct SpatialEdgeIndexLike <: AbstractPreparationTrait end
"Capability: a point-in-polygonal-area locator."
struct PointInAreaLike <: AbstractPreparationTrait end

"""
    preptrait(prep::AbstractPreparation)::AbstractPreparationTrait

The capability of `prep`. Every concrete preparation implements this.
"""
function preptrait end

"""
    PreparationSpec

Supertype of the *configs* passed to [`prepare`](@ref)'s `preps` tuple. A spec declares
which GeoInterface trait it attaches to via `appliesto(spec)` and how to build via
`buildprep(spec, manifold, geom) -> Union{AbstractPreparation, Nothing}` (`nothing` means
"nothing to index here", e.g. an empty geometry). Both are extended (not exported) via
`import GeometryOpsCore: appliesto, buildprep`.
"""
abstract type PreparationSpec end

function appliesto end
function buildprep end

"""
    Prepared(geom, manifold::Manifold, preps::Tuple, extent::Extents.Extent)

A geometry `geom` carrying preparations `preps`. Forwards GeoInterface, so it works
anywhere a geometry does. `extent` is always cached (manifold-aware: 3D on `Spherical`).
The wrapper and its preparation set are immutable; a miss on `get` means the caller takes
its unindexed path. Children of a `Prepared` built by `prepare` are themselves `Prepared`.
"""
struct Prepared{T <: GI.AbstractTrait, M <: Manifold, G,
                P <: Tuple{Vararg{AbstractPreparation}}, E <: Extents.Extent}
    geom::G
    manifold::M
    preps::P
    extent::E
    function Prepared(geom::G, m::M, preps::P, extent::E) where
            {M <: Manifold, G, P <: Tuple{Vararg{AbstractPreparation}}, E <: Extents.Extent}
        t = GI.trait(geom)
        t === nothing && throw(ArgumentError("`Prepared` requires a GeoInterface geometry; got $(typeof(geom))"))
        new{typeof(t), M, G, P, E}(geom, m, preps, extent)
    end
end

Base.parent(p::Prepared) = p.geom
manifold(p::Prepared) = p.manifold

"""
    get(p::Prepared, t::AbstractPreparationTrait) -> Union{AbstractPreparation, Nothing}
    get(geom, t::AbstractPreparationTrait) -> nothing

The first preparation on `p` with capability `t`, or `nothing`. Plain geometries always
return `nothing`, so call sites need no special-casing.
"""
Base.get(p::Prepared, t::AbstractPreparationTrait) = _findprep(t, p.preps)
Base.get(@nospecialize(_), ::AbstractPreparationTrait) = nothing

_findprep(::AbstractPreparationTrait, ::Tuple{}) = nothing
_findprep(t::AbstractPreparationTrait, preps::Tuple) =
    preptrait(first(preps)) === t ? first(preps) : _findprep(t, Base.tail(preps))

"""
    getprep(m::Manifold, geom, t::AbstractPreparationTrait) -> Union{AbstractPreparation, Nothing}

Manifold-checked retrieval: like `get(geom, t)` but a manifold mismatch is a miss —
preparations bake in manifold-specific state (kernel point types, 3D extents), so a
`Planar` edge tree must never be served to a `Spherical` algorithm.
"""
getprep(m::Manifold, p::Prepared, t::AbstractPreparationTrait) =
    manifold(p) === m ? get(p, t) : nothing
getprep(::Manifold, @nospecialize(_), ::AbstractPreparationTrait) = nothing

"""
    prepare(...)

Generic entry point for prepared-geometry optimizations. GeometryOps implements
`prepare(geom; preps, manifold)` (the wrapper builder) and `prepare(alg, geom)` methods
for algorithms (e.g. `RelateNG`).
"""
function prepare end
