# # Prepared geometry

export prepare, Prepared, getprep, RingEdgeTrees

#=
## What is prepared geometry?

When you run many operations against the same geometry — thousands of
point-in-polygon tests against one country border, say — most of the work each
call does is rediscovering structure of that geometry: its extent, where its
edges live in space, and so on.  A *prepared* geometry computes that structure
once and carries it along, so every subsequent operation can reuse it.

```julia
import GeometryOps as GO, GeoInterface as GI

poly = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
prep = GO.prepare(poly)

GO.contains(prep, (5.0, 5.0))   # uses the edge index and cached extent
```

`Prepared` is a transparent GeoInterface wrapper: it behaves exactly like the
geometry it wraps everywhere in GeometryOps (and in any GeoInterface-compatible
package), and functions that know about preparations can look them up with
[`getprep`](@ref).

## Consuming preparations

There are three idioms for algorithm code, distinguished by what a miss costs:

```julia
# 1. Choose-path: use the index if present, else the plain algorithm.
#    Right when building the index costs as much as one unindexed query.
trees = GO.getprep(poly, GO.AbstractRingEdgeTrees)
isnothing(trees) ? plain_path(poly) : indexed_path(poly, trees)

# 2. Get-or-create: build ephemerally on a miss.
#    Right when the index pays for itself within a single operation.
trees = GO.getprep(poly, GO.AbstractRingEdgeTrees) do
    GO.RingEdgeTrees(poly)
end

# 3. Query-only, by concrete or abstract type.
idx = GO.getprep(poly, GO.RingEdgeTrees)
```

## Extensibility

Every joint here is a function you can overload:

- `buildprep(spec, geom)` — how a spec in `prepare(geom; preps = (...,))` becomes a
  preparation.  Defaults to `spec(geom)`, so types and closures already work.
- `default_preparations(trait, geom)` — what `prepare` builds when you don't say.
- `build_edge_tree(backend, ring)` — how an edge tree is built for a ring; add a
  method to plug in a new spatial tree backend.
- `exterior_tree(p)` / `hole_trees(p)` — accessors consumers use, so a custom
  `AbstractRingEdgeTrees` subtype may store its trees however it likes.

A user-defined preparation is just `struct MyPrep <: GO.AbstractPreparation`
plus a constructor `MyPrep(geom)`.  It flows through `prepare(geom; preps = (MyPrep,))`
and `getprep(geom, MyPrep)` with no registration step, because retrieval matches
by `isa` — which also means consumers can query by an abstract kind and accept
any subtype a user hooks in.
=#

"""
    AbstractPreparation

Supertype for all geometry preparations — precomputed acceleration structures
stored in a [`Prepared`](@ref) wrapper and retrieved with [`getprep`](@ref).

Subtype this (or one of the abstract kinds below it, like
[`AbstractRingEdgeTrees`](@ref)) and provide a constructor `MyPrep(geom)` to
define your own preparation.
"""
abstract type AbstractPreparation end

"""
    Prepared{T, G, P, E}

A geometry `geom` bundled with a tuple of preparations `preps` and a cached
`extent`.  Construct one with [`prepare`](@ref).

`Prepared` implements GeoInterface by forwarding to the parent geometry, so it
can be passed to anything that accepts a geometry.  `GI.extent` returns the
cached extent.  Children (rings, member geometries, …) are returned **plain** —
preparations live only on the wrapper they were built for.

Retrieve a preparation with [`getprep`](@ref); get the parent geometry back
with `Base.parent`.
"""
struct Prepared{T <: GI.AbstractGeometryTrait, G, P <: Tuple, E}
    geom::G
    preps::P
    extent::E
end

function Prepared(geom, preps::Tuple, extent)
    geom isa Prepared && throw(ArgumentError(
        "geometry is already `Prepared`; use `prepare(p; preps)` to add preparations"))
    trait = GI.trait(geom)
    trait isa GI.AbstractGeometryTrait || throw(ArgumentError(
        "`Prepared` requires a geometry (an object with a GeoInterface geometry trait), got $(typeof(geom))"))
    return Prepared{typeof(trait), typeof(geom), typeof(preps), typeof(extent)}(geom, preps, extent)
end

Base.parent(p::Prepared) = p.geom

function Base.show(io::IO, p::Prepared{T}) where T
    print(io, "Prepared", "{", nameof(T), "}")
    names = join((nameof(typeof(prep)) for prep in p.preps), ", ")
    print(io, " with cached extent", isempty(p.preps) ? "" : " and preparations ($names)")
end

# ## GeoInterface forwarding
#
# `Prepared` acts as its parent everywhere.  We forward the core accessors that
# everything else (npoint, getring, gethole, coordinates, …) falls back to.
GI.isgeometry(::Type{<:Prepared}) = true
GI.geomtrait(::Prepared{T}) where T = T()

GI.ngeom(t::GI.AbstractGeometryTrait, p::Prepared) = GI.ngeom(t, p.geom)
GI.getgeom(t::GI.AbstractGeometryTrait, p::Prepared) = GI.getgeom(t, p.geom)
GI.getgeom(t::GI.AbstractGeometryTrait, p::Prepared, i) = GI.getgeom(t, p.geom, i)
GI.ncoord(t::GI.AbstractGeometryTrait, p::Prepared) = GI.ncoord(t, p.geom)
GI.getcoord(t::GI.AbstractPointTrait, p::Prepared, i) = GI.getcoord(t, p.geom, i)
GI.is3d(t::GI.AbstractGeometryTrait, p::Prepared) = GI.is3d(t, p.geom)
GI.ismeasured(t::GI.AbstractGeometryTrait, p::Prepared) = GI.ismeasured(t, p.geom)
GI.crs(t::GI.AbstractGeometryTrait, p::Prepared) = GI.crs(p.geom)
# Point-trait disambiguators (GeoInterface defines these on `Any` for points):
GI.ngeom(::GI.AbstractPointTrait, ::Prepared) = 0
GI.getgeom(::GI.AbstractPointTrait, ::Prepared) = nothing
GI.getgeom(::GI.AbstractPointTrait, ::Prepared, i) = nothing
# The cached extent is authoritative:
GI.extent(::GI.AbstractGeometryTrait, p::Prepared) = p.extent
Extents.extent(p::Prepared) = p.extent

# ## Retrieval

"""
    getprep(geom, P::Type)
    getprep(f, geom, P::Type)

Return the first preparation stored in `geom` that `isa P`, or `nothing` if
there is none.  Plain (un-prepared) geometries always return `nothing`, so
consumer code needs no special-casing.

`P` may be abstract — e.g. `getprep(geom, AbstractRingEdgeTrees)` finds any
edge-tree preparation regardless of backend — or concrete.

The two-function form `getprep(f, geom, P)` returns `f()` on a miss, which
gives get-or-create when `f` builds:

```julia
trees = getprep(poly, AbstractRingEdgeTrees) do
    RingEdgeTrees(poly)    # built ephemerally when `poly` wasn't prepared
end
```
"""
getprep(geom, ::Type{P}) where P = nothing
getprep(p::Prepared, ::Type{P}) where P = _first_prep(P, p.preps)

function getprep(f, geom, ::Type{P}) where P
    prep = getprep(geom, P)
    return isnothing(prep) ? f() : prep
end

@inline _first_prep(::Type{P}, preps::Tuple) where P =
    first(preps) isa P ? first(preps) : _first_prep(P, Base.tail(preps))
@inline _first_prep(::Type{P}, ::Tuple{}) where P = nothing

# ## Building

"""
    buildprep(spec, geom)

Build one preparation for `geom` from a spec listed in `prepare`'s `preps`
tuple.  The default is `spec(geom)`, so a preparation type (`RingEdgeTrees`) or
a closure (`g -> RingEdgeTrees(g; tree = STRtree)`) both work as specs.
Overload this to make other spec objects buildable.
"""
buildprep(spec, geom) = spec(geom)

"""
    default_preparations(trait, geom)

The tuple of preparation specs [`prepare`](@ref) builds for a geometry when
none are given.  Defaults to [`RingEdgeTrees`](@ref) for polygons and nothing
else (every `Prepared` caches its extent regardless).  Overload on a trait —
or on your geometry type — to change the default.
"""
default_preparations(trait, geom) = ()
default_preparations(::GI.PolygonTrait, geom) = (RingEdgeTrees,)

"""
    prepare(geom; preps = default_preparations(GI.trait(geom), geom))
    prepare(p::Prepared; preps = ())

Build a [`Prepared`](@ref) wrapper around `geom`: compute its extent, build
each preparation in `preps` (via [`buildprep`](@ref)), and bundle everything.

On an already-`Prepared` input, build the given `preps` against the parent
geometry and prepend them (so newly added preparations win [`getprep`](@ref)
lookups); the cached extent is reused.

```julia
prep = prepare(poly)                                      # defaults
prep = prepare(poly; preps = ())                          # extent cache only
prep = prepare(poly; preps = (g -> RingEdgeTrees(g; tree = STRtree),))
```
"""
function prepare(geom; preps = default_preparations(GI.trait(geom), geom))
    built = map(spec -> buildprep(spec, geom), preps)
    return Prepared(geom, built, GI.extent(geom))
end
function prepare(p::Prepared; preps = ())
    isempty(preps) && return p
    built = map(spec -> buildprep(spec, p.geom), preps)
    return Prepared(p.geom, (built..., p.preps...), p.extent)
end

# ## Edge-tree preparations

"""
    AbstractRingEdgeTrees <: AbstractPreparation

The preparation *kind* for "one spatial index over the edges of each ring of a
polygon".  Consumers (like the planar point-in-polygon test) query this
abstract type and access trees through [`exterior_tree`](@ref) and
[`hole_trees`](@ref), so any subtype with any SpatialTreeInterface-compatible
tree hooks in automatically.
"""
abstract type AbstractRingEdgeTrees <: AbstractPreparation end

"""
    exterior_tree(p::AbstractRingEdgeTrees)

The edge tree of the polygon's exterior ring.  Defaults to `p.exterior`.
"""
exterior_tree(p::AbstractRingEdgeTrees) = p.exterior

"""
    hole_trees(p::AbstractRingEdgeTrees)

An indexable collection of the edge trees of the polygon's holes, in
`GI.gethole` order.  Defaults to `p.holes`.
"""
hole_trees(p::AbstractRingEdgeTrees) = p.holes

"""
    RingEdgeTrees(polygon; tree = NaturalIndex)

One spatial index over the edge extents of each ring of `polygon` — the
default [`AbstractRingEdgeTrees`](@ref).  Edge `i` of an `n`-point ring runs
from point `i` to point `i + 1`, except that an unclosed ring's last edge wraps
back to point `1`, matching the implicit-closure semantics of the plain
point-in-polygon algorithm.

`tree` picks the backend via [`build_edge_tree`](@ref): `NaturalIndex`
(default), `STRtree`, or any callable `ring -> spatial tree`.
"""
struct RingEdgeTrees{TE, TH} <: AbstractRingEdgeTrees
    exterior::TE
    holes::Vector{TH}
end

function RingEdgeTrees(polygon; tree = NaturalIndex)
    trait = GI.trait(polygon)
    trait isa GI.PolygonTrait || throw(ArgumentError(
        "`RingEdgeTrees` requires a polygon, got $(typeof(trait)); use it inside `prepare` on each polygon"))
    exterior = build_edge_tree(tree, GI.getexterior(polygon))
    holes = [build_edge_tree(tree, hole) for hole in GI.gethole(polygon)]
    return RingEdgeTrees(exterior, holes)
end

"""
    build_edge_tree(backend, ring)

Build a SpatialTreeInterface-compatible spatial index over the edges of `ring`,
where the tree's leaf indices are edge indices (edge `i` runs from point `i` to
point `i + 1`, wrapping to point `1` from the last point of an unclosed ring).

Methods exist for `NaturalIndex` and `STRtree`; the fallback calls
`backend(ring)`, so any callable works.  Add a method to plug in a new tree
type — it only needs to implement SpatialTreeInterface.
"""
build_edge_tree(backend, ring) = backend(ring)
build_edge_tree(::Type{<:NaturalIndex}, ring) = NaturalIndex(_ring_edge_extents(ring))
build_edge_tree(::Type{<:STRtree}, ring) = STRtree(_ring_edge_extents(ring))

# Extents of the ring's edges, adding the implicit closing edge when the ring
# is unclosed — mirrors how `_point_filled_curve_orientation` walks edges.
function _ring_edge_extents(ring)
    n = GI.npoint(ring)
    closed = equals(GI.getpoint(ring, 1), GI.getpoint(ring, n))
    nedges = closed ? n - 1 : n
    return [begin
        p1 = GI.getpoint(ring, i)
        p2 = GI.getpoint(ring, i == n ? 1 : i + 1)
        Extents.Extent(
            X = minmax(Float64(GI.x(p1)), Float64(GI.x(p2))),
            Y = minmax(Float64(GI.y(p1)), Float64(GI.y(p2))),
        )
    end for i in 1:nedges]
end
