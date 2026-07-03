# # Prepared geometry

export prepare, Prepared, getprep, EdgeTree, EdgeTrees

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

## Materialization

`prepare` does not wrap your geometry in place — it *materializes* it into
GeometryOps' native layout: rings become vectors of coordinate tuples, and
every level of the geometry (each ring of a polygon, each polygon of a
multipolygon, …) becomes its own `Prepared` node with its own preparations and
cached extent.  Two consequences:

- Performance is uniform no matter where the input came from: a prepared
  ArchGDAL polygon and a prepared GeoJSON polygon query at the same speed,
  because both were converted to tuple storage up front.
- Preparedness survives decomposition.  `GI.getgeom` of a prepared
  multipolygon returns *prepared* polygons; `GI.getexterior` of a prepared
  polygon returns a *prepared* ring carrying its edge tree.  Any algorithm
  that tears a geometry apart through GeoInterface accessors keeps its
  acceleration all the way down.

Coordinate *number types are preserved* (a `Float32` polygon stays `Float32`);
only the memory layout changes.  Measure (`m`) coordinates are dropped, like
`GO.tuples`.  `Base.parent(prep)` returns the converted geometry — the
original object is not kept.

## Consuming preparations

There are three idioms for algorithm code, distinguished by what a miss costs:

```julia
# 1. Choose-path: use the index if present, else the plain algorithm.
#    Right when building the index costs as much as one unindexed query.
prep = GO.getprep(ring, GO.AbstractEdgeTree)
isnothing(prep) ? plain_path(ring) : indexed_path(ring, GO.edge_tree(prep))

# 2. Get-or-create: build ephemerally on a miss.
#    Right when the index pays for itself within a single operation.
prep = GO.getprep(ring, GO.AbstractEdgeTree) do
    GO.EdgeTree(ring)
end

# 3. Query-only, by concrete or abstract type.
idx = GO.getprep(ring, GO.EdgeTree)
```

## Extensibility

Every joint here is a function you can overload:

- `buildprep(spec, geom)` — how a spec becomes a preparation.  Defaults to
  `spec(geom)`, so types and closures already work.
- `default_preparations(trait, geom)` — what `prepare` builds at each node of
  the recursion when you don't say.
- `build_edge_tree(backend, ring)` — how an edge tree is built for a ring; add
  a method to plug in a new spatial tree backend.
- `edge_tree(p)` — accessor consumers use, so a custom `AbstractEdgeTree`
  subtype may store its tree however it likes.

A user-defined preparation is just `struct MyPrep <: GO.AbstractPreparation`
plus a constructor `MyPrep(geom)`.  It flows through
`prepare(geom; preps = (MyPrep,))` and `getprep(geom, MyPrep)` with no
registration step, because retrieval matches by `isa` — which also means
consumers can query by an abstract kind and accept any subtype a user hooks in.
=#

"""
    AbstractPreparation

Supertype for all geometry preparations — precomputed acceleration structures
stored in a [`Prepared`](@ref) wrapper and retrieved with [`getprep`](@ref).

Subtype this (or one of the abstract kinds below it, like
[`AbstractEdgeTree`](@ref)) and provide a constructor `MyPrep(geom)` to
define your own preparation.
"""
abstract type AbstractPreparation end

"""
    Prepared{T, G, P, E}

A geometry node bundled with a tuple of preparations `preps` and a cached
`extent`.  Construct one with [`prepare`](@ref).

The stored geometry `geom` is in GeometryOps' native layout: coordinate-tuple
storage whose *children are themselves `Prepared` nodes* — `GI.getgeom` of a
prepared multipolygon yields prepared polygons, `GI.getexterior` of a prepared
polygon yields a prepared ring carrying its edge tree.  `Prepared` implements
GeoInterface by forwarding to that storage, so it can be passed to anything
that accepts a geometry; `GI.extent` returns the cached extent.

Retrieve a preparation with [`getprep`](@ref); `Base.parent` returns the
converted geometry (the original input object is not kept).
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
# `Prepared` acts as its stored geometry everywhere.  We forward the core
# accessors that everything else (npoint, getring, gethole, coordinates, …)
# falls back to.  Because the storage's children are themselves `Prepared`
# nodes, forwarding `getgeom` is what makes preparedness survive decomposition.
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

`P` may be abstract — e.g. `getprep(ring, AbstractEdgeTree)` finds any
edge-tree preparation regardless of backend — or concrete.

The two-function form `getprep(f, geom, P)` returns `f()` on a miss, which
gives get-or-create when `f` builds:

```julia
prep = getprep(ring, AbstractEdgeTree) do
    EdgeTree(ring)    # built ephemerally when `ring` wasn't prepared
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

Build one preparation for `geom` from a spec selected for a node of the
`prepare` recursion.  The default is `spec(geom)`, so a preparation type
(`EdgeTree`) or a closure (`g -> EdgeTree(g; backend = STRtree)`) both work as
specs.  Overload this to make other spec objects buildable.
"""
buildprep(spec, geom) = spec(geom)

"""
    default_preparations(trait, geom)

The tuple of preparation specs [`prepare`](@ref) builds at a node of the
recursion when none are given.  Defaults to [`EdgeTree`](@ref) for linear
rings and nothing else (every `Prepared` node caches its extent regardless).
Overload on a trait — or on your geometry type — to change the default.
"""
default_preparations(trait, geom) = ()
default_preparations(::GI.LinearRingTrait, geom) = (EdgeTree,)

"""
    prepare(geom; preps = nothing)
    prepare(p::Prepared; preps::Tuple = ())

Materialize `geom` into GeometryOps' native layout and build a tree of
[`Prepared`](@ref) nodes over it: every ring becomes a vector of coordinate
tuples (number type preserved; `m` coordinates dropped), and every level —
ring, polygon, multi-geometry member — gets its own `Prepared` node with its
own preparations and cached extent.

`preps` controls what gets built at each node:

- `nothing` (default): [`default_preparations`](@ref) at every node — rings
  get a `NaturalIndex` [`EdgeTree`](@ref).
- a function `(trait, geom) -> Tuple` of specs, called at every node of the
  recursion.  [`EdgeTrees`](@ref) is a ready-made one for choosing the
  edge-tree backend: `prepare(poly; preps = EdgeTrees(HPR()))`.
- a tuple of specs: applied to the **top node only** (children still get
  defaults).  Each spec is built via [`buildprep`](@ref), so preparation
  types and closures both work.

On an already-`Prepared` input, build the given `preps` tuple against the
stored geometry and prepend them (so newly added preparations win
[`getprep`](@ref) lookups); nothing is re-materialized.

```julia
prep = prepare(poly)                                  # defaults
prep = prepare(poly; preps = EdgeTrees(STRtree))      # pick the tree backend
prep = prepare(poly; preps = (t, g) -> ())            # extent caches only
prep = prepare(poly; preps = (MyPrep,))               # custom prep on the top node
```
"""
function prepare(geom; preps = nothing)
    trait = GI.trait(geom)
    trait isa GI.AbstractGeometryTrait || throw(ArgumentError(
        "`prepare` requires a geometry (an object with a GeoInterface geometry trait), got $(typeof(geom))"))
    return _prepare(trait, geom, preps, GI.crs(geom), true)
end

function prepare(p::Prepared; preps::Tuple = ())
    isempty(preps) && return p
    built = map(spec -> buildprep(spec, p.geom), preps)
    return Prepared(p.geom, (built..., p.preps...), p.extent)
end

# Which specs apply at a node: `nothing` = defaults everywhere; a tuple =
# top node only; anything callable = `preps(trait, geom)` at every node.
_node_preps(::Nothing, trait, geom, istop) = default_preparations(trait, geom)
_node_preps(specs::Tuple, trait, geom, istop) = istop ? specs : default_preparations(trait, geom)
_node_preps(f, trait, geom, istop) = f(trait, geom)

# Build the preps for a node and close the `Prepared` shell over it.
function _wrap(trait, geom, extent, preps, istop)
    built = map(spec -> buildprep(spec, geom), _node_preps(preps, trait, geom, istop))
    return Prepared{typeof(trait), typeof(geom), typeof(built), typeof(extent)}(geom, built, extent)
end

# Leaf storage: coordinate tuples with the input's number types.  Measures are
# dropped (a 3-tuple must mean x/y/z to GeoInterface), matching `GO.tuples`.
function _tuple_points(geom)
    if GI.is3d(geom)
        return [(GI.x(p), GI.y(p), GI.z(p)) for p in GI.getpoint(geom)]
    else
        return [(GI.x(p), GI.y(p)) for p in GI.getpoint(geom)]
    end
end

# Points: stored as a bare coordinate tuple.
function _prepare(trait::GI.PointTrait, geom, preps, crs, istop)
    pt = GI.is3d(geom) ? (GI.x(geom), GI.y(geom), GI.z(geom)) : (GI.x(geom), GI.y(geom))
    ext = Extents.Extent(X = (pt[1], pt[1]), Y = (pt[2], pt[2]))
    return _wrap(trait, pt, ext, preps, istop)
end

# Curves and multipoints: a GeoInterface wrapper over tuple storage.
function _prepare(trait::Union{GI.AbstractCurveTrait, GI.MultiPointTrait}, geom, preps, crs, istop)
    pts = _tuple_points(geom)
    T = GI.geointerface_geomtype(trait)
    ext = GI.extent(T(pts; crs))
    return _wrap(trait, T(pts; crs, extent = ext), ext, preps, istop)
end

# Everything with geometry children (polygons, multi-geometries, collections):
# recurse, so the stored children are themselves `Prepared` nodes.  The
# `map(identity, …)` tightens the child vector's eltype (heterogeneous
# collections get a small union).
function _prepare(trait::GI.AbstractGeometryTrait, geom, preps, crs, istop)
    children = map(identity, [_prepare(GI.trait(c), c, preps, crs, false) for c in GI.getgeom(geom)])
    ext = isempty(children) ? GI.extent(geom) :
        mapreduce(c -> c.extent, Extents.union, children)
    T = GI.geointerface_geomtype(trait)
    return _wrap(trait, T(children; crs, extent = ext), ext, preps, istop)
end

# ## Edge-tree preparations

"""
    AbstractEdgeTree <: AbstractPreparation

The preparation *kind* for "a spatial index over the edges of a curve".  It
lives on the prepared **ring**, so any consumer holding a ring — the planar
point-in-polygon test, line/curve processes, clipping — can discover it with
`getprep(ring, AbstractEdgeTree)`.  Consumers access the tree through
[`edge_tree`](@ref), so any subtype with any SpatialTreeInterface-compatible
tree hooks in automatically.
"""
abstract type AbstractEdgeTree <: AbstractPreparation end

"""
    edge_tree(p::AbstractEdgeTree)

The spatial index stored by an edge-tree preparation.  Defaults to `p.tree`.
"""
edge_tree(p::AbstractEdgeTree) = p.tree

"""
    EdgeTree(curve; backend = NaturalIndex)

A spatial index over the edge extents of a curve — the default
[`AbstractEdgeTree`](@ref), built for every linear ring by `prepare`.  Edge
`i` of an `n`-point ring runs from point `i` to point `i + 1`, except that an
unclosed ring's last edge wraps back to point `1`, matching the
implicit-closure semantics of the plain point-in-polygon algorithm.

`backend` picks the tree via [`build_edge_tree`](@ref): `NaturalIndex`
(default), `STRtree`, a `FlexibleRTrees` bulk-load algorithm (`STR()`,
`HPR()`, `Unsorted()`), or any callable `curve -> spatial tree`.
"""
struct EdgeTree{T} <: AbstractEdgeTree
    tree::T
end

function EdgeTree(geom; backend = NaturalIndex)
    GI.trait(geom) isa GI.AbstractCurveTrait || throw(ArgumentError(
        "`EdgeTree` requires a curve (linear ring or line string), got $(typeof(GI.trait(geom)))"))
    tree = build_edge_tree(backend, geom)
    return EdgeTree{typeof(tree)}(tree)
end

"""
    EdgeTrees(backend = NaturalIndex)

A ready-made `preps` selector for [`prepare`](@ref) that puts an
[`EdgeTree`](@ref) with the given `backend` on every linear ring:

```julia
prep = prepare(poly; preps = EdgeTrees(STRtree))
prep = prepare(poly; preps = EdgeTrees(FlexibleRTrees.HPR()))
```
"""
struct EdgeTrees{B}
    backend::B
end
EdgeTrees() = EdgeTrees(NaturalIndex)
(s::EdgeTrees)(trait, geom) =
    trait isa GI.LinearRingTrait ? (g -> EdgeTree(g; backend = s.backend),) : ()

"""
    build_edge_tree(backend, curve)

Build a SpatialTreeInterface-compatible spatial index over the edges of
`curve`, where the tree's leaf indices are edge indices (edge `i` runs from
point `i` to point `i + 1`, wrapping to point `1` from the last point of an
unclosed ring).

Methods exist for `NaturalIndex`, `STRtree`, and `FlexibleRTrees` bulk-load
algorithms; the fallback calls `backend(curve)`, so any callable works.  Add a
method to plug in a new tree type — it only needs to implement
SpatialTreeInterface.
"""
build_edge_tree(backend, ring) = backend(ring)
build_edge_tree(::Type{<:NaturalIndex}, ring) = NaturalIndex(_ring_edge_extents(ring))
build_edge_tree(::Type{<:STRtree}, ring) = STRtree(_ring_edge_extents(ring))
build_edge_tree(alg::FlexibleRTrees.BulkLoadAlgorithm, ring) =
    FlexibleRTrees.RTree(alg, _ring_edge_extents(ring))

# Extents of the ring's edges, adding the implicit closing edge when the ring
# is unclosed — mirrors how `_point_filled_curve_orientation` walks edges.
# Coordinate number types are preserved.
function _ring_edge_extents(ring)
    n = GI.npoint(ring)
    closed = equals(GI.getpoint(ring, 1), GI.getpoint(ring, n))
    nedges = closed ? n - 1 : n
    return [begin
        p1 = GI.getpoint(ring, i)
        p2 = GI.getpoint(ring, i == n ? 1 : i + 1)
        Extents.Extent(
            X = minmax(GI.x(p1), GI.x(p2)),
            Y = minmax(GI.y(p1), GI.y(p2)),
        )
    end for i in 1:nedges]
end
