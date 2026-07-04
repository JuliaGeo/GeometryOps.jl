# # Prepared geometry

export prepare, Prepared, getprep, hasprep, EdgeTree

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

Coordinate *number types are preserved* (a `Float32` polygon stays `Float32`),
and so are points already in a native representation (a curve of
`UnitSphericalPoint`s keeps them); only the memory layout changes.  Measure
(`m`) coordinates are dropped, like `GO.tuples`.  Materialized linear rings
are always **closed** (the first point is repeated at the end if the input
was unclosed), so preparations built against materialized storage never need
to handle an implicit closing edge.  `Base.parent(prep)` returns the
converted geometry — the original object is not kept.

`prepare` takes the manifold as its first argument, like other GeometryOps
functions: `prepare(geom)` means `prepare(Planar(), geom)`.  The manifold
decides what gets built where — e.g. edge trees are a planar default —
and flows into every preparation via [`buildprep`](@ref) and
[`build_edge_tree`](@ref).

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

- `buildprep(manifold, spec, geom)` — how a spec becomes a preparation.
  Defaults to `spec(geom)`, so types and closures already work; overload the
  manifold-taking form for manifold-aware preparations.
- `appliesto(spec, trait, istop)` — where a spec given in the `preps` tuple
  applies during the recursion (default: the top node only; `EdgeTree`
  declares itself for every curve).
- `default_preparations(manifold, trait, geom)` — what `prepare` builds at
  each node of the recursion when you don't say.
- `build_edge_tree(manifold, backend, ring)` — how an edge tree is built for
  a ring; add a method to plug in a new spatial tree backend.
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

Supertype for geometry preparations — precomputed acceleration structures
stored in a [`Prepared`](@ref) wrapper and retrieved with [`getprep`](@ref).
Subtype it (or an abstract kind like [`AbstractEdgeTree`](@ref)) and provide
a constructor `MyPrep(geom)` to define your own.
"""
abstract type AbstractPreparation end

"""
    Prepared{T, G, P, E}

A geometry node bundled with a tuple of preparations and a cached extent.
Construct with [`prepare`](@ref); retrieve preparations with
[`getprep`](@ref).

The stored geometry is in GeometryOps' native layout and its children are
themselves `Prepared` nodes, so preparedness survives decomposition —
`GI.getexterior` of a prepared polygon is a prepared ring.  `Prepared`
implements GeoInterface by forwarding to that storage; `GI.extent` returns
the cached extent; `Base.parent` returns the converted geometry (the
original input object is not kept).

Invariant consumers may rely on: every preparation stored in a `Prepared`
node was built against the stored (materialized) geometry, and materialized
linear rings are closed.  An edge tree retrieved from a prepared ring
therefore indexes exactly the `GI.npoint(ring) - 1` explicit edges.
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
# Forward the core accessors everything else falls back to.  The storage's
# children are themselves `Prepared`, so `getgeom` yields prepared nodes.
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

The lookup is resolved at compile time from the preparation tuple's type —
the generated body is a constant field access (or `nothing`), so there is
no runtime search.
"""
getprep(geom, ::Type{P}) where P = nothing
@generated function getprep(p::Prepared{T, G, PT, E}, ::Type{P}) where {T, G, PT, E, P}
    i = findfirst(t -> t <: P, collect(PT.parameters))
    return isnothing(i) ? :(nothing) : :(p.preps[$i])
end

function getprep(f, geom, ::Type{P}) where P
    prep = getprep(geom, P)
    return isnothing(prep) ? f() : prep
end

"""
    hasprep(geom, P::Type)::Bool

Whether `geom` stores a preparation that `isa P` — the boolean companion to
[`getprep`](@ref).  Like `getprep`, this looks at `geom`'s own node only:
`hasprep(prepared_polygon, AbstractEdgeTree)` is `false` because edge trees
live on the polygon's *rings*.
"""
hasprep(geom, ::Type{P}) where P = !isnothing(getprep(geom, P))

# Strip a `Prepared` shell (a no-op on anything else); hot kernels walk raw
# point storage directly after retrieving the preparations they need.
_unwrap_prepared(g) = g
_unwrap_prepared(p::Prepared) = parent(p)

# ## Building

"""
    buildprep(manifold, spec, geom)
    buildprep(spec, geom)

Build one preparation for `geom` from a spec.  Defaults to `spec(geom)`, so
a preparation type (`EdgeTree`) or a closure both work as specs; overload to
make other spec objects buildable.  Manifold-aware preparations overload the
three-argument form — the two-argument form is the manifold-oblivious
fallback it reaches by default.
"""
buildprep(m::Manifold, spec, geom) = buildprep(spec, geom)
buildprep(spec, geom) = spec(geom)

"""
    default_preparations(manifold, trait, geom)

The tuple of preparation specs [`prepare`](@ref) builds at a node of the
recursion when none are given: [`EdgeTree`](@ref) for curves (linear rings
and line strings) on the planar manifold, nothing else.  Overload on the
manifold and/or trait to change the default.
"""
default_preparations(m::Manifold, trait, geom) = ()
default_preparations(::Planar, ::GI.AbstractCurveTrait, geom) = (EdgeTree,)

"""
    prepare([manifold::Manifold], geom; preps = nothing)
    prepare([manifold::Manifold], p::Prepared; preps::Tuple = ())

Materialize `geom` into GeometryOps' native layout — coordinate-tuple rings
(closed, number type preserved, `m` coordinates dropped) — and build a tree
of [`Prepared`](@ref) nodes over it, one per level (ring, polygon,
multi-geometry member), each with its own preparations and cached extent.

The `manifold` defaults to `Planar()` and flows into every preparation —
it decides both what gets built by default and how (see
[`default_preparations`](@ref), [`buildprep`](@ref),
[`build_edge_tree`](@ref)).

`preps` controls what gets built at each node:

- `nothing` (default): [`default_preparations`](@ref) at every node — on
  `Planar()`, rings and line strings get a `NaturalIndex`
  [`EdgeTree`](@ref).
- a tuple of specs: each spec applies at the nodes it declares via
  [`appliesto`](@ref) — the top node only by default, every curve for
  `EdgeTree` (bare or curried, `EdgeTree(HPR())` picking the backend) —
  and nodes where no given spec applies still get defaults.  Each spec is
  built via [`buildprep`](@ref), so preparation types and closures both
  work.
- a function `(trait, geom) -> Tuple` of specs, called at every node of the
  recursion, overriding defaults everywhere.

On an already-`Prepared` input, build the given `preps` tuple against the
stored geometry and prepend them (so newly added preparations win
[`getprep`](@ref) lookups); nothing is re-materialized.

```julia
prep = prepare(poly)                                   # defaults
prep = prepare(poly; preps = (EdgeTree(STRtree),))     # pick the tree backend
prep = prepare(poly; preps = (t, g) -> ())             # extent caches only
prep = prepare(poly; preps = (MyPrep,))                # custom prep on the top node
```
"""
prepare(geom; preps = nothing) = prepare(Planar(), geom; preps)

function prepare(m::Manifold, geom; preps = nothing)
    trait = GI.trait(geom)
    trait isa GI.AbstractGeometryTrait || throw(ArgumentError(
        "`prepare` requires a geometry (an object with a GeoInterface geometry trait), got $(typeof(geom))"))
    return _prepare(m, trait, geom, preps, GI.crs(geom), true)
end

prepare(p::Prepared; preps::Tuple = ()) = prepare(Planar(), p; preps)

function prepare(m::Manifold, p::Prepared; preps::Tuple = ())
    isempty(preps) && return p
    built = map(spec -> buildprep(m, spec, p.geom), preps)
    return Prepared(p.geom, (built..., p.preps...), p.extent)
end

# Which specs apply at a node: `nothing` = defaults everywhere; a tuple =
# each spec where it declares itself via `appliesto`, with defaults filling
# the nodes where none applies; anything callable = `preps(trait, geom)` at
# every node.
_node_preps(m, ::Nothing, trait, geom, istop) = default_preparations(m, trait, geom)
function _node_preps(m, specs::Tuple, trait, geom, istop)
    applicable = _filter_specs(specs, trait, istop)
    return isempty(applicable) ? default_preparations(m, trait, geom) : applicable
end
_node_preps(m, f, trait, geom, istop) = f(trait, geom)

# Tuple-recursive filter, so the applicable-spec tuple stays concretely typed.
_filter_specs(specs::Tuple, trait, istop) =
    appliesto(first(specs), trait, istop) ?
        (first(specs), _filter_specs(Base.tail(specs), trait, istop)...) :
        _filter_specs(Base.tail(specs), trait, istop)
_filter_specs(::Tuple{}, trait, istop) = ()

# Build the preps for a node and close the `Prepared` shell over it.
function _wrap(m::Manifold, trait, geom, extent, preps, istop)
    built = map(spec -> buildprep(m, spec, geom), _node_preps(m, preps, trait, geom, istop))
    return Prepared{typeof(trait), typeof(geom), typeof(built), typeof(extent)}(geom, built, extent)
end

# Leaf storage: coordinate tuples with the input's number types — except
# points already in a native representation (`UnitSphericalPoint`), which are
# stored as-is.  Measures are dropped (a 3-tuple must mean x/y/z to
# GeoInterface), matching `GO.tuples`.
function _materialize_points(geom)
    if GI.npoint(geom) > 0 && first(GI.getpoint(geom)) isa UnitSpherical.UnitSphericalPoint
        return collect(GI.getpoint(geom))
    elseif GI.is3d(geom)
        return [(GI.x(p), GI.y(p), GI.z(p)) for p in GI.getpoint(geom)]
    else
        return [(GI.x(p), GI.y(p)) for p in GI.getpoint(geom)]
    end
end

# Points: stored as a bare coordinate tuple (or kept as a `UnitSphericalPoint`).
function _prepare(m::Manifold, trait::GI.PointTrait, geom, preps, crs, istop)
    pt = geom isa UnitSpherical.UnitSphericalPoint ? geom :
        GI.is3d(geom) ? (GI.x(geom), GI.y(geom), GI.z(geom)) : (GI.x(geom), GI.y(geom))
    ext = pt isa UnitSpherical.UnitSphericalPoint ?
        Extents.Extent(X = (pt[1], pt[1]), Y = (pt[2], pt[2]), Z = (pt[3], pt[3])) :
        Extents.Extent(X = (pt[1], pt[1]), Y = (pt[2], pt[2]))
    return _wrap(m, trait, pt, ext, preps, istop)
end

# Curves and multipoints: a GeoInterface wrapper over materialized point storage.
function _prepare(m::Manifold, trait::Union{GI.AbstractCurveTrait, GI.MultiPointTrait}, geom, preps, crs, istop)
    pts = _materialize_points(geom)
    # Materialized rings are always closed — the invariant preparations and
    # their consumers rely on (see the `Prepared` docstring).
    if trait isa GI.LinearRingTrait && !isempty(pts) && first(pts) != last(pts)
        push!(pts, first(pts))
    end
    T = GI.geointerface_geomtype(trait)
    ext = GI.extent(T(pts; crs))
    return _wrap(m, trait, T(pts; crs, extent = ext), ext, preps, istop)
end

# Polygons: the children are rings *by construction*, even when the backend
# types them as line strings (GeoJSON does) — materialize them as linear
# rings so they pick up ring defaults like `EdgeTree`.
function _prepare(m::Manifold, trait::GI.PolygonTrait, geom, preps, crs, istop)
    children = map(identity, [_prepare(m, GI.LinearRingTrait(), r, preps, crs, false) for r in GI.getring(geom)])
    ext = isempty(children) ? GI.extent(geom) :
        mapreduce(c -> c.extent, Extents.union, children)
    return _wrap(m, trait, GI.Polygon(children; crs, extent = ext), ext, preps, istop)
end

# Multi-geometries and collections: recurse, so the stored children are
# themselves `Prepared`.  `map(identity, …)` tightens the child vector's
# eltype (heterogeneous collections get a small union).
function _prepare(m::Manifold, trait::GI.AbstractGeometryTrait, geom, preps, crs, istop)
    children = map(identity, [_prepare(m, GI.trait(c), c, preps, crs, false) for c in GI.getgeom(geom)])
    ext = isempty(children) ? GI.extent(geom) :
        mapreduce(c -> c.extent, Extents.union, children)
    T = GI.geointerface_geomtype(trait)
    return _wrap(m, trait, T(children; crs, extent = ext), ext, preps, istop)
end

# ## Edge-tree preparations

"""
    AbstractEdgeTree <: AbstractPreparation

The preparation *kind* for a spatial index over a curve's edges.  It lives on
prepared rings and line strings; consumers discover it with
`getprep(curve, AbstractEdgeTree)` and read the tree through
[`edge_tree`](@ref), so any subtype with any SpatialTreeInterface tree works.
"""
abstract type AbstractEdgeTree <: AbstractPreparation end

"""
    edge_tree(p::AbstractEdgeTree)

The spatial index stored by an edge-tree preparation.  Defaults to `p.tree`.
"""
edge_tree(p::AbstractEdgeTree) = p.tree

"""
    EdgeTree(curve; backend = NaturalIndex, manifold = Planar())

A spatial index over the edge extents of a curve — the default
[`AbstractEdgeTree`](@ref), built by `prepare` for every linear ring and line
string on the planar manifold.  The index space is trait-keyed as described
in [`build_edge_tree`](@ref), which `backend` also selects the tree through:
`NaturalIndex` (default), `STRtree`, a `FlexibleRTrees` bulk-load algorithm,
or any callable `curve -> spatial tree`.
"""
struct EdgeTree{T} <: AbstractEdgeTree
    tree::T
    # Explicit inner constructor: the geometry-taking outer constructor below
    # would otherwise overwrite the default `EdgeTree(tree)` method.
    EdgeTree{T}(tree) where T = new{T}(tree)
end

function EdgeTree(geom; backend = NaturalIndex, manifold::Manifold = Planar())
    GI.trait(geom) isa GI.AbstractCurveTrait || throw(ArgumentError(
        "`EdgeTree` requires a curve (linear ring or line string), got $(typeof(GI.trait(geom)))"))
    tree = build_edge_tree(manifold, backend, geom)
    return EdgeTree{typeof(tree)}(tree)
end

# `prepare` reaches `EdgeTree` through this seam so the tree is built on the
# manifold being prepared for, not the planar default.
buildprep(m::Manifold, ::Type{EdgeTree}, geom) = EdgeTree(geom; manifold = m)

"""
    EdgeTree(backend)

The curried spec form: applying the `EdgeTree` constructor to a backend
(instead of a geometry) returns a spec for [`prepare`](@ref)'s `preps`
tuple.  Like the bare `EdgeTree` type, it applies to every curve of the
recursion (see [`appliesto`](@ref)), building
`EdgeTree(curve; backend, manifold)` there:

```julia
prep = prepare(poly; preps = (EdgeTree(STRtree),))
prep = prepare(poly; preps = (EdgeTree(FlexibleRTrees.HPR()),))
```
"""
EdgeTree(backend::Union{Base.Callable, FlexibleRTrees.BulkLoadAlgorithm}) = _EdgeTreeSpec(backend)

struct _EdgeTreeSpec{B}
    backend::B
end
buildprep(m::Manifold, s::_EdgeTreeSpec, geom) = EdgeTree(geom; backend = s.backend, manifold = m)

"""
    appliesto(spec, trait, istop)::Bool

Where a spec given in `prepare`'s `preps` tuple applies during the
recursion.  The fallback is the top node only, so an unadorned custom spec
means "this preparation, on the geometry I called `prepare` on".
`EdgeTree` (bare or curried) declares itself for every curve instead.
Nodes where no given spec applies fall back to
[`default_preparations`](@ref).
"""
appliesto(spec, trait, istop) = istop
appliesto(::Type{EdgeTree}, trait, istop) = trait isa GI.AbstractCurveTrait
appliesto(::_EdgeTreeSpec, trait, istop) = trait isa GI.AbstractCurveTrait

"""
    build_edge_tree(manifold, backend, curve)

Build a SpatialTreeInterface-compatible spatial index over the edges of
`curve`.  Leaf `i` is the edge from point `i` to point `i + 1`; an unclosed
*ring* gets one extra leaf for its implicit closing edge (last point back to
point `1`), while a line string indexes exactly its consecutive point pairs.
Consumers rely on this trait-keyed index space.  (Curves inside a `Prepared`
geometry never exercise the unclosed-ring case — materialized rings are
closed — so it only matters for trees built over raw geometry.)

Only `Planar()` methods exist so far, for `NaturalIndex`, `STRtree`, and
`FlexibleRTrees` bulk-load algorithms; the fallback calls `backend(curve)`,
so any callable works.  Spherical edge trees (arc-extent leaves over
`UnitSphericalPoint` storage) are the intended extension point.
"""
build_edge_tree(::Planar, backend, curve) = backend(curve)
build_edge_tree(::Planar, ::Type{<:NaturalIndex}, curve) = NaturalIndex(_edge_extents(curve))
build_edge_tree(::Planar, ::Type{<:STRtree}, curve) = STRtree(_edge_extents(curve))
build_edge_tree(::Planar, alg::FlexibleRTrees.BulkLoadAlgorithm, curve) =
    FlexibleRTrees.RTree(alg, _edge_extents(curve))
build_edge_tree(backend, curve) = build_edge_tree(Planar(), backend, curve)

# Extents of a curve's edges, in the trait-keyed index space described in the
# `build_edge_tree` docstring: only an unclosed *ring* gets the extra
# wrap-around edge.  Coordinate number types are preserved.
function _edge_extents(curve)
    n = GI.npoint(curve)
    wrap = GI.trait(curve) isa GI.LinearRingTrait &&
        !equals(GI.getpoint(curve, 1), GI.getpoint(curve, n))
    return [begin
        p1 = GI.getpoint(curve, i)
        p2 = GI.getpoint(curve, i == n ? 1 : i + 1)
        Extents.Extent(
            X = minmax(GI.x(p1), GI.x(p2)),
            Y = minmax(GI.y(p1), GI.y(p2)),
        )
    end for i in 1:(wrap ? n : n - 1)]
end
