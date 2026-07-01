# Prepared geometry (F4): `Prepared{Geom, Preparations}` design

**Status: DRAFT — decisions 1–2 ratified by Anshul, 3–8 provisional (made in his absence, listed in §Decisions for ratification).**

This is the design for follow-up **F4** from `2026-06-10-relateng-design.md` ("Generalize
`prepare` into a package-wide prepared-geometry mechanism"). It consolidates the ideas
recorded in [PR #278](https://github.com/JuliaGeo/GeometryOps.jl/pull/278),
[issue #87](https://github.com/JuliaGeo/GeometryOps.jl/issues/87), and
[issue #103](https://github.com/JuliaGeo/GeometryOps.jl/issues/103), and uses the RelateNG
port's index inventory as the requirements list.

## Problem

GeometryOps is accreting one-off preparation state with no unifying layer:

- `PreparedRelate` (`src/methods/geom_relations/relateng/relate_ng.jl:596-632`) — an
  algorithm-side struct: RelateNG + `RelateGeometry` + segment strings + edge tree.
- `NaturallyIndexedRing` / `prepare_naturally` (`src/utils/NaturalIndexing.jl:199-240`) —
  a ring wrapper carrying a `NaturalIndex`, explicitly marked as a throwaway experiment.
- The extent-cached GI wrapper tree `_relate_cache_extents`
  (`relate_geometry.jl:74,119-139`) — rebuilt privately inside every `RelateGeometry`.
- Clipping's `AutoAccelerator` docstring promises it "will also consider the existing
  preparations on the geoms" — with no mechanism to find them.

Each new accelerated algorithm reinvents caching, caches cannot be shared across
algorithms, and a geometry cannot carry more than one preparation.

## Design summary

A **GeoInterface-forwarding wrapper** that carries a **tuple of preparations**, retrieved
**by capability trait** ("I only care that you *have* an edge index, not how you got it" —
PR #278). Preparations are built **eagerly** at `prepare` time; the wrapper and its
preparation set are **immutable**. Sub-geometry preparations are handled by **recursive
wrapping**: `prepare` rebuilds the GI tree once, and a preparation always describes the
geometry it directly wraps.

```julia
p = prepare(multipoly; preps = (ChildTree(), RingEdgeIndex()), manifold = Planar())
# Prepared{MultiPolygonTrait, …, Tuple{ChildTree{…}}}        (tree over polygon extents)
#  └─ Prepared{PolygonTrait, …, Tuple{}}                     (extent-cached pass-through)
#      └─ Prepared{LinearRingTrait, …, Tuple{RingEdgeIndex{…}}}

tree = get(p, SpatialIndexLike())        # ⇒ the ChildTree
get(p, PointInAreaLike())                # ⇒ nothing — caller uses the unindexed path
GI.getgeom(p, 1)                         # ⇒ the *prepared* polygon; GI code sees a normal geometry
```

## Core interface (GeometryOpsCore)

```julia
"A built, immutable acceleration structure attached to exactly one geometry."
abstract type AbstractPreparation end

"Capability traits: what a preparation can do, not what it is."
abstract type AbstractPreparationTrait end
struct SpatialIndexLike     <: AbstractPreparationTrait end  # extent tree over child elements
struct SpatialEdgeIndexLike <: AbstractPreparationTrait end  # extent tree over segments/edges
struct PointInAreaLike      <: AbstractPreparationTrait end  # point-in-polygonal-area locator

"Every concrete preparation declares its capability."
function preptrait end   # preptrait(::RingEdgeIndex) = SpatialEdgeIndexLike()

"Configs passed to `prepare`; declare where they attach and how to build."
abstract type PreparationSpec end
function appliesto end   # appliesto(::RingEdgeIndex) = GI.LinearRingTrait
function build end       # build(spec, manifold, geom) -> AbstractPreparation
```

The wrapper:

```julia
struct Prepared{T <: GI.AbstractTrait, M <: Manifold, Geom,
                Preps <: Tuple{Vararg{AbstractPreparation}}, E}
    geom::Geom      # children are themselves Prepared (recursive rebuild)
    manifold::M     # preparations are manifold-specific (kernel point types, 3D extents)
    preps::Preps
    extent::E       # always cached — subsumes `_relate_cache_extents`
end

GI.trait(::Prepared{T}) where {T} = T()
GI.extent(p::Prepared) = p.extent
# + full GeoInterface forwarding (ngeom/getgeom/npoint/getpoint/coordinates/crs/is3d/ismeasured),
#   following the 2024 `src/prepared.jl` (commit 225997c65) and `_relate_cache_extents`.
```

Retrieval:

```julia
# Hit: first preparation in `preps` whose preptrait matches (compile-time tuple recursion).
# Miss: `nothing`. Plain geometries always miss, so call sites are uniform.
Base.get(p::Prepared, t::AbstractPreparationTrait) -> Union{AbstractPreparation, Nothing}
Base.get(geom, ::AbstractPreparationTrait) = nothing

# Manifold-checked variant for algorithms (a Planar edge tree is useless to a Spherical kernel):
getprep(m::Manifold, p, t) = (q = get(p, t); q !== nothing && manifold(p) == m ? q : nothing)
```

Construction is one bottom-up `apply`-style pass — children are prepared first, so a
parent's preparation (e.g. a tree over ring extents) can index already-prepared children:

```julia
prepare(geom; preps = (), manifold = Planar(), tuples = true)
```

Every node is wrapped (empty `preps` allowed) so extents are cached at every level,
exactly like `_relate_cache_extents` today. `tuples = true` converts coordinates to tuple
form during the rebuild (the PR #278 default; `tuples = false` preserves the original
point representation).

## v1 preparations (GeometryOps)

| Spec | Attaches to (`appliesto`) | Capability | Built structure | Replaces / serves |
|---|---|---|---|---|
| `SegmentIndex()` | any lineal/polygonal geometry (top level) | `SpatialEdgeIndexLike` | flattened `NaturalIndex` over all segments + `owners` table | `PreparedEdgeIndex` reuse in RelateNG |
| `RingEdgeIndex()` | `LinearRingTrait` | `SpatialEdgeIndexLike` | `NaturalIndex` over the ring's segment extents | `NaturallyIndexedRing` (delete it) |
| `ChildTree()` | `Multi*`, `GeometryCollection`, `PolygonTrait` (rings) | `SpatialIndexLike` | STRtree or `NaturalIndex` over child extents | clipping `AutoAccelerator` promise |
| `PointInArea()` | `PolygonTrait` | `PointInAreaLike` | `IndexedPointInAreaLocator` (sorted leaves) | RelateNG per-polygon locators |

Notes:
- All specs take the index implementation as a config (default `NaturalIndex`; anything
  implementing SpatialTreeInterface works — the swap point documented at
  `edge_intersector.jl:294-299`).
- `PointInArea` is Planar-only in v1; on `Spherical`, `build` is not defined and `prepare`
  raises a clear error (the spherical point-in-area index was deliberately time-boxed out
  of the relateng plan — Task 16).
- Third-party providers (e.g. a LibGEOS prepared geometry) can ship their own
  `AbstractPreparation` + `preptrait` method — the issue #103 goal. Nothing in the core
  interface names a concrete index type.

## Consumer seams (v1 validation)

1. **RelateNG** (`PreparedRelate` stays; it is algorithm-side state, not geometry state):
   - `RelateGeometry` constructor: if the input is `Prepared`, skip
     `_relate_cache_extents` — extents are already on the wrapper.
   - `_build_prepared_edge_index`: `getprep(m, a, SpatialEdgeIndexLike())` hit ⇒ reuse
     instead of rebuilding.
   - `RelatePointLocator`: a polygonal element that is `Prepared` with `PointInAreaLike`
     ⇒ use it as the indexed locator for that element.
2. **`NaturallyIndexedRing`** is deleted; its two call sites move to
   `prepare(ring; preps = (RingEdgeIndex(),))`.
3. **Clipping `AutoAccelerator`**: when an input `isa Prepared` with `SpatialIndexLike`,
   prefer the existing tree over building one.

Success criterion: the relateng prepared-vs-unprepared equality suite
(`test/methods/relateng/relate_ng.jl`, `PREPARED_FIXTURES`) passes with `Prepared` inputs
substituted, plus GI-equivalence tests (`GI.coordinates`, `GI.extent`, trait round-trips)
for wrapped vs plain geometries.

## What this design deliberately does not do

- **No laziness, no mutation.** `get` never builds; a miss means the caller takes the
  unindexed path. Laziness is the property of the *unprepared* regime (e.g.
  `RelatePointLocator`'s second-query indexing stays as-is). Consequence: a
  `PointInArea()` prep builds locators for *all* polygonal elements up front.
- **No survival across transformations.** `apply`/`transform`/corrections return plain
  geometries; preparations are dropped, never patched. Documented, not enforced.
- **No Symbol-keyed access.** The 2024 `getprep(p, ::Symbol)` NamedTuple approach is
  superseded by capability traits.
- **No serialization story.**

## Decisions for ratification

1. **Eager + immutable cache model** — *ratified* (chosen over lazy-inside and
   JTS-style mutable get-or-cache).
2. **Design F4 now, `PreparedRelate` kept with a sourcing seam** — *ratified*.
3. **Recursive wrapping for sub-geometry preparations** — *provisional (recommended)*.
   Alternatives considered: flat store with level tags + path bookkeeping (leaks
   `owners`-style bookkeeping into every consumer; sub-geometries can't travel with their
   preparations); top-level-only v1 (punts on the core of the idea). Flipping this
   reshapes §Core interface but not the trait/retrieval design.
4. **Every node wrapped, extent always cached** — provisional. Unifies with
   `_relate_cache_extents`; cost is deep (but regular) wrapper types.
5. **Spec/build split** (`PreparationSpec` + `build(spec, m, geom)` vs preparations
   doubling as their own configs) — provisional; chosen so built state and config don't
   share a struct.
6. **Manifold recorded on the wrapper; mismatch = miss** (silent fallback to the
   unindexed path rather than an error) — provisional.
7. **Placement: interface in GeometryOpsCore, concrete preparations in GeometryOps** —
   provisional; enables third-party providers without a GO dependency.
8. **`Base.get` overload as the retrieval verb** (vs a new exported `getprep` for
   everything) — provisional; `get` matches the PR #278 sketch, and the manifold-checked
   `getprep` is the algorithm-facing entry.

## References

- PR #278 (`as/prepare` branch): `Prepared{Pa,Pr}`, `AbstractPreparation`,
  `AbstractPreparationTrait`, `preptrait`, `SpatialIndex{T}`, `SpatialEdgeIndex{T}`.
- Issues #87 (closed → #278/#103), #103 (open).
- Precedents in-tree: `_relate_cache_extents` (`relate_geometry.jl:119-139`),
  `NaturallyIndexedRing` (`src/utils/NaturalIndexing.jl:199-240`), commit `225997c65`
  (rafaqz's `src/prepared.jl`).
- Requirements inventory: `PreparedRelate` internals (`relate_ng.jl:596-676`),
  `RelatePointLocator` (`point_locator.jl:235-270`), `AutoAccelerator`
  (`clipping_processor.jl:25-40`).
