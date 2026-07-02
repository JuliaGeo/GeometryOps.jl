# Prepared geometry (F4): `Prepared{Geom, Preparations}` design

**Status: ACCEPTED — implemented on branch `prepared-geometry`; deviations below.**

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
| `PointInArea()` | `Union{PolygonTrait, MultiPolygonTrait}` | `PointInAreaLike` | `IndexedPointInAreaLocator` (sorted leaves) | RelateNG per-polygon locators |

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

## Decisions (ratified)

Decisions 1–2 were ratified by Anshul up front. Decisions 3–8 were provisional in the
draft and are now resolved: recursive wrapping (3) was ratified verbally by Anshul, and
4–8 were ratified by implementation on the `prepared-geometry` branch. The per-decision
statuses below reflect that; see the Amendments section for the deviations encountered.

1. **Eager + immutable cache model** — *ratified* (chosen over lazy-inside and
   JTS-style mutable get-or-cache).
2. **Design F4 now, `PreparedRelate` kept with a sourcing seam** — *ratified*.
3. **Recursive wrapping for sub-geometry preparations** — *ratified (verbally, by Anshul)
   and implemented*.
   Alternatives considered: flat store with level tags + path bookkeeping (leaks
   `owners`-style bookkeeping into every consumer; sub-geometries can't travel with their
   preparations); top-level-only v1 (punts on the core of the idea). Flipping this
   reshapes §Core interface but not the trait/retrieval design.
4. **Every node wrapped, extent always cached** — *ratified by implementation*. Unifies
   with `_relate_cache_extents`; cost is deep (but regular) wrapper types.
5. **Spec/build split** (`PreparationSpec` + `build(spec, m, geom)` vs preparations
   doubling as their own configs) — *ratified by implementation* (built as `buildprep`,
   not `build` — see Amendments); chosen so built state and config don't share a struct.
6. **Manifold recorded on the wrapper; mismatch = miss** (silent fallback to the
   unindexed path rather than an error) — *ratified by implementation*; the relateng seam
   refines this (Amendments (e)): a matched-manifold `Prepared` is trusted as-is, a
   mismatched one is stripped and its extents rebuilt from coordinates.
7. **Placement: interface in GeometryOpsCore, concrete preparations in GeometryOps** —
   *ratified by implementation*; enables third-party providers without a GO dependency.
8. **`Base.get` overload as the retrieval verb** (vs a new exported `getprep` for
   everything) — *ratified by implementation*; `get` matches the PR #278 sketch, and the
   manifold-checked `getprep` is the algorithm-facing entry.

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

## Amendments (discovered during implementation, 2026-07-01)

The design was implemented on branch `prepared-geometry` (Tasks 1–11, all green). The
following deviations and clarifications were recorded during API fact-finding and code
review. The first group narrows the scope/API from the draft; the second records details
that surfaced in review of the individual tasks.

### Scope and API deviations (from API fact-finding)

- **`SegmentIndex` dropped from v1.** A flattened whole-geometry edge tree carries an
  `owners` table whose indices are only meaningful relative to relateng's
  `extract_segment_strings` traversal order; building it independently risks silent
  owner mismatches. RelateNG keeps building its own `PreparedEdgeIndex` inside
  `prepare(::RelateNG, a)`. Revisit when a second consumer needs flat segment trees.
  (The v1 table above still lists `SegmentIndex` as the original design; it is not built.)
- **Clipping `AutoAccelerator` seam deferred** for the same reason (clipping indexes
  flattened edge tables, not child extents).
- **`tuples` keyword deferred** — compose as `prepare(GO.tuples(geom); ...)`.
- **`build` is spelled `buildprep`**; `appliesto`/`buildprep` are unexported (providers
  extend via `import GeometryOpsCore: appliesto, buildprep`).
- **Topmost-wins attachment:** a spec is consumed at the highest matching level and not
  offered further down (so `ChildTree` on a MultiPolygon indexes polygons at the top
  node only). Multi-level trees = call `prepare` with the spec on the sub-geometry.
- **Mismatched-manifold `Prepared` inputs to relateng** take the ordinary rewrap path
  (correct extents rebuilt; preparations shadowed) rather than erroring.

### Clarifications from code review

- **(a) Points and multipoints pass through `prepare` unwrapped.** A bare point or
  multipoint is returned as-is, never wrapped in `Prepared`. The point-trait GeoInterface
  forwarding that `Prepared` carries (in `GeometryOpsCore/src/types/preparations.jl`)
  exists purely for method-table hygiene — to forward the point-trait accessors that would
  otherwise collide with the `AbstractGeometryTrait` forwarding loop, and to stay correct
  if a third party ever constructs a `Prepared` point directly — not because `prepare` ever
  produces one. (Task 2 review.)
- **(b) The legacy `_union_stored_extents` cache is dimensionally collapsed for a spherical
  GC with a bare point — but contained.** For a `Spherical` GeometryCollection that
  contains a bare point, the legacy plain-input extent path (`_union_stored_extents`)
  unions the point's manifold-blind 2D lon/lat box with the polygon's 3D unit-sphere box;
  `Extents.union` keeps only shared keys, so the cached extent collapses to a
  coordinate-mixed 2D box (verified `X = (0.604…, 12.0)` — a unit-sphere coordinate unioned
  with a raw degree value). This is **contained**: `_relate_extent` recomputes the driving
  extent in both the plain and prepared paths (`relate_geometry.jl:75`), so the `relate`
  output is unaffected; `prepare`'s `_union_prepared_extents` avoids the garbage cache
  entirely (it uses `rk_interaction_bounds`); and the Task 10 sweep pins the divergence via
  a 2D-vs-3D extent-keys tripwire. (Task 4 + Task 10 reviews.)
- **(c) `PointInArea` attaches to `Union{PolygonTrait, MultiPolygonTrait}`, not
  `PolygonTrait`** (the v1 table above is updated to match). This mirrors relateng's
  `addPolygonal` element granularity: `_extract_elements!` (`point_locator.jl:322-328`)
  keeps a whole polygonal geometry as ONE element, so a MultiPolygon input needs the
  locator attached at the MultiPolygon node (not per-polygon). GeometryCollection inputs,
  by contrast, decompose into per-polygon elements, and because the spec is offered
  top-down it flows past the GC node and attaches per-polygon. (Task 7 review
  adjudication.)
- **(d) `Geodesic` is unsupported by `prepare`.** `prepare` fails on the `Geodesic`
  manifold even with no preps, because it reaches `rk_interaction_bounds`, which has
  methods for `Planar` and `Spherical` only. This is a prepare-level limitation (not
  specific to any spec).
- **(e) Mismatched-manifold handling is part of the relateng seam contract.** When a
  `Prepared` input's manifold matches the algorithm's, the wrapper (and its cached extents)
  is trusted as-is; when it mismatches, the wrapper is **stripped** and extents are rebuilt
  from coordinates (commit `130e278cf`). Planar→spherical was already safe because the
  spherical path recomputes, but spherical→planar leaked cached 3D extents into the planar
  `GI.extent` fast path until the strip fix. (Task 8 review.)
- **(f) relateng `PointInArea` reuse requires matching `exact` types.** The relate point
  locator reuses a `Prepared` `PointInAreaLike` locator only when the algorithm's `exact`
  type matches the preparation's; both default to `True()`, so the common case reuses.
  (Task 9 review.)
