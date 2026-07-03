# Minimal prepared geometry — design

**Status:** implemented & benchmarked (2026-07-03) — see Results at the bottom

A clean, minimal reimplementation of prepared geometry, off `main`, entirely inside
GeometryOps (no GeometryOpsCore changes). Supersedes the earlier
`prepared-geometry` branch design for the purpose of understanding structure and
benchmarking; deliberately drops manifold parameters, capability traits, spec
objects, and recursive child wrapping.

## Goals

1. **Fundamentals only**: a wrapper that carries preparations alongside a geometry,
   a way to build it, and a way to query it — nothing else.
2. **One real consumer**: edge-index trees wired into the planar point-in-polygon
   test via `SpatialTreeInterface.depth_first_search` with a custom ray-strip
   predicate. This is the profiling target.
3. **Extensible at every joint**: users can add new preparation types, new tree
   backends, and new consumers without touching GeometryOps internals.
4. **Multiple instantiations per preparation kind**: e.g. edge trees may be backed
   by `NaturalIndex` *or* `STRtree` (or a user tree) — consumers stay generic via
   SpatialTreeInterface.

## Types

```julia
abstract type AbstractPreparation end

struct Prepared{T <: GI.AbstractGeometryTrait, G, P <: Tuple, E}
    geom::G      # the wrapped geometry, unchanged
    preps::P     # tuple of AbstractPreparation instances (heterogeneous, type-searched)
    extent::E    # cached extent — every Prepared has one
end
```

`Prepared` forwards the GeoInterface of its parent (`geomtrait`, `ngeom`,
`getgeom`, `ncoord`, `getcoord`, `npoint`, `is3d`, `ismeasured`, `crs`), except
`GI.extent`, which returns the cached extent. Children come back **plain** — no
recursive wrapping. A preparation that needs per-child structure holds it
internally.

## API

```julia
prepare(geom; preps = default_preparations(GI.trait(geom), geom)) -> Prepared
prepare(p::Prepared; preps = ())   # add preparations to an existing Prepared (prepended, so they win lookup)

getprep(geom, P::Type)             # -> the first prep `isa P`, or `nothing`. Plain geometries always give `nothing`.
getprep(f, geom, P::Type)          # get-or-else: the prep if present, otherwise `f()`
```

The three usage idioms, in consumer code:

```julia
# 1. Choose-path: use the index if it exists, otherwise the plain algorithm.
#    Right when building would cost as much as one unindexed query (e.g. point-in-polygon).
tree = getprep(poly, AbstractRingEdgeTrees)
isnothing(tree) ? plain_path(poly) : indexed_path(poly, tree)

# 2. Get-or-create: build ephemerally on a miss.
#    Right when the index pays for itself within one operation (e.g. clipping).
trees = getprep(poly, AbstractRingEdgeTrees) do
    RingEdgeTrees(poly)
end

# 3. Query-only: hand the prep to something else.
idx = getprep(poly, RingEdgeTrees)   # concrete type also works
```

### Extensibility hooks

| Hook | Default | Overload to… |
|---|---|---|
| `buildprep(spec, geom)` | `spec(geom)` | make a custom spec object buildable; types and closures already work |
| `default_preparations(trait, geom)` | `()`; `(RingEdgeTrees,)` for polygons | change what `prepare` builds by default (method on *your* trait/type) |
| `build_edge_tree(backend, ring)` | `backend(ring)`; special-cased for `NaturalIndex` & `STRtree` | add a new tree backend |
| `exterior_tree(p)` / `hole_trees(p)` | field access | supply your own `AbstractRingEdgeTrees` subtype with different storage |

A user preparation is: `struct MyPrep <: GO.AbstractPreparation ... end` plus a
constructor `MyPrep(geom)`. It then flows through `prepare(geom; preps = (MyPrep,))`
and `getprep(geom, MyPrep)` with no further registration. Retrieval matches by
`isa`, so querying by an abstract kind (e.g. `AbstractRingEdgeTrees`) finds any
subtype — that is the whole "capability" mechanism.

## v1 preparation: `RingEdgeTrees`

```julia
abstract type AbstractRingEdgeTrees <: AbstractPreparation end

struct RingEdgeTrees{T} <: AbstractRingEdgeTrees   # T = tree type, e.g. NaturalIndex{…} or STRtree
    exterior::T
    holes::Vector{T}
end

RingEdgeTrees(polygon; tree = NaturalIndex)   # tree = STRtree, or any callable ring -> spatial tree
```

One spatial tree per ring, over **edge extents**. Edge `i` of a ring with `n`
points is `(point i, point i+1)`, except when the ring is unclosed, where edge
`n` wraps to `(point n, point 1)` — so prepared results match the plain
algorithm's implicit-closure semantics exactly.

## The seam: indexed planar point-in-polygon

`_point_filled_curve_orientation` (Hao & Sun 2018) is per-edge independent: every
edge either proves the point **on** the boundary, increments the ray-crossing
counter, or contributes nothing. Edges whose bounding box misses the rightward
ray strip `{(x′, y′) : x′ ≥ x, y′ = y}` can do none of these, so:

1. The 26-case edge test is extracted into `_hao_sun_edge(x, y, p_start, p_end; exact)`
   — shared verbatim by the plain loop and the indexed path (one source of truth).
2. A tree-accepting method runs `depth_first_search` with the custom predicate
   `ext -> ext.X[2] ≥ x && ext.Y[1] ≤ y ≤ ext.Y[2]`, accumulates crossings, and
   early-exits the whole traversal on an *on* hit via `Action(:full_return, on)`.
3. `_point_polygon_process` looks up `getprep(polygon, AbstractRingEdgeTrees)`
   once and passes the per-ring tree (or `nothing`) down. This accelerates
   `contains`, `within`, `covers`, `coveredby`, `intersects`, `disjoint`, and
   `touches` for point-vs-polygon whenever the polygon is prepared.

The cached extent also makes the existing `_maybe_skip_disjoint_extents`
short-circuit free for prepared geometries (plain geometries recompute their
extent per call).

## Non-goals for v1 (deliberate)

- No manifold parameter on `Prepared` — the indexed seam is planar-only, matching
  main's predicates. A spherical prep can be added later as its own type.
- No recursive child wrapping, no `MultiPolygon` PIP index (prepare(multipolygon)
  still caches the extent; a flat all-rings tree à la JTS is the natural v2).
- No line/polygon–polygon edge-pair acceleration (that loop carries sequential
  state; RelateNG territory).
- No clipping rewiring — but `getprep`'s get-or-create form is shaped for it.

## Results (2026-07-03, Julia 1.12.6, Apple Silicon)

`GO.contains(poly, pt)` median time per query, 1000 uniform points over the
bounding box of a regular n-gon (~78% inside):

| n vertices | plain | prepared (NaturalIndex) | prepared (STRtree) | speedup (Nat) |
|---:|---:|---:|---:|---:|
| 64     | 0.77 µs | 73 ns  | 44 ns  | 11× |
| 256    | 3.3 µs  | 94 ns  | 49 ns  | 36× |
| 1 024  | 13.6 µs | 165 ns | 131 ns | 82× |
| 4 096  | 54.3 µs | 240 ns | 209 ns | 226× |
| 16 384 | 218 µs  | 236 ns | 363 ns | 925× |
| 65 536 | 873 µs  | 300 ns | 779 ns | 2 916× |

- **There is no crossover size** — prepared wins from n = 64 up. The plain path
  is O(n) per query even for trivially-rejectable points, because
  `_maybe_skip_disjoint_extents` recomputes the polygon extent per call; a
  far-outside point costs 808 µs plain vs **10.6 ns** prepared at n = 65 536.
- **Build cost** (`prepare`, NaturalIndex): ≈ one plain query — 16 µs / 34 KB at
  n = 1 024, 0.93 ms / 2.2 MB at n = 65 536. STRtree builds ~7× slower and ~6×
  bigger; it queries faster below ~8 k vertices, NaturalIndex faster above.
  `prepare` pays for itself by the second query.
- **Profile** (n = 65 536, NaturalIndex): zero runtime dispatch in the hot path
  (`@inferred` clean); self-cost is tree traversal (`depth_first_search`,
  `NaturalIndexNode` materialization, generator iteration). The only per-query
  allocation is the 16-byte `Ref` crossing counter (1/query).

Possible follow-up optimizations (not applied): a fold/reduce variant in
SpatialTreeInterface to eliminate the `Ref`; a level-iterating query
specialized to `NaturalIndex` to shave traversal overhead.
