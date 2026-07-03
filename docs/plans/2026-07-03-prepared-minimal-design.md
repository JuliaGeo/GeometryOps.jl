# Minimal prepared geometry — design

**Status:** in progress (2026-07-03)

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

## Benchmark plan (profile-performance skill)

Regular n-gon ("circle") polygons with n ∈ 2⁶…2¹⁶ vertices; point queries in /
out / near-boundary. Compare: plain `contains`, prepared with `NaturalIndex`,
prepared with `STRtree`; plus `prepare` construction cost and per-query
allocations. Report the crossover size where preparation wins and the asymptotic
speedup.
