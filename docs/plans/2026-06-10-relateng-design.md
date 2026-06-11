# RelateNG for GeometryOps.jl — Design

**Date:** 2026-06-10
**Status:** Validated design, pre-implementation
**Reference implementation:** JTS `org.locationtech.jts.operation.relateng` (Martin Davis, JTS 1.20)

## Context and motivation

GeometryOps' DE-9IM predicates (`intersects`, `within`, `touches`, …) are implemented as
~2,500 lines of hand-written per-trait-pair processors in `src/methods/geom_relations/`.
They are planar-only, have no `relate()` returning a full DE-9IM matrix, no prepared
(build-once, query-many) mode, and known gaps in degenerate cases.

JTS's RelateNG engine is a single unified evaluator for all DE-9IM predicates with
aggressive short-circuiting, prepared-geometry support, pluggable boundary node rules,
and well-tested handling of dark corners (GeometryCollections, empty/zero-length inputs,
invalid geometry tolerance). This document specifies a Julia implementation of RelateNG
in GeometryOps, designed for:

1. **Exactness** — answers computed via exact geometric predicates, with *no constructed
   (rounded) intersection coordinates* in the topology logic, unlike JTS.
2. **Manifold generality** — the engine is generic over `Manifold`; `Planar` ships first,
   `Spherical` is a planned follow-up implementing the same kernel contract.
3. **Performance** — predicate-specialized compilation, allocation-free hot paths,
   reuse of GeometryOps' spatial-index acceleration.

## Key decisions

### D1. Split at the layer boundary: faithful topology layer, redesigned geometry layer

RelateNG is two things glued together, and we treat them differently:

- **Topology layer** (predicate contract, IM bookkeeping, evaluation phases,
  node-labelling semantics, GC/empty/boundary-rule handling): purely combinatorial and
  manifold-independent. **Ported faithfully** — same phase structure, same predicate
  requirement flags, same labelling rules, file-by-file traceable to the Java source,
  validated by the JTS XML test suite. This is where JTS's hard-won correctness lives;
  we keep it diffable against the reference and able to absorb upstream fixes.
- **Geometry layer** (every question asked about coordinates): **redesigned** as a small
  manifold-parameterized exact kernel with symbolic node identity. This is exactly where
  JTS's floating-point construction and planar assumptions live.

The interface between the layers is the `RelateKernel` API (below) — it is also the
contract a future `Spherical` kernel implements.

### D2. Exact predicates; no constructed intersection points

JTS constructs intersection coordinates in floating point
(`EdgeSegmentIntersector.java:71`, via `RobustLineIntersector`) and makes them
load-bearing in two places: node identity (`TopologyComputer` groups `NodeSection`s in a
`Map<Coordinate, NodeSections>`) and angle ordering (the constructed point is the apex of
`PolygonNodeTopology.compareAngle`). JTS's own comments acknowledge the containing-segment
check "is not reliable" near proper intersections due to roundoff.

We eliminate constructed points. Every segment-pair intersection is first classified
*combinatorially* by exact orientation/on-segment tests:

- **At a vertex** of A and/or B → the node is that input vertex; its coordinate is exact.
- **Proper interior crossing** → the node is represented *symbolically* as the segment
  pair. Its DE-9IM contribution (interior∩interior) and the cyclic order of the four
  incident half-edges are decidable by orientation tests on the original endpoints; the
  coordinate is never needed.

This is possible because relate produces no output geometry: every answer is a finite
set of sign computations on input coordinates.

> **Amendment (2026-06-11, Task 18).** One bounded deviation: locating a *crossing
> node* against a **lineal or GeometryCollection target** uses a representative
> Float64 point, correctly rounded from the exact rational crossing point
> (`_crossing_locate_point` in `topology_computer.jl`). Polygonal targets still need
> no coordinate (a node of a polygonal geometry is on its boundary, exactly).
> The failure window of the rounded point is only when the exact crossing lies
> within half an ULP of another element's endpoint (lineal targets — and only a
> false BOUNDARY answer is possible there; a false INTERIOR is impossible, since a
> crossing that is exactly representable rounds to itself) or within half an ULP of
> a non-parent polygon's boundary (GC targets). Correct rounding means this matches
> or beats JTS's `RobustLineIntersector`-constructed coordinate in every such
> window. A fully exact alternative (rational on-segment / rational point-in-polygon
> tests over the target's edges) is recorded as **Follow-up F6**.

### D3. Exact crossing-node coincidence via a slow path (follow-up: make it fast)

Whether two *proper crossings* coincide (only relevant for self-intersecting/invalid
input or predicates with `require_self_noding`) is decided by exact comparison of
rational intersection parameters (extended-precision/rational arithmetic). This is
accepted as a slow path for now. **Follow-up F1** investigates a fast filter
(interval arithmetic before the rational fallback) or a proof that the case is
unreachable for valid inputs with non-self-noding predicates.

### D4. Integration: `Algorithm` type, opt-in first

The engine is exposed as `RelateNG{M<:Manifold} <: Algorithm{M}` (mirroring `GEOS()`,
`TG()`, `FosterHormannClipping`). The new `relate(...)` always uses it. Existing named
predicates gain `intersects(RelateNG(), a, b)` etc. methods while their current
implementations remain the default. The default flips in a later release once the JTS
XML suite passes (**Follow-up F3**) — mirroring JTS's own `-Djts.relate=ng` migration.

## Public API

```julia
relate(a, b)::DE9IM                    # full DE-9IM matrix
relate(a, b, pattern::String)::Bool    # pattern match ("T*F**FFF*", with 0/1/2/T/F/*)
relate(alg::RelateNG, a, b, ...)       # explicit algorithm form

RelateNG(; boundary_rule = Mod2Boundary(), accelerator = AutoAccelerator(), exact = True())
RelateNG(m::Manifold; ...)             # Planar() default

prepared = prepare(RelateNG(), a)      # PreparedRelate: cached locators + edge index
relate(prepared, b); intersects(prepared, b); ...

# Named predicates, opt-in:
intersects(RelateNG(), a, b), within(RelateNG(), a, b), contains, covers,
coveredby, touches, crosses, overlaps, disjoint, equals
```

`DE9IM` is an immutable struct wrapping the 9 entries (dimension codes `F/0/1/2` packed
compactly, e.g. 2 bits × 9 in a `UInt32` or `NTuple{9,Int8}`), with `Base.show` printing
the standard string form (`"212101212"`), `matches(im, pattern)`, and indexed access
`im[Interior, Boundary]`.

Boundary node rules are zero-size structs: `Mod2Boundary` (OGC default),
`EndpointBoundary`, `MultivalentEndpointBoundary`, `MonovalentEndpointBoundary`.

## Code layout

```
src/methods/geom_relations/relateng/
├── kernel.jl                  # RelateKernel API definition (the layer contract)
├── kernel_planar.jl           # Planar implementation
├── de9im.jl                   # DE9IM type, Location/Dimension codes, DimensionLocation packing
├── topology_predicate.jl      # predicate framework: TopologyPredicate, BasicPredicate, IMPredicate
├── relate_predicates.jl       # named predicates + IMPatternMatcher + RelateMatrixPredicate
├── relate_geometry.jl         # RelateGeometry: input facade (dims, bounds, points, edges)
├── point_locator.jl           # RelatePointLocator, LinearBoundary, AdjacentEdgeLocator
├── node_sections.jl           # NodeSection, NodeSections, symbolic NodeId
├── polygon_node_converter.jl  # minimal→maximal ring rewriting at nodes
├── relate_node.jl             # RelateNode, RelateEdge: edge wheel, label propagation
├── topology_computer.jl       # TopologyComputer: IM accumulation, early exit
├── edge_intersector.jl        # edge enumeration via SpatialTreeInterface + classification
└── relate_ng.jl               # RelateNG algorithm type, evaluation phases, prepare
```

One file per ported JTS concept on the topology layer, so each diffs against its Java
counterpart.

## The `RelateKernel` API (geometry layer)

Every coordinate-level question the engine may ask. Each function takes the manifold and
GeometryOps' `exact::True/False` flag, and returns **only discrete classifications —
never constructed coordinates**:

| Function | Returns | Planar implementation |
|---|---|---|
| `orient(m, exact, a, b, c)` | `-1/0/+1` | `Predicates.orient` (AdaptivePredicates); spherical later: sign of scalar triple product, `RobustCrossProduct` robust path |
| `classify_intersection(m, exact, a0, a1, b0, b1)` | symbolic: `disjoint`, `proper_cross`, `touch(endpoint incidences)`, `collinear_overlap(endpoint ordering)` | orientation signs + on-segment tests; replaces `RobustLineIntersector` |
| `point_on_segment(m, exact, p, q0, q1)` | `Bool` (and endpoint/interior distinction) | orient == 0 + range check |
| `point_in_ring(m, exact, p, ring)` | in/on/out | existing Hao–Sun machinery; spherical later: arc-crossing count from a reference exterior point |
| `edge_cmp_around_node(m, exact, node, d1, d2)` | cyclic ordering | JTS `compareAngle` logic, quadrant + orient, but apex is a *symbolic* node |
| `node_id(...)` / `nodes_coincide(...)` | node key / `Bool` | `VertexNode(coord)` keys exactly by coordinate; `CrossingNode(segA, segB)` compares by exact rational intersection parameters (D3 slow path) |
| `interaction_bounds(m, geom)`, `bounds_disjoint`, `bounds_covers` | bounds / `Bool` | `Extents`; spherical later: lon-wrapped extents or caps |

The topology layer may **only** call these. A kernel-conformance testset is written
against this API and instantiated for `Planar`; the future `Spherical` kernel must pass
the same suite on great-circle analogues.

## Topology layer (faithful port)

JTS-class → Julia mapping, with idiom changes:

- **Predicates as concrete structs.** Each named predicate is a small mutable struct
  holding its tri-state value / partial IM. JTS interface methods become functions:
  `init_dims!`, `init_bounds!`, `update_dim!`, `finish!`, `is_known`. The declarative
  flags (`require_interaction`, `require_covers`, `require_exterior_check`,
  `require_self_noding`) are pure functions of the predicate *type*, so the whole
  evaluation specializes per predicate and unused checks are dead-code-eliminated —
  a performance win unavailable to JVM interface dispatch.
- **`RelateGeometry`**: real dimension (zero-length-line demotion to P), emptiness,
  bounds, lazily-built unique point set, line ends, lazily-built `RelatePointLocator`.
  GeometryCollection *union semantics* ported exactly (highest dimension wins;
  `AdjacentEdgeLocator` resolves boundary-vs-interior for points on multiple polygon
  boundaries).
- **`RelatePointLocator`**: location + dimension packed as `DimensionLocation` codes
  (`Int8`); `LinearBoundary` valence counting parameterized by boundary rule.
- **`TopologyComputer`**: same entry points (`add_point_on_geometry!`,
  `add_line_end!`, `add_area_vertex!`, `add_intersection!`), node-sections dictionary
  keyed by symbolic `NodeId`, same exterior-dimension initialization, early-exit checks
  after every update.
- **Node analysis ported verbatim in logic**: `NodeSections` grouping →
  `PolygonNodeConverter` (minimal→maximal ring rewriting when shell/holes meet at a
  node) → `RelateNode`/`RelateEdge` building the CCW edge wheel, merging collinear
  edges (area-over-line override), propagating left/right locations around the wheel —
  with every geometric comparison routed through the kernel and the angle-ordering apex
  being the symbolic node.
- **Evaluation phases** match `RelateNG.evaluate` (RelateNG.java:222–268) exactly:
  1. required-bounds interaction screen (`require_covers` / `require_interaction`)
  2. `init_dims!` — exit if known (incl. empty-geometry handling)
  3. `init_bounds!` — exit if known
  4. point–point fast path when both inputs are puntal
  5. B points located against A; A line-ends/area-vertices located against B; exit
     checks throughout
  6. edge phase: mutual edge intersection enumeration + node topology analysis
  7. `finish!` — fill remaining IM entries from dimension defaults; return value

## Edge enumeration, indexing, prepared mode

- **No monotone-chain port.** Edge sets from A and B each get an extent tree (STRtree or
  `NaturalIndex` behind `SpatialTreeInterface`); the dual depth-first traversal yields
  candidate segment pairs, each handed to `classify_intersection`. This mirrors the
  clipping `IntersectionAccelerator` pattern (`AutoAccelerator` picks nested-loop below
  a size threshold), is user-selectable via a `RelateNG` field, and works unchanged on
  the sphere with cap/wrapped extents — avoiding monotone chains' planar x-sortedness.
- **Prepared mode.** `prepare(RelateNG(), a)` builds once: A's edge tree, indexed
  point-in-area locators per polygon, unique-point set, `LinearBoundary`. Each
  evaluation against a B only indexes B's side (or nested-loops if B is small).
  Same caveat as JTS: predicates requiring self-noding bypass parts of the cache.
  This should seed GeometryOps' broader prepared-geometry story (**Follow-up F4**).

## Performance disciplines

- Points are coordinate tuples; `NodeSection` immutable; node dictionary maps
  `NodeId → Vector{NodeSection}`.
- Hot path (per-discovery predicate update) non-allocating: IM packed in `UInt32`,
  locations/dimensions as `Int8` codes.
- Whole-evaluation specialization on `{Manifold, predicate type, exact}`; verified with
  `@allocated`/JET-adjacent tests.
- Benchmarks under `benchmarks/`: RelateNG vs existing GO predicates vs LibGEOS
  (incl. prepared GEOS), across polygon-scaling providers, for both early-exit-friendly
  (`intersects`) and full-matrix (`relate`) workloads.

## Testing and validation

**Extend the existing harness** at `test/external/jts/jts_testset_reader.jl` (XML.jl +
WellKnownGeometry), which is currently overlay-specific:

- Parse `expected_result` by operation kind: geometry (overlay, as now), **boolean**
  (named predicates), **DE-9IM string** (`relate`, pattern in `arg3`). Read
  `<precisionModel>` headers.
- Vendor relate-relevant XML files (`TestRelate*.xml`, `TestFunctionPL/LL/LA/AA*.xml`,
  boundary-node-rule files; EPL/EDL-licensed) into `test/data/jts/` so CI has no JTS
  checkout dependency. Replace the hardcoded path.
- Wire a `relate` dispatcher into the case runner; register in `runtests.jl`.

Cross-validation layers:

1. **Vendored JTS XML suite** = planar ground truth. Intentional divergences go on an
   explicit, documented skip-list — never silent.
2. Existing 63-pair predicate suite (`test/methods/geom_relations.jl`) run through
   `RelateNG()` via `@testset_implementations`; must agree with current implementations
   and LibGEOS (GEOS ≥ 3.13 itself runs RelateNG — an independent check).
3. **Differential fuzzing** vs LibGEOS: random pairs from `benchmarks/geometry_providers.jl`
   plus adversarial near-degenerate cases (shared vertices, collinear edges,
   ulp-perturbed crossings), with `exact=True()`. Exact-vs-floating divergences are
   triaged and become documented exactness wins, not failures.
4. **Per-file unit tests** ported from JTS JUnit suites (`RelateNGTest`,
   `RelatePointLocatorTest`, `PolygonNodeConverterTest`, `LinearBoundaryTest`, …),
   written test-first alongside each stage.
5. **Kernel conformance testset** against the abstract kernel API (the spec for the
   future spherical kernel).

## Phasing

Each stage lands reviewable and green:

1. **Foundations** — `DE9IM`, location/dimension codes, boundary-rule structs;
   predicate framework + all named predicates + pattern matcher, with ported unit
   tests. Purely combinatorial.
2. **Kernel** — `RelateKernel` API + `Planar` implementation + conformance testset
   (incl. the rational crossing-coincidence slow path).
3. **Point location** — `RelateGeometry`, `RelatePointLocator`, `LinearBoundary`,
   `AdjacentEdgeLocator`; XML-harness generalization lands here so point-only predicate
   paths validate immediately.
4. **Node topology** — `NodeSection(s)`, `PolygonNodeConverter`, `RelateNode`,
   `RelateEdge`, `TopologyComputer`.
5. **Engine** — `RelateNG` phases, edge enumeration via dual tree traversal, prepared
   mode; full XML suite + differential fuzzing pass. (GeometryCollection support may
   split into its own stage if this grows large.)
6. **Surface & perf** — `relate()` export, predicate methods on `RelateNG()`, literate
   docs, benchmarks, allocation checks.

## Follow-up register (non-blocking)

- **F1.** Faster exact crossing-node coincidence: interval-arithmetic filter before the
  rational fallback, or prove unreachability for valid inputs + non-self-noding
  predicates.
- **F2.** `Spherical` kernel: `RobustCrossProduct` orientation, great-circle-arc
  classification, cap/wrapped bounds; must pass the kernel conformance suite.
- **F3.** Flip named-predicate defaults to RelateNG once parity is proven; deprecate the
  per-pair processors in `geom_relations/`.
- **F4.** Generalize `prepare` into a package-wide prepared-geometry mechanism.
- **F5.** GeometryCollection edge cases beyond the JTS XML suite's coverage (mixed-dim
  GCs with overlapping components) — extra fuzzing.
- **F6.** Exact crossing-node location against lineal/GC targets: replace the
  correctly-rounded representative point of `_crossing_locate_point` with rational
  on-segment / rational point-in-polygon tests over the target's edges (see the
  D2 amendment of 2026-06-11 for the half-ULP failure window this would close).

## References

- JTS sources: `jts/modules/core/src/main/java/org/locationtech/jts/operation/relateng/`
- JTS XML tests: `jts/modules/tests/src/test/resources/testxml/`
- Martin Davis, "JTS Topological Relationships — the Next Generation" (lin-ear-th-inking
  blog, 2024) and the RelateNG design notes in the JTS repository (`doc/`).
- GEOS 3.13+ `RelateNG` port (used here as an independent differential-testing oracle
  via LibGEOS.jl).
