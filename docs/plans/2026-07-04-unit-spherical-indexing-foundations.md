# Unit-spherical indexing foundations (2026-07-04)

**Branch: `unit-spherical-indexing`** (off `prepared-minimal`).  Two
foundation stones for spherical preparations: dimension-agnostic bulk-loaded
trees, and robust intersection tests between [`SphericalCap`](@ref)s and 3D
`Extents.Extent`s in unit-spherical space, so caps can drive tree queries
over spherical geometry.

## Findings from the audit

- `FlexibleRTrees` is **already dimension-agnostic**: `STR`'s `_str_tile!`
  recurses over `_ndims(E)` and the Hilbert encoder is Skilling's transpose
  algorithm (any N).  Tests already sweep N ∈ (2, 3).  What remains:
  coverage at N ∈ (1, 4), 3D type-stability checks, and a 3D test for
  `NaturalIndexing.NaturalIndex` (generic in code, only exercised in 2D).
- Still 2D by design, not touched: `NaturallyIndexedRing` (a geographic ring
  wrapper) and the external SortTileRecursiveTree.jl `STRtree` backend
  (2D-hardwired upstream; the in-repo `STR()` loader is the N-D equivalent).
- `SpatialTreeInterface.sanitize_predicate` is the documented hook for
  querying trees with new predicate objects — caps plug in there.

## Cap–extent intersection design

Semantics: the cap is the **closed** region `{p ∈ 𝕊² : p ⋅ center ≥
radiuslike}` — defined by the stored `radiuslike` float, which makes
exactness meaningful.  `Extents.intersects(cap, ext)` is **conservative**:
`false` is a proof of disjointness; `true` means "possibly intersecting".
That is exactly the contract a tree-filter predicate needs (no false
negatives ⇒ no dropped candidates); exact tightness would require comparing
`max {p ⋅ c : p ∈ box ∩ 𝕊²}` against `radiuslike` — a case analysis over
faces/edges of the box with nested radicals — left as future work.

The test is a cascade of three separations, each monotone under box
enlargement (so pruning at union'd node extents never drops a leaf):

1. **Half-space**: the cap lies in `H = {p : p ⋅ c ≥ k}`.  Disjoint if the
   box's support corner `s` (per-axis `hi`/`lo` by `sign(cᵢ)`) has
   `s ⋅ c < k`.  Exact sign.
2. **Sphere shell**: the cap lies on `𝕊²`.  Disjoint if the box's nearest
   corner-to-origin distance² exceeds 1 or its farthest is below 1.  Exact
   signs.
3. **Cap bounding box**: disjoint if `Extents.extent(cap)` (per-axis
   `cos(θᵢ ∓ r)` bounds, outward-rounded) misses the box.  Conservative by
   rounding, catches boxes that straddle the cap plane far from the cap.

**Exactness mechanism**: ExactPredicates' `@genpredicate` does not fit —
its adaptive filter assumes multihomogeneous polynomials, and `s ⋅ c − k`
mixes degrees (the `k` term is degree 0 in the coordinates).  Instead the
two sign kernels (`dot3 − k`, `‖v‖² − 1`) use the same two-stage structure
by hand: float evaluation with a `γ₄`-style forward error bound, falling
back to exact `Rational{BigInt}` arithmetic when inconclusive.  Float
inputs convert to rationals exactly, so the fallback is exact by
construction and testable against big-rational ground truth.

## API additions

- `Extents.extent(cap::SphericalCap)` — conservative 3D AABB of the cap.
- `Extents.intersects(cap, ext)` / `(ext, cap)` — the cascade above.
- `SpatialTreeInterface.sanitize_predicate(::SphericalCap)` — closes over a
  precomputed cap bbox so tree traversal pays no per-node `sqrt`.
- `FlexibleRTrees.query` routed through `sanitize_predicate` (one method),
  so extents, geometries, callables, and caps all work.

## Test plan

- FlexibleRTrees: extend the query ≡ brute-force sweep to N ∈ (1, 2, 3, 4);
  3D `@inferred` construction/query; 3D `NaturalIndex` query ≡ brute force.
- Sign kernels: exact agreement with `Rational{BigInt}` ground truth on
  random and adversarial inputs (support corner dotted to `k` within ±1 ulp).
- `Extents.extent(cap)`: sampled cap points always fall inside, across cap
  sizes (point, small, hemisphere, > hemisphere, full sphere) and centers.
- `intersects` soundness: boxes built around sampled cap points must test
  true; claimed-disjoint pairs verified against a dense spherical lattice.
- Integration: 3D `RTree`/`NaturalIndex` over arc boxes on the sphere,
  queried with caps through `SpatialTreeInterface.query` ≡ linear scan of
  the same predicate; `prepared`/`SpatialTreeInterface`/`UnitSpherical`
  suites re-run.

## Results

Implemented in commits b7561f816 (trees) and ad234a804 (caps); all suites
green on the first run.

- Tree coverage: query ≡ brute force passes at N ∈ (1, 2, 3, 4) for
  `STR`/`HPR`/`Unsorted` across 7 collection sizes; construction and query
  are inferrable in 3D; Hilbert bijectivity + unit-step adjacency hold at
  N ∈ (1, 2, 3, 4); `NaturalIndex` answers 3D queries correctly.  No code
  fixes were needed — the loaders really were dimension-agnostic already.
  `FlexibleRTrees.query` became a re-export of `SpatialTreeInterface.query`
  (they were duplicate implementations), which routes all predicates
  through `sanitize_predicate`.
- Cap–extent: 14 003 sign-kernel assertions match exact rational ground
  truth, including k within ±ulps of the float dot product; 4 002
  witness-soundness checks (no false negative ever); 300 cap/box pairs
  against a 200k-point spherical lattice (margin-certified witnesses all
  detected, 100+ disjoint pairs proven); one-ulp boundary discrimination
  around the shell and the cap plane behaves exactly; cap queries through
  `RTree{STR}`, `RTree{HPR}`, and `NaturalIndex` equal a linear scan of
  the same predicate.
- Deviation from the plan: rather than routing `FlexibleRTrees.query`
  through `sanitize_predicate` as its own method, the duplicate function
  was deleted outright and STI's re-exported.
- Not done here (future work): an exactly *tight* cap–box test
  (`max {p ⋅ c : p ∈ box ∩ 𝕊²}` vs `radiuslike` via face/edge case
  analysis), robustified cap–cap intersection, and arc AABB production
  code (`_edge_extents` analog on the sphere) — that lands with spherical
  `prepare`.

## Arc extents and robust cap predicates (same day, follow-up)

Review feedback: `arc_extent` is a UnitSpherical primitive, not a
`prepare` concern (only the trait-keyed edge enumeration belongs there),
and the angle-space cap–cap/cap–point predicates were unsound — no error
control on `atan`-space comparisons, ~`√eps` error near identical or
antipodal centers, and `radius`-based decisions disagreeing with the
`radiuslike`-based ones at the rim.

- `arc_extent(a, b)`: endpoints' box padded by the arc's sagitta
  `1 − cos(θ/2)` plus rounding slack; antipodal endpoints degrade to a
  cover of every candidate arc.  Exact axis-extrema tightening deferred.
- Cap predicates re-decided in cosine space over the one canonical
  definition (the raw half-space `{p : ‖p‖ = 1, p ⋅ c ≥ k}`, shared with
  the cap–box filter): a float screen at `1e-7` (worst rounding of these
  expressions is ~4e-8, from `√` of cancellation-prone products), then
  exact `Rational{BigInt}` case analysis.  Radicals reduce to one
  squaring level; `rx + ry ≥ π` and `d + rs ≤ π` guards use dedicated
  radical-pair sign helpers.
- Sharp edge found by the consistency tests (and asserted, not hidden):
  `radiuslike = −1` with a center of float norm just over 1 is *not* the
  full sphere — it excludes a ~`√eps` disk at the antipode.  Robust full
  caps need `k ≤ −‖c‖` by margin.
- Benchmarks (M-series laptop, per call over 1000 random pairs):
  cap–cap intersects old 14.6 ns → **new 5.8 ns** (no `atan`); cap–cap
  contains old 14.8 ns → new ~9 ns; cap–point old 14.2 ns → new 7.1 ns.
  Exact fallback costs 4.3 µs (intersects) / 6.8 µs (contains) and is hit
  on ~1e-4 of random pairs (only via the `|k| ≈ 1` degenerate-cap screen
  edge) — sub-ns amortized.  Cap–box filter: 19.5 ns as a sanitized tree
  predicate (25.9 ns self-contained); `arc_extent` 3.4 ns;
  `Extents.extent(cap)` 10.3 ns.
- S2 testing methodology notes (from the reference checkout) recorded in
  `2026-07-04-s2-cap-testing-notes.md`; adopted: two-sided per-axis
  tightness for the cap AABB, adversarial axis-crushed coordinates,
  ±ulp tangency sweeps solved at 512-bit precision against an independent
  angle-space ground truth, log-uniform radii, explicit arc-midpoint
  (max sagitta) checks, degenerate/exactly-proportional arc endpoints.
  ~25.5k assertions added; all suites green.
