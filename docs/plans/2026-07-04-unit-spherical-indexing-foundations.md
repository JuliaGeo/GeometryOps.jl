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

(filled in after implementation)
