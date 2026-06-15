# RelateNG on the Spherical manifold — implementation design

2026-06-15. Builds on the spike `2026-06-12-relateng-spherical-spike.md` /
`.jl`. This document records the validated design decisions; the bite-sized
task list lives in `2026-06-15-relateng-spherical-kernel.md`.

## Goal

Make `relate(RelateNG(; manifold = Spherical()), a, b)` work end-to-end —
Phases 1+2+3 of the spike: the spherical kernel, the acceleration, and the
exact paths — verified by a standalone spherical kernel-conformance suite and
end-to-end tests.

## Decisions (from the 2026-06-15 brainstorm)

1. **Scope: full implementation**, not a stub scaffold. All 13 `rk_*`
   functions for `::Spherical`, including the un-prototyped angle-ordering
   cluster and the `Rational{BigInt}` exact paths.
2. **Conformance: a separate standalone suite**
   (`test/methods/relateng/kernel_conformance_spherical.jl`). The planar
   `kernel_conformance.jl` is left untouched. The spherical suite uses
   exactly-representable integer-xyz inputs in general position (so sign
   predicates stay exact) and carries its own `Rational{BigInt}` differential
   reference for the crossing-apex property.
3. **Boundary: full end-to-end** (Phases 1+2+3). Out of scope: Phase 4's
   broad differential-testing program (JTS-XML port, densified-geodesic
   oracle); we keep validation to the conformance suite plus targeted
   end-to-end smoke tests.
4. **Spherical kernel point type is `UnitSphericalPoint{Float64}`** (not a
   plain `NTuple{3,Float64}`): isbits, `GI.x/y/z`-accessible, and the kernel's
   cross/dot math runs on it directly.

## Architecture

The engine is **already manifold-generic in its control flow**: `RelateNG`,
`RelateGeometry`, `TopologyComputer`, both point locators,
`RayCrossingCounter` and `edge_intersector` all take `m::Manifold`, and the
accelerator selection already dispatches `Planar` vs generic `Manifold`. Two
things are missing:

1. **Point element type `P` is pinned to `Tuple{Float64,Float64}`** at
   construction (`topology_computer.jl:53`, the `Set`/`Dict`/struct fields in
   `relate_geometry.jl`, `point_locator.jl`, `indexed_point_in_area.jl`). The
   structs are *already* parameterized on `P`; only the construction fixes it.
   The enabling refactor derives `P` from the manifold and converts lon/lat →
   xyz once at ingest.

2. **The spherical kernel itself** (`kernel_spherical.jl`) — the 13 `rk_*`
   methods on `::Spherical`.

Correctness vs acceleration are cleanly separable: at `point_locator.jl:525`
non-`Planar` already falls through to `rk_point_in_ring`, and
`edge_intersector` already falls back to `NestedLoop()` for non-`Planar`. So
once the kernel and `P`-threading exist, `relate(Spherical(), …)` is
*correct* (if unaccelerated). Acceleration (3D edge index, lon-interval
point locator, STR, prepared mode) is layered on after.

### Point-type plumbing

- `_kernel_point_type(::Planar) = NTuple{2,Float64}`,
  `_kernel_point_type(::Spherical) = UnitSphericalPoint{Float64}`.
- `_to_kernel_point(m, geopoint)`: identity `(x,y)` on `Planar`; lon/lat → xyz
  via `UnitSphereFromGeographic`, **renormalized** (`normalize`) and
  signed-zero-normalized per component, on `Spherical`.
- `_node_point` / `_node_points` / `_canonical_segment` / `crossing_node` in
  `kernel.jl` generalize from `(x,y)` to N components.
- Thread the derived `P` through `RelateGeometry` →
  `TopologyComputer`/`RelateSegmentString`/`RelatePointLocator`.

### The spherical kernel

Every predicate reduces to a sign of `det(u,v,w) = (u×v)·w`:

- `rk_orient` → `sign((a×b)·c)`; exact via
  `ExactPredicates.orient(tup(a),tup(b),tup(c),(0,0,0))`.
- `rk_point_on_segment` → coplanarity + arc-span dot tests.
- `rk_classify_intersection` → full `SegSegClass`: `SS_PROPER` (candidate
  direction `d = n_a×n_b` strictly interior to both arcs), `SS_COLLINEAR`
  (same great circle, overlapping spans), `SS_TOUCH`/endpoint flags.
- `rk_point_in_ring` → meridian-arc crossing parity to a pole reference;
  pole-insideness derived from the ring's signed area (S2 convention,
  interior on the left of the directed ring).
- angle cluster (`rk_quadrant`, `rk_compare_edge_dir`, `rk_crossing_dirs_ccw`,
  `rk_is_crossing`, `rk_is_interior_segment`) → tangent-plane hemisphere
  split with a reference direction `r` (a pole; fall back when the apex is
  that pole).
- `rk_nodes_coincide` → `Rational{BigInt}` on xyz components.
- `rk_interaction_bounds` → `arc_extent` (spike-proven); area elements
  extended by the six ±eᵢ axis-point test.
- antipodal degeneracy → informative `ArgumentError` naming `AntipodalEdgeSplit`.

### Acceleration

- 3D `Extent{(:X,:Y,:Z)}` flow through `_relate_edge_index` unchanged
  (`NaturalIndex` is dimension-generic).
- `edge_intersector.jl`: `_select_edge_set_accelerator(::Spherical, …)` →
  tree accelerator; `_segment_envs_disjoint(::Spherical, …)` → 3D arc-extent
  disjoint test.
- 3D STR ordering for prepared geometries.
- Lon-interval indexed point locator (`SortedPackedIntervalRTree` over
  longitude intervals; antimeridian crossers split) — an *optimization* over
  the already-correct fall-through.

### Antipodal edges

Kernel throws on exactly-antipodal consecutive vertices. The opt-in
`AntipodalEdgeSplit` correction (in `transformations/correction/`) inserts the
lon/lat midpoint; this is the final, independently-deferrable task.

### Ring/direction containment (`_ring_contains_dir`)

Whether a direction `N` is interior to a ring (S2 interior-on-left) is decided
by the ring's **winding** about `N`: project each edge onto the plane ⊥ `N` and
sum the signed turn. Each edge contributes a turn of magnitude < π, so the sum
has no atan2 branch ambiguity — a per-triangle signed-solid-angle sum
(Van Oosterom–Strackee) does *not* have this property and reports spurious
interiors for reference directions far from the loop's winding axis (e.g. an
equatorial axis against a polar-cap ring). A single interior-on-left enclosure
sweeps +2π, so `N` is interior iff the total exceeds π.

This is exact for rings smaller than a hemisphere — the common geographic case,
and all rings the conformance/end-to-end suites build. **Limitation:** a ring
whose interior is *larger* than a hemisphere (interior = the bigger region)
would under-report a non-encircled interior axis (winding 0 but interior).
`_ring_contains_pole` (used for `rk_point_in_ring`'s pole reference) and the
area-bounds axis extension (`_widen_area_axes`) both rest on this assumption.
The axis extension only over-prunes when wrong, so it stays conservative-safe;
super-hemisphere rings in `rk_point_in_ring` are explicitly out of scope here.

## Out of scope

Phase 4 differential-testing program (JTS-XML spherical port,
densified-geodesic oracle). The conformance suite + end-to-end smoke tests are
the proof for this work.

Super-hemisphere rings (interior larger than a hemisphere) in
`_ring_contains_dir` — see the limitation above.
