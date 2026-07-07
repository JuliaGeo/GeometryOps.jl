# Prepared v2: recursive materialization — decisions (2026-07-03)

**Status: implemented & benchmarked** (commits 9dc75133b, 61d49bb5b,
3727e40f8) — results at the bottom. This amends
`2026-07-03-prepared-minimal-design.md`.

## Direction (agreed)

`prepare` **materializes the geometry into GeometryOps' native layout**
(tuple-vector rings, à la `GO.tuples`) and builds an **eager tree of prepared
nodes** — recursion lives in the *data*, not in access-time wrapping (no
F4-style lazy wrapping / topmost-wins / stripping):

- `Prepared{MultiPolygon}` *stores* `Vector{Prepared{Polygon}}` children;
  `GI.getgeom(p, i)` returns the stored child. Stable identity, zero
  per-access work.
- `Prepared{Polygon}` stores prepared rings ⇒ **the edge index becomes a
  property of the ring** (per-ring `EdgeTree` prep). `RingEdgeTrees` and its
  `hole_trees(trees)[i]` parallel-array pairing dissolve; the PIP seam
  becomes "ask the ring you hold for its tree" (`getprep(ring, ...)`).
- Preparedness survives GI decomposition ⇒ MultiPolygon/GC acceleration falls
  out of the existing polygon seam; ring-level consumers
  (`_line_filled_curve_interactions`, clipping) can discover trees.
- `getprep` API, hooks (`buildprep`, `default_preparations`,
  `build_edge_tree`), and "plain geometry ⇒ `nothing`" contract unchanged.
- Rationale/evidence: the Natural Earth benchmark needed a manual `GO.tuples`
  before `prepare` to be fair — real users would hit that trap.

## Ratified decisions

1. **`Base.parent` returns the converted geometry**, not the original object.
   Documented; no second reference kept. — *agreed as proposed.*
2. **Arbitrary number types, NOT forced Float64.** Tupleization must preserve
   the input's coordinate number type (e.g. `GO.tuples(geom, T)` keyed off the
   input coord type, not hardcoded `Float64`). Do not widen/convert
   coordinates; exactness discussion is moot because we keep the user's
   numbers. Extent/edge-tree code must be generic over the coord type
   (watch: `_ring_edge_extents` currently does `Float64(...)` — remove that;
   Hilbert quantization already normalizes so any Real works).
3. **Heterogeneous GeometryCollections**: abstractly-typed / small-union child
   vectors are fine — take the natural eltype from tupleization or tighten via
   `map(identity, children)`. No special machinery.
4. **No opt-out kwarg** — `prepare` always materializes. No `convert=false`.

## Implementation sketch (next session)

- Rewrite `src/prepared.jl`: `prepare` = recursive build: points/rings →
  tuple-vector storage + per-node preps + cached extent; polygon → prepared
  rings; multi/GC → vector of prepared children (`map(identity, ...)` to
  tighten eltype). GI forwarding now serves STORED children.
- Per-ring prep type (e.g. `EdgeTree <: AbstractPreparation`, kind
  `AbstractEdgeTree`) replaces `RingEdgeTrees`; `default_preparations` on
  `LinearRingTrait` (in polygon context) builds it. Backend selection
  (`NaturalIndex` default / `STRtree` / `FlexibleRTrees` algs / callable)
  carries over via `build_edge_tree`.
- Seam update in `_point_polygon_process` /
  `_point_filled_curve_orientation`: `tree = getprep(ring, AbstractEdgeTree)`
  per ring — simpler than today's `_exterior_tree`/`_hole_tree` adapters
  (delete those).
- Tests: existing equivalence suites already compare original-layout plain vs
  prepared — they guard the materialization exactly. Add: MultiPolygon
  acceleration test (prepared MP children are Prepared), non-Float64 coord
  test (Float32, Int), GC small-union test, `parent` semantics test.
- Re-run Natural Earth benchmark WITHOUT the manual `GO.tuples` on the
  prepared side (keep plain side as-is) — expect the un-tupled prepared path
  to speed up dramatically; that's the headline validation. Update script so
  prepared variants take the RAW geometry.
- Re-run: test/prepared.jl, test/utils/FlexibleRTrees.jl, geom_relations,
  clipping suites (fresh julia, `--project=test`, background w/ EXIT_CODE
  markers).

## Session context for pickup

- Worktree: `/Users/anshul/temp/GO_jts/GeometryOps.jl/.worktrees/prepared-minimal`
  (branch `prepared-minimal` off origin/main, local-only, 8 commits).
- Julia: standard juliaup/homebrew on PATH; if "failed to find source" errors,
  `rm Manifest.toml` and re-instantiate (Dyad-vs-standard mismatch — see
  memory `go-jts-julia-session-quirks`).
- MCP julia session needs restart after adding new `include`d files.
- v1 results (keep for comparison): PIP 41–178 ns prepared vs 0.77–873 µs
  plain; Natural Earth 5.2×/10.9×/19.1×; FlexibleRTrees: HPR wins random
  input ~1000× over unsorted, natural order wins ring edges.

## Results (implemented 2026-07-03)

API as ratified: `prepare(geom; preps)` where `preps` is `nothing`
(defaults everywhere), a `(trait, geom) -> tuple` selector called at every
node (`EdgeTrees(backend)` is the shipped one), or a tuple applied to the top
node only. `EdgeTree`/`AbstractEdgeTree`/`edge_tree` replace `RingEdgeTrees`;
the PIP seam is one per-ring `getprep` (`_point_ring_orientation`), and the
`_exterior_tree`/`_hole_tree` adapters are gone.

- **Micro parity with v1**: same-session comparison shows the v2 machinery
  costs ~10 ns/query over the bare kernel, same as v1. (Caveat for absolute
  numbers: the benchmarking machine ran ~1.55× slower on latency-bound tree
  traversal during the v2 session — verified on untouched code paths — so
  only same-session ratios are meaningful.)
- **MultiPolygon acceleration (new)**: 32 × 1024-gon multipolygon, 1000 pts —
  plain 301 µs/q vs prepared 301 ns/q = **998×**, from forwarding alone; no
  new algorithm code. `prepare` cost 0.49 ms ≈ 1.6 plain queries.
- **Natural Earth from RAW GeoJSON** (no manual `GO.tuples`; the v2 point):
  whole-dataset speedups **5.0× / 10.3× / 18.6×** at 110m/50m/10m with flat
  ~74–88 ns/q prepared; by-size buckets at 10m: 3.2× / 30× / 141× / **591×**
  (>10k verts). Break-even ~4–10 queries/polygon. All correctness gates pass;
  prepared storage is Float32 end-to-end (GeoJSON's number type, preserved
  per decision 2) — the bigger >10k speedup vs v1 (591× vs 481×) is
  consistent with halved memory traffic in the kernel.

## Clipping wired to `AutoAccelerator` (2026-07-03)

The clipping accelerators operate on *rings* (`_build_ab_list` receives
`getexterior` outputs), so prepared rings' `EdgeTree`s slot in directly — no
new preparation, no index composition. `_edge_tree_and_coords(ring, T)` is
the whole seam: reuse the ring's tree as-is when the ring is closed (leaf
indices match `eachedge` one-to-one), else build ephemerally as before.
Genericity: **any** SpatialTreeInterface tree is reused — the accelerators no
longer assume traversal order (candidates are collected and sorted into
nested-loop order, the `SingleSTRtree` pattern), so Hilbert-sorted or
foreign-library trees behave identically to `NaturalIndex`; tests use
HPR-backed ring trees as the out-of-order case. `AutoAccelerator` prefers
tree paths whenever an input ring is prepared (its long-standing TODO).

Measured (1024-vertex circle pair, intersection): `NestedLoop` 22 730 µs;
`AutoAccelerator` plain 65.5 µs; one side prepared 57 µs; both prepared
**48.4 µs** — reuse removes the per-clip edge-list + tree build, ~26 % of the
accelerated clip.

Found by the new tests: union's hole assembly pushed raw *input* rings into
output polygons (three sites; intersection/difference already normalized via
`tuples`) — broken for any non-plain wrapper, `Prepared` included. Fixed by
normalizing through `tuples`/`_linearring` consistently.

## Ring-aware accelerators + line-string preps (2026-07-03, follow-up)

`eachedge` numbers a geometry's edges as concatenated per-curve numberings,
so `_edge_parts` decomposes any geometry into per-curve (tree, coords,
offset) parts and the existing candidate collect-and-sort merges queries
across curves — composition at the *loop* level, no synthetic tree types, no
traversal-order assumptions (the earlier collect+sort refactor is what made
this ~free). Whole polygons/multipolygons through `foreach_pair` (i.e.
`intersection_points`, resolving its acceleration TODO) now reuse prepared
trees and prune curve pairs by extents; curve inputs keep the static fast
path. Line strings get `EdgeTree`s from `prepare` too — indexing exactly
their `eachedge` pairs, no wrap leaf, so reusable even unclosed — with a PIP
seam guard so a wrap-less tree never serves filled-curve orientation. Added
public `hasprep` (node-level boolean mirror of `getprep`); `AutoAccelerator`
uses it instead of a recursive helper — whole geometries fall to the size
heuristic, whose tree paths reuse preps anyway.

Measured (1280-vert holed donut × 1024-vert blob, `intersection_points` on
whole polygons): NestedLoop 14 475 µs; DoubleNaturalTree plain 71.8 µs;
DoubleNaturalTree prepared **23.6 µs** (3× over plain-tree — no edge-list or
tree builds). Boolean intersection on the same pair: prepared 693 µs ≤
plain 709 µs (hole-pass dominated; no regression).

Gotcha for equivalence tests: `_intersection_point` is not bit-symmetric in
its operands — ground truth must be computed with matching argument order.

### Tightening pass (same day)

`_edge_parts` now returns a concretely-typed 1-tuple for curve inputs, so
the parts loops (`_single_parts_loop`, `_dual_pairs_loop`) serve curves and
whole geometries alike; the duplicate `_single_tree_loop`/`_dual_tree_loop`
are gone, along with the single-use `_prep_matches_eachedge` (inlined) and
the `_ring_edge_extents`/`_line_edge_extents` split (merged into
`_edge_extents`).  Kept helpers are function barriers: `_query_part!` /
`_collect_dual_pairs!` hoist per-part dynamic dispatch out of per-edge work.
Benchmarks unchanged (clip 43.4 µs prepared / 62.0 plain; intersection_points
22.8 µs prepared / 70.0 plain-tree, same-session).

Bugs found on the way (all fixed + regression-tested):
- `EdgeTree(geom; backend)` overwrote the default struct constructor and
  broke precompilation → explicit inner constructor.
- The PIP seam handed `Prepared` rings to the kernels, paying a forwarding
  layer per point → `_unwrap_prepared` strips the shell after prep lookup.
- **GeoJSON polygons type their rings as LineStrings**, so ring nodes missed
  the `LinearRingTrait` defaults and lost their edge trees — correct results,
  silently sequential (>10k bucket collapsed to 14×). Polygon children are
  now materialized as `LinearRingTrait` unconditionally. Lesson: a
  correctness gate cannot catch a silent index loss; the per-size timing
  buckets did.
