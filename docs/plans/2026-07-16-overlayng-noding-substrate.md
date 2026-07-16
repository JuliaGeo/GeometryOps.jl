# OverlayNG port — exact-arrangement design, and the phase 1 (noding substrate) implementation spec

Status: design settled 2026-07-16, after a four-agent survey of JTS `overlayng`/`noding`,
s2geometry's boolean-operation stack, and this repo's foundation, followed by four validation
spikes (S1–S4, verdicts in the appendix). Phase 1 is specified here at implementation depth;
phases 2–3 are outlined so phase-1 authors know their consumer. This document is written to be
consumed by implementation subagents: where it names an existing function, treat the name as a
strong pointer but verify the signature against the code before use.

Reference checkouts (read-only): `/Users/anshul/temp/GO_jts/jts` (JTS, see
`modules/core/src/main/java/org/locationtech/jts/operation/overlayng/` and `…/noding/`),
`/Users/anshul/temp/GO_jts/s2geometry`. This repo: branch `relateng` (tip ≥ `69e416484`).

---

## 0. The one governing decision

Overlay (`intersection` / `union` / `difference` / `symdifference`) is computed as an **exact
arrangement with symbolic nodes, rounded once at emission** — the extension of the RelateNG
port's design D2 to constructive output.

- **No constructed coordinate ever enters a decision.** Node identity, ordering, angular sort,
  labeling, and result assembly are all decided by exact kernel predicates over input vertices
  and symbolic crossing keys. Float64 values appear only as *filters* (certified bounds,
  escalate on uncertainty) and at final *emission*.
- **There is no snapping, no snap-rounding, no precision model, no retry ladder — permanently.**
  JTS's `SnappingNoder`/`SnapRoundingNoder`/`OverlayNGRobust` machinery exists to patch an
  inexact-noding substrate (constructed Float64 points feeding later float decisions). This
  substrate has no such defect, so that machinery is dead weight; do not port it, do not
  reintroduce tolerances. If a fixed-precision *output* mode is ever wanted, it is a
  post-process on the exact result, never part of the computation.
- The pipeline keeps JTS OverlayNG's shape — **node → half-edge graph → label → build → emit**
  — so the engine port stays file-by-file diffable against Java, exactly like the RelateNG port.

Motivation: the current Foster–Hormann clipping's recurring bug class (order-dependent wrong
results and `TracingError`s from inexact constructed intersections + ent/exit tracing —
issues #114, #281, #193, #335, #191, #337, #417, #129). The exact arrangement makes that class
unrepresentable, and the spikes showed it is *faster* than GEOS at the noding stage, not slower.

Evidence base (details in appendix): S1 measured the symbolic noder at **0.03×–0.17× of
LibGEOS's entire intersection op** with zero exact-predicate fallbacks on real data; S2
validated the graph/label/build structural claims end-to-end against LibGEOS bit-for-bit on
planar cases; S3 audited ~930k noded edges and found **zero rounding artifacts on valid
input**, and proved a certified double-double emission path (100% certified, 273× faster than
rational); S4 (landed, `8d938e832`) brought spherical exact classification to planar parity.

## 1. Architecture and phasing

File layout (names mirror JTS where a counterpart exists, for diffability):

```
src/methods/clipping/overlayng/
├── noding/                      # ← PHASE 1 (this spec)
│   ├── noded_arrangement.jl     #   types: NodedArrangement, NodedEdge, node table
│   ├── collect.jl               #   stage 1: candidate enumeration + classification
│   ├── node_identity.jl         #   stage 3: two-tier node grouping
│   ├── split.jl                 #   stage 4: noded-edge emission (topology only)
│   └── emit.jl                  #   coordinate realization (certified fast path + exact)
├── edge_source.jl               # ← PHASE 2: EdgeSourceInfo/depth_delta consumption
├── overlay_label.jl             #   OverlayLabel (plain struct, per JTS)
├── overlay_graph.jl             #   half-edge graph, EdgeMerger by node pair
├── overlay_labeller.jl          #   multi-pass label propagation
├── maximal_edge_ring.jl         #   result ring linking (+ minimal-ring split)
├── polygon_builder.jl           #   shells/holes via rk_point_in_ring
├── line_builder.jl              #   result lines (skip JTS's dead merged path)
├── intersection_point_builder.jl
└── overlay_ng.jl                # ← PHASE 3: isResultOfOp, driver, ops, public opt-in API
```

Kernel additions (phase 1) go in the existing kernel files
(`src/methods/geom_relations/relateng/kernel*.jl`), not in overlayng — they are manifold
predicates, not overlay logic. GeometryOps is a single module, so placement is file
organization only; graduating `noding/` to a top-level `src/noding/` happens if and when a
second consumer (buffer) materializes — not before.

Branch/PR mechanics (decided): branch `overlayng` off `relateng`, one stacked PR per phase.
Phase 1 is internal-only — every new name is un-exported and `_`-prefixed or clearly internal;
no public API or docs-page changes (avoids the docstring `@ref` trap entirely; if any docstring
is written, `@ref` may target **exported names only** — this killed CI twice on relateng).

- **Phase 1 — noding substrate** (§2): geometries → `NodedArrangement` + emission. ~800–1,000
  SLOC + tests.
- **Phase 2 — engine core** (§3): half-edge graph, labels, builders. ~4–4.5k SLOC, faithful
  JTS port with the amendments in §3.
- **Phase 3 — ops + API** (§4): the four ops, mixed dimensions, `OverlayNG{M}` opt-in
  algorithm, differential validation harnesses, PrecompileTools workload, benchmarks.

## 2. Phase 1: the noding substrate

### 2.1 Contract: `NodedArrangement`

```julia
struct NodedEdge
    string_idx :: Int32    # index into segstrings
    seg_idx    :: Int32    # segment within the parent string
    node_lo    :: Int32    # node id at sub-segment start (in traversal order of the parent)
    node_hi    :: Int32    # node id at sub-segment end
end

struct NodedArrangement{P}                 # P = kernel point type; exactly two instantiations
    segstrings :: Vector{RelateSegmentString{P}}   # ingested inputs (NOT text — JTS "SegmentString")
    nodes      :: NodeTable{P}             # see §2.4
    seg_nodes  :: …                        # per-segment ordered interior node-id lists
    edges      :: Vector{NodedEdge}
end
```

Constructor-guaranteed invariants (each one is a test):

1. Every proper crossing between any two indexed segments appears as one shared node id on
   both parents.
2. Node ids are geometrically unique — coincident symbolic keys are merged (§2.4).
3. Per-segment interior node lists are exactly ordered along the parent (§2.5) and deduped;
   no `NodedEdge` is zero-length (S1 finding: one geometric node is routinely reported by
   several candidate pairs — dedup before splitting is mandatory).
4. No Float64-constructed coordinate influenced any of 1–3. Floats appear only inside
   certified filters whose uncertain cases escalated to exact predicates.
5. A `NodedEdge` carries no geometry: its shape is a lookup into its parent segment. All
   source metadata (owner, ring id, dimension) is reached through `string_idx` — nothing is
   copied, so nothing can desynchronize. (Stronger than JTS, which copies a `data` payload
   onto split substrings.)

Naming note: the field is `segstrings`, never `strings` — these are JTS SegmentStrings
(polyline chains of kernel points), and the bare name collides confusingly with `Base.String`.

### 2.2 Ingest

Reuse `RelateSegmentString{P}` unchanged (fields incl. `is_a`, `dim`, `id`, `ring_id`,
`pts::Vector{P}`; built by `_rss_create_ring`/`_rss_create_line`, repeated points removed with
kernel `==` at ingest). Do **not** add fields: `is_hole` is derivable (`ring_id > 1` for
polygon rings — verify the convention in `relate_geometry.jl`), and `depth_delta` belongs to
phase 2's edge-source step, not to the string. Kernel-point conversion happens once here (the
S1 lesson that extraction dominates sparse pairs → keep arrangement construction separable
from ingest so a future prepared overlay converts once; phase 1 just keeps the two steps as
separate functions).

Validity contract: inputs are contractually valid, same as relate — spherical ring
self-crossings are caught by the existing `prepare(...; validate)` path (which uses the same
enumeration machinery), and invalid input is the user's problem to `fix` (`CrossingEdgeSplit`
/ `AntipodalEdgeSplit`). The noder itself performs no validation and no self-noding of A
against A. Input spikes (a–b–a) need no special handling here: they surface as collinear
coincident edges and phase 2's merger resolves them (JTS `DIM_COLLAPSE` semantics).

### 2.3 Stage 1 — collect (`collect.jl`)

- Candidate enumeration: per side, `_relate_edge_index(m, ss_list)` (the existing payload
  `RTree(Unsorted(), owners; extents…)` over per-segment extents — planar boxes /
  `spherical_arc_extent` 3-D boxes), then
  `SpatialTreeInterface.dual_depth_first_search(Extents.intersects, tree_a, tree_b)` for the
  A×B join. Accept optional caller-supplied prebuilt trees (a `PreparedRelate` carries exactly
  this index) — an optional argument, **no** new prepare type, no auto-preparing.
- Per candidate pair: `rk_classify_intersection(m, a0,a1,b0,b1; exact=True())` (post-S4 this
  has certified Float64 fast paths; do not re-filter above it).
  - `SS_PROPER` → `crossing_node(...)` recorded on **both** parent segments.
  - `SS_TOUCH` / `SS_COLLINEAR` → `vertex_node` copies recorded on each segment in whose
    *interior* the flagged vertex lies (the classification's `*_on_*` flags say which). S1
    verified across ~150k real classifications that every touch/collinear intersection is an
    input vertex; keep the defensive assertion anyway — it is cheap and the claim is
    load-bearing.
- Single-threaded (decided). The collect loop is embarrassingly parallel if ever needed; S1
  says it is not where time goes.

### 2.4 Stage 3 — node identity (`node_identity.jl`): two tiers, no canonical key

Node ids are assigned by a two-tier scheme (this supersedes the earlier
`rk_canonical_node_key` idea from spike S2 — computing a canonical rational key costs
rational arithmetic *per node* (same order as emission, 6–18 µs) to optimize a path S1
measured firing **zero** times on real data; the float sweep costs ~1 ms total):

- **Tier 1 (free):** a `Dict{NodeKey,Int32}` — egal-equal keys merge (kernel points normalize
  signed zeros, so vertex-node bit-equality is already canonical).
- **Tier 2 (rare):** geometric coincidence across *different* keys (two distinct segment
  pairs crossing at one point; a crossing landing exactly on a third string's vertex).
  Compute a throwaway Float64 approximation of each crossing node with a certified error
  radius; sort by one coordinate; sweep for overlapping intervals; confirm candidates with
  `rk_nodes_coincide(...; exact=True())`; merge confirmed pairs with a union-find over
  **candidates only** (bounded by candidate count — zero on all real data measured, a handful
  on constructed degree-6 cases).

The approximate position + radius used here is the same computation as the emission fast path
(§2.6) minus the final certification — share the code.

### 2.5 Stage 2 — ordering: new kernel predicate `rk_compare_along_segment`

Added to the `rk_` contract in `kernel.jl` with per-manifold implementations, documented in
the same style as its siblings:

```
rk_compare_along_segment(m, s0, s1, na, nb; exact) -> -1 / 0 / +1
    Order of two nodes (NodeKeys) along the oriented segment (s0, s1).
```

- Float stage: approximate along-segment parameter (planar) / signed discriminant against the
  crossing directions (spherical), each with a certified error bound; adjacent pairs whose
  float gap exceeds the summed bounds are decided in float.
- Exact fallback (lazy — computed only for pairs the filter cannot separate): planar —
  compare exact `Rational{BigInt}` parameters derived from `_exact_crossing_point`; spherical
  — triple-product sign tests against `_sph_crossing_dir` directions. Returning 0 is legal
  and means the nodes coincide — but by construction stage 3 ran first, so equal-order nodes
  have already merged; assert this.
- S1 measured the float filter resolving 100% of 124,500 adjacent comparisons in the dense
  regime (240–590× vs exact-always) and the whole stage a no-op on real data (no real segment
  acquired ≥2 interior nodes). The implementation must be zero-cost for the <2-node case.

### 2.6 Emission (`emit.jl`) — the only lossy step, memoized in the node table

`node_point(arr, id)::P_out` realizes a node's output coordinate on demand and caches it:

- Vertex nodes: the input vertex, bit-exact pass-through.
- Planar crossing nodes: **certified double-double fast path** (S3-proven): TwoSum on endpoint
  differences → compensated 2×2 determinants → dd division → dd recombination; the result is
  *accepted as correctly rounded* iff `|residual| + dd_error_bound < ½·ulp` of the candidate
  (with the determinant-conditioning term, so near-parallel pairs fail the certificate).
  Fallback: `_exact_crossing_point` (`Rational{BigInt}`) rounded to Float64. S3: 64,982/64,982
  real crossings certified, 0 disagreements with the rational answer, 0.024 µs vs 6.44 µs.
- Spherical crossing nodes: dd-certified `na×nb` direction, normalized; lon/lat conversion via
  Float64 `atan2` (documented: the conversion itself is not certified — worst observed
  deviation 1.4e-14° ≈ 1.5 nm; full certification through trig is not worth it because **no
  decision ever consumes emitted coordinates**). Fallback on certificate failure:
  `_sph_crossing_dir` exact direction → high-precision normalize → lon/lat (the
  `CrossingEdgeSplit` precedent, `src/transformations/correction/crossing_edge_split.jl`).

Substrate refactor bundled here (S2 finding g): make `_exact_crossing_point` /
`_sph_crossing_dir` the single shared exact-crossing authority with one call shape, consumed
by (a) this fallback, (b) §2.5's exact ordering keys, (c) `rk_nodes_coincide` — today their
call shapes differ slightly. Grow them in place in the kernel files; no wrappers.

### 2.7 Substrate fix — material-interior authority (for phase 2, landed in phase 1)

`_ring_interior_on_left(m, pts, is_hole; exact)` answers with the **denoted-region** side —
for a hole, the cavity. Overlay labeling needs the polygon's **material-interior** side (the
flip for holes). Add the sibling next to the authority (in `relate_geometry.jl` /
`kernel_spherical.jl`, wherever `_ring_interior_on_left` lives):

```julia
_ring_material_interior_on_left(m, pts, is_hole; exact) =
    is_hole ? !_ring_interior_on_left(m, pts, is_hole; exact)
            :  _ring_interior_on_left(m, pts, is_hole; exact)
```

so the flip exists exactly once and relate/overlay/extents agree by construction. Phase 2
derives JTS's per-edge `depth_delta::Int8` from it (`material_interior_on_left ? -1 : +1`,
matching JTS `Edge.locationLeft/Right`: positive delta ⇒ Left=EXTERIOR, Right=INTERIOR).
Getting this wrong produced a wrong-area hole intersection instantly in spike S2 — write the
regression test from that case.

### 2.8 Explicitly out of scope for phase 1

Half-edge graph, `OverlayLabel`, merger, builders (phase 2). Any op, export, or docs page
(phase 3). Snapping/precision machinery (never). Self-noding invalid inputs (contract, §2.2).
Threading. Prepared-overlay types (optional prebuilt-tree *argument* only).

### 2.9 Tests and gates

New test files under `test/methods/clipping/overlayng/` (one file per julia process when run
locally; they will be wired into the runtests tree):

- **Invariant tests** for §2.1 items 1–5 on constructed cases: two crossing quads; the
  degree-6 node (three segments through one exact point — Tier 2 must merge and the assertion
  in §2.5 must hold); crossing-exactly-on-a-third-string's-vertex; collinear shared boundary
  (GADM-style vertex-identical border → zero phantom crossings, zero interior splits); a–b–a
  spike input.
- **Ordering cross-check**: dense synthetic (comb/zigzag, ~10³ crossings on few segments):
  float-filtered order == exact-always order, elementwise; both manifolds.
- **Emission certificate audit**: for every crossing node in the test corpus, certified
  result == rational result (asserted, not sampled); spherical direction deviation bound
  checked against the exact direction.
- **Rounded-arrangement audit** (S3's census, small): emit all nodes for a small NE-110m
  subset (a few shifted-self pairs), re-classify **crossing-incident** edges on the rounded
  coordinates (the S3 invariant: only crossing nodes move, so only their incident edges need
  auditing), assert zero introduced proper crossings.
- **Classification census invariants** on the NE subset: every touch/collinear carries ≥1
  vertex flag; planar and spherical produce the identical proper/touch/collinear multiset on
  identical data (S1 observed this; encode it).
- Slow/optional (env-gated like the existing realdata benchmarks, not CI): the full NE-110m
  sweep and one GADM pair.

Performance budgets (regression bars from S1/S3, generous ×1.5 for production overhead —
record in the PR, enforce by judgment not CI): NE110 planar sweep noding ~1 ms; GADM CAN×USA
~82 ms planar / ~674 ms spherical; synthetic 250×250 ~19 ms; planar emission ≤0.05 µs/node
fast path.

### 2.10 Conventions for implementers (repo-specific, learned the hard way)

- Floats filter, exact predicates decide; every filter must have a written error-bound
  justification in a comment at the constant's definition (see S4's landed style in
  `kernel_spherical.jl`).
- Engine types stay type-erased over input geometry types (parameterize on kernel point `P`
  only — the Julia 1.12 first-call compile blowup was caused by geometry-typed engine
  internals; do not reintroduce).
- Match JTS names where a counterpart exists; do not invent jargon. Comments state
  constraints, not narration.
- No adapters: if two pieces don't fit, grow the existing abstraction at its source (the §2.6
  and §2.7 refactors are the sanctioned examples).
- Docstring `@ref` targets exported names only. Never `git add -A`. One test file per julia
  process locally; redirect julia output to log files.

## 3. Phase 2 outline (engine core) — for context; port faithfully from JTS with these amendments

JTS sources: `edgegraph/HalfEdge.java` (~500 lines, port required), then the overlayng package
(`Edge`/`EdgeMerger`/`EdgeKey` → `OverlayGraph`/`OverlayEdge` → `OverlayLabel` (plain fields,
NOT bit-packed) → `OverlayLabeller` (5 passes) → `MaximalEdgeRing`/`OverlayEdgeRing`/
`PolygonBuilder` → `LineBuilder`/`IntersectionPointBuilder` → point/mixed paths). Skip:
`ElevationModel`, `FastOverlayFilter` (dead code), strict-mode branches, `LineBuilder`'s
merged path, `RingClipper`/`LineLimiter` (replaced by construct-free whole-ring extent pruning
+ one PIP per pruned ring).

Spike-S2-validated amendments (each was tested end-to-end):

1. Half-edge direction points are the parent segment's **far endpoint** (an original vertex);
   all angular sorting goes through `rk_compare_edge_dir` with a `NodeKey` apex. At a Tier-2
   coincidence-merged node, directions foreign to the representative key take the kernel's
   exact-rational slow path — it exists and is correct; do not special-case.
2. `EdgeMerger` keys on the **unordered node-id pair** (uniqueness of the segment/minor-arc
   between two nodes; antipodal ingest guard covers the sphere) — never on coordinates, and on
   node ids (post-Tier-2), never raw `NodeKey`s.
3. `depth_delta` from `_ring_material_interior_on_left` (§2.7), so `oriented` semantics flow
   through with no further plumbing.
4. Porting trap: `linkResultAreaMaxRingAtNode` must be called **ungated** for every in-result
   edge (JTS's unfulfilled `// TODO: skip already-linked` is deliberate); gating empties
   degree-2 nodes.
5. `PolygonBuilder` containment via `rk_point_in_ring` / the indexed locators — never planar
   even-odd on emitted coordinates.
6. Spherical empty-vs-full: a boundaryless nonempty result is disambiguated by locating one
   input vertex under the op's semantics (topology requirement of the closed manifold, not a
   robustness measure).
7. Labeling's disconnected-edge PIP runs against original inputs via the existing indexed
   locators; symbolic endpoints locate exactly (rare path).

Known future workload for phase 2/3 from the spikes: full minimal-ring split +
hole-nesting (`buildMinimalRings`) — the S3 prototype's union imprecision on 355-island
France and antimeridian Russia came from faking exactly this.

## 4. Phase 3 outline (ops + API)

`isResultOfOp` + driver; the four ops incl. new `symdifference`; mixed-dimension results per
`OverlayUtil` rules; `OverlayNG{M} <: Algorithm{M}` **opt-in** (`intersection(OverlayNG(...),
a, b)`); Foster–Hormann remains the default until differential validation completes.
Validation: JTS overlay JUnit/XML suites; differential fuzz vs LibGEOS (OverlayNG is GEOS's
default engine ≥3.9); spherical differential vs compiled s2 `S2BooleanOperation` (the
Sudan-experiment C++ harness pattern); NE/GADM sweeps with area conservation
(`area(A∪B)+area(A∩B) ≈ area(A)+area(B)` — the signed-area fix `69e416484` makes this a
machine-precision gate). PrecompileTools workload extension; `benchmarks/` additions.

## Appendix — spike verdicts (2026-07-15/16, condensed)

Spike code lives in the session scratchpad
(`/private/tmp/claude-501/-Users-anshul-temp-GO-jts/36b4dc87-418b-4e8c-8a1f-b2b36a1ca475/scratchpad/spikes/{s1_noder,s2_skeleton,s3_census,s4_prefilter,s5_area}/`)
— **ephemeral**: consult if still present, but this document stands alone.

- **S1 (noder throughput)**: symbolic noder on the relateng substrate ran at 0.03×–0.17× of
  LibGEOS's full intersection op (GADM CAN×USA 5.9M edges: 82 ms vs 1.41 s). Zero exact
  fallbacks on real data in ordering and dedup; float ordering filter mandatory in dense
  regimes (240–590×). Touch/collinear = input vertex: zero counterexamples in ~150k. GADM
  shared borders are vertex-identical → zero phantom crossings. Emission was the per-node cost
  center (6/18 µs rational) → S3 fixed. Extraction dominates sparse pairs → keep ingest
  separable.
- **S2 (walking skeleton)**: all structural claims validated; planar results area- and
  vertex-exact vs LibGEOS; crossing emission bit-exact on an integer grid; degree-6 and
  spherical tangent angular ordering correct. Produced amendments §3.1–5 and the
  material-interior fix (§2.7).
- **S3 (rounding census)**: ~930k noded edges audited across NE/GADM/graticule/synthetic/
  spherical: **0** rounding-introduced crossings, **0** float collisions, order inversions only
  on antimeridian-wrapping Russia; area conservation ≤8e-16 (valid outputs); worst LibGEOS
  differential 1.7e-9 on a 1e-4 sliver (their snap-rounding). Only artifacts sat inside
  already-invalid input (NE110 Sudan), where our result matched makeValid-then-overlay while
  raw LibGEOS threw. Certified-dd emission: 100.0000% certified on 64,982 crossings, 0
  soundness violations, 273×; spherical float direction 704× at ≤1.4e-14°.
- **S4 (kernel prefilter — landed `8d938e832`)**: spherical exact classification reduced to
  four `rk_orient` signs (proper-crossing branch, exact path only) + certified span triage +
  repeated-vertex orient short-circuit; 1.65M-pair audit, 0 disagreements; spherical classify
  now ≈0.08 µs/pair (below planar); GADM CAN×USA spherical collect 3.12 s → 674 ms.
- **S5 (signed area — landed `69e416484`)**: `area(Spherical())` Eriksson fan restored to the
  native signed Van Oosterom–Strackee `atan2` form; concave overcounts (Norway +174%, Chile
  +432%) fixed; small-polygon accuracy bit-identical where denom>0; latent denom<0 misbranch
  confirmed and fixed. Relevant here because it makes area conservation usable as a
  machine-precision overlay gate.
