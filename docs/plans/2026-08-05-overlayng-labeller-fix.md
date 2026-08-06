# OverlayNG labeller robustness on near-coincident inputs (F0 / F1 / F2 / F3)

**Date** 2026-08-05/06 · **Branch** `ng-port` · **Files**
`src/methods/clipping/overlayng/overlay_labeller.jl`,
`src/methods/geom_relations/relateng/kernel_spherical.jl`,
`test/methods/clipping/overlayng/labeller_robustness.jl`,
`test/methods/relateng/kernel.jl`

Two independent defects, both needed for the reproducer to pass: three labeller
changes (F0/F1/F2, below) and one kernel change (F3, "The comparator fix").

## The defect

Valid spherical inputs whose boundaries coincide to within ±1 ulp produce genuine
hair-width sliver faces. The exact arrangement resolves them correctly. The
labeller then produced an in-result edge set in which some node had more
result-area boundary arriving than leaving, and the only report of it was
`_OverlayTopologyError("Ring edge missing at max-ring build")` raised several
steps downstream in `maximal_edge_ring.jl` — a symptom with no pointer to its
cause.

### Why degree balance is the right invariant

Once `_mark_result_area_edges!` and `_unmark_duplicate_edges_from_result_area!`
have run, the marked half-edges are exactly the boundary of the result region,
oriented with the region on the right. Walk a node's star as a cycle and let
`b_k` be "the sector between star edge `k` and `k+1` is in the result". Then

    out-degree = #{k : b_{k-1} ∧ ¬b_k}      in-degree = #{k : b_k ∧ ¬b_{k-1}}

which are the enters and the leaves of the same set of sectors around a closed
cycle, hence equal — **provided adjacent star edges agree about the sector they
share**. So an imbalance is never anything but two edges at one node carrying
inconsistent labels, and the ring builder's throw is that inconsistency observed
three data structures later.

## The three changes

### F0 — the balance check itself (`_check_result_area_balance`)

Called at the end of `_unmark_duplicate_edges_from_result_area!`, so it runs on
every overlay. O(#half-edges). On violation it throws naming the node id, its
emitted coordinate, and its whole star with per-half-edge out/in marks, e.g.

    result-area degree imbalance at node 8 (0.9999999999999998, 49.75):
    out=0 in=1; star (edge => dest, out/in): 10=>1 0/1 11=>9 0/0

Truncated nodes are skipped, for the same reason pass 1 skips them: clip pruning
(`split.jl`) thinned their star, so a boundary genuinely can enter one and not
leave. They lie strictly outside the clip box.

### F1 — pass 5 locates on kernel points (`_node_kernel_point`)

JTS has one coordinate per node, so "emitted coordinate" versus "kernel point"
is a distinction that does not exist in `OverlayLabeller`. Here it does, and
design §0 is explicit that no decision may consume a constructed coordinate.
`node_point` is the emission — the substrate's one lossy step — and for this use
it is doubly lossy on `Spherical`: the kernel unit vector goes to lon/lat by
`atan`/`asin` here, and the locator turns it straight back into a unit vector by
`cos`/`sin`. Measured on real reprojected data that round trip moves a
pass-through vertex by up to 14 ulps — an order of magnitude more than the
separations that create slivers in the first place. It is also biased towards the
wrong answer: an emitted coordinate that lands bit-identically on an input vertex
makes the locator report `LOC_BOUNDARY`, which `locateEdgeBothEnds`'
`!= LOC_EXTERIOR` test reads as INTERIOR.

`_node_kernel_point` is dispatched on the node key's point type, which *is* the
manifold's kernel point type:

| node kind | planar | spherical |
|---|---|---|
| vertex | `k.pt` — bit-identical to `node_point` | `k.pt`, the ingested `UnitSphericalPoint` |
| crossing | `node_point` (see below) | `_sph_crossing_dir(True(), k)`, scaled by its largest component then normalized |

**The planar arm is a deliberate no-op.** `RayCrossingCounter` stores its query
point as a `Tuple{Float64,Float64}`, so the exact `Rational{BigInt}` crossing
cannot be handed to it at all; and the emitted coordinate is the *certified*
correctly-rounded image of exactly that rational (emit.jl accepts it only when
the residual plus the dd error bound is below ½ ulp), so it is the best `Float64`
representation that exists. Planar behaviour therefore does not change at all —
asserted directly in the test file, node for node.

The spherical crossing arm scales by `max|dᵢ|` in rational arithmetic before
converting, so a large-numerator exact direction cannot overflow to `Inf` on the
way to `Float64` (a latent hazard `_dir_to_lonlat` in emit.jl also carries).

### F2 — propagation, with PIP demoted to a seed (`_label_disconnected_area!`)

**This is a deliberate divergence from JTS `OverlayLabeller.labelDisconnectedEdges`.**

JTS asks the locator, once per edge and once per input, where that edge lies
relative to the other input's area, and takes the answer as final. Two edges
meeting at one node get two independent answers. At sub-ulp separations those
answers are effectively coin flips, and two of them disagreeing across a single
degree-2 node is precisely an unbalanced node. Measured on the real
reprojected-watershed reproducer: 28 unbalanced nodes out of 2,936. A per-edge
point-in-area rule cannot be made self-consistent by making it more accurate,
because the quantity it is sampling is discontinuous exactly where the sampling
happens.

The divergence is to notice that the answers are **not independent**. A node
carrying no boundary edge of input `gi` is a node at which `gi`'s area location
cannot change: crossing it one stays on the same side of `gi`'s boundary, because
`gi`'s boundary is not there. Every edge incident on such a node therefore has
the *same* `gi` location, and one known value determines all of them — exactly,
not approximately. Propagating it is not a heuristic trading accuracy for
consistency; it reads a value the arrangement already fixed.

Structure, per area input `gi`:

* **Phase A** floods a known location out of every `NotPart`-for-`gi` edge that
  already has one (pass 1 labelled it at its far node) across *transparent*
  nodes. A node is transparent when its star is whole (not truncated) and
  `_find_propagation_start_edge(…, gi) == 0`.
* **Phase B** gives each remaining chain ONE point-in-area verdict, by the same
  rule JTS uses per edge (`locateEdgeBothEnds`: interior iff no endpoint is
  EXTERIOR) generalized to every node of the chain.

A chain can still land on the wrong side of a ±1 ulp separation — irreducible,
and no worse than JTS — but it can no longer land on two sides at once, which is
what broke the builder. It is also strictly less work than JTS: one locator hit
per distinct node of an unresolved chain, against two per unresolved edge.

Two scope decisions inside F2:

* **Non-area inputs keep the old per-edge behaviour** (`set_location_all!(…,
  LOC_EXTERIOR)`). Line inputs already have their own propagation in passes 2/4
  with a `LOC_EXTERIOR`-only rule for `Line` parents; folding them into F2 would
  change line semantics, which this fix has no business doing.
* **A `gi`-collapse edge at a transparent node does not block propagation but is
  never a seed.** A collapse is `gi` linework with zero depth delta, so it does
  not separate — but pass 3 gives it its parent ring's role (hole ⇒ INTERIOR,
  shell ⇒ EXTERIOR), which is a boundary marker, not the location of a
  neighbouring face.

## F3 — the comparator fix (`kernel_spherical.jl`)

**The synthetic reproducer was only half a labeller defect.** With F0/F1/F2 in
place three cells (`n = 7, 31, 63` at `k = 1`) still threw, and the reason was
upstream of the labeller.

Measured on `n = 7, k = 1`, where the sweep reduces to 23 nodes:

* Exact point-in-area of A's own kernel vertices against B says A's shared edge
  is EXTERIOR to B from lat 49.0 to 49.75 and INTERIOR above 49.875 — one sign
  change, on the one A segment that has a proper crossing with B. Fully
  consistent with the exact segment-pair classification (three proper crossings,
  odd parity, matching A's endpoints).
* The graph nonetheless labelled A's segment `[49.25, 49.375]` INTERIOR,
  EXTERIOR, INTERIOR across its three sub-edges where the truth is EXTERIOR,
  INTERIOR, EXTERIOR — because the two crossing nodes on that segment were
  **ordered backwards along it**.
* That order made pass 1's own side locations mutually contradictory across a
  chain of nodes carrying no B edge at all. **No labelling rule can repair it**:
  forcing the chain to one value relocates the imbalance to the crossing node,
  where B's boundary genuinely does separate the two A edges and they *must*
  differ. (Worked through analytically for `n = 7`: forcing gives node 2 three
  marked half-edges, which cannot balance.)

### The defect: a tolerance derived from the quantity it is meant to police

`rk_compare_along_segment(::Spherical, …)` takes the sign of
`disc = (da × db) · N` on the FLOAT node directions and certified it against
`64 * eps * mag`, `mag = |da||db||N|` — computed from those same float
directions. A crossing node's direction is `±(na × nb)` with `na = a0×a1`,
`nb = b0×b1`, which is ill-conditioned in two independent ways:

* `na` is a cross product of two nearly equal unit vectors, so its **absolute**
  error stays at ~ulp while its magnitude shrinks with the segment: a 1e-6 rad
  segment gives `|na| ≈ 1e-6` and ~1e-10 relative error;
* `na × nb` is a cross product of two nearly **parallel** normals when the arcs
  are near-tangent — the near-coincident-boundary case — losing accuracy as
  `1/sin θ`.

`mag` is built from the degraded vectors, so it shrank *with* the accuracy
instead of against it, and the filter certified noise. The old comment claimed
the bound was "amplified for crossing nodes by their arc geometry"; it was not.

### The fix: carry the direction's own relative error

`_float_node_dir_err(k)` returns the float direction together with a bound on its
relative error. From `|Δd| ≤ |Δna||nb| + |na||Δnb| + Δ(na×nb)` and
`|d| = |na||nb| sin θ`:

    ε = (|Δna|·|nb| + |na|·|Δnb| + Δ(na×nb)) / |d|

with each `Δ` a rounding bound from `_cross3_err` (a cross product returning its
own error bound, 1-norm of the per-component `2u·(|aⱼbₖ| + |aₖbⱼ|)`). Both
conditioning failures fall out of the one expression, with no threshold to pick;
a vertex node's direction is its stored coordinate and has `ε = 0`. The filter
then trusts itself only when

    |disc| > (ε_a + ε_b + |ΔN|/|N| + 16·eps) · mag

**This subsumes a `_SPH_TANGENT_GATE`-style hard gate rather than needing one
alongside it.** `|disc| ≤ mag` always, so `ε ≥ 1` makes the bound `≥ mag ≥ |disc|`
and the float path unreachable — and `ε ≥ 1` *is* "no significant digits left",
which is what near-tangency produces. A fixed gate alone would not have sufficed:
a crossing just above a 1e-9 gate still carries ~1e-7 relative error, which the
old `64·eps·mag` bound would have certified regardless.

Measured on the reproducer pair: `ε_a = 285.8`, `ε_b = 483.1`, `|disc| = 8.7e-45`
vs old `tol = 2.2e-58` (certified, wrong) and new `tol = 1.2e-41` (escalates,
correct). Pinned in `test/methods/relateng/kernel.jl`, reconstructed from the
four defining segments so the test depends on nothing but the kernel; verified to
fail on the pre-fix comparator (`1 == -1`, both argument orders, plus the
reversed segment).

### Sibling-filter audit (`kernel_spherical.jl`)

No other float path in the file has this pattern. Checked:

| site | verdict |
|---|---|
| `_on_arc_span_filter` (`_SPAN_ERR_C`) | **Sound** — a Higham running-error bound built from the abs-magnitudes of the terms *before* cancellation, so it grows relative to a cancelled result exactly as it should. The correct shape; the opposite of the defect. |
| `_rk_classify_intersection(True(), …)` | **Sound** — no tolerance at all; the four-orient reduction defers every sign to `rk_orient`'s ExactPredicates filter→exact ladder. |
| `_sph_classify` / `_arcs_cross_properly` / `_sph_compare_around` / `_sph_quadrant3` | **Sound** — pure determinant signs, `Rational{BigInt}` on the `True()` path, no threshold. The `False()` path is approximate by the caller's explicit opt-out. |
| `_ring_is_ccw` / `_spherical_loop_curvature` (`11.25·eps·n`) | **Sound** — S2's `GetCurvatureMaxError`, a per-vertex bound on a Kahan-compensated sum, and `_sph_turn_angle` already uses `robust_cross_product`, which is the standard mitigation for exactly the short-segment cross-product problem the comparator lacked. |
| `_sph_interaction_extent` / `_widen` | **Sound** — extents widened by 4 ulp, conservative in the safe direction. |

Reported, not fixed (different file, different contract): `emit.jl`'s
`_SPH_TANGENT_GATE` tests `|d|² ≥ g²·|na|²·|nb|²` on the same float `na`, `nb`,
so its measured `sin θ` is itself inaccurate for short segments. It is not a
soundness bug — it gates an emitted *coordinate*, which by design no decision
consumes, and its fallback is exact — but the gate is looser than it reads.

### Performance and escalation rate

Arrangement build, old filter vs new, min-of-3 per case:

| case | comparator calls | escalations old → new | build old → new |
|---|---|---|---|
| NE110 neighbours (24 pairs) | 0 | 0 → 0 | 6.65 → 6.61 ms |
| NE110 shifted self-overlaps | 276 | 0 → 0 | 4.26 → 4.27 ms |
| NE10 neighbours (12 pairs) | 0 | 0 → 0 | 65.6 → 65.3 ms |
| NE10 shifted self-overlaps | 6,728 | 0 → 0 | 271.7 → 270.6 ms |
| synthetic sliver sweep (12) | 80 | 24 → 80 | 36.3 → 37.8 ms |

**No escalation-rate change on real data**: across 7,004 comparator calls on
Natural Earth the conditioned bound escalates exactly as often as the old one —
never. Real crossings are transversal and well-conditioned (`ε < 1e-13`, asserted
in the kernel test). The only measurable cost is the pathological synthetic case,
where every call now correctly escalates to `Rational{BigInt}`: +4% on the
arrangement build, +1.5 ms across twelve builds. Country pairs never call the
comparator at all — it fires only when a single segment carries two or more
nodes.

### The sweep, before and after (spherical union, same tree)

| n | k | before (JTS pass 5, old comparator) | after (F0+F1+F2+F3) |
|---|---|---|---|
| 1 | 1 | ok, 1 ring | ok, 1 ring |
| 1 | 2 | ok, 2 | ok, 2 |
| 3 | 1 | ok, 2 | ok, 2 |
| 3 | 2 | ok, 2 | ok, 2 |
| 7 | 1 | THROW `Ring edge missing` | **ok, 3** |
| 7 | 2 | ok, 2 | ok, 2 |
| 15 | 1 | THROW `Ring edge missing` | **ok, 5** |
| 15 | 2 | ok, 6 | ok, 6 |
| 31 | 1 | THROW `Ring edge missing` | **ok, 7** |
| 31 | 2 | ok, 9 | ok, 9 |
| 63 | 1 | THROW `Ring edge missing` | **ok, 11** |
| 63 | 2 | THROW `Ring edge missing` | **ok, 16** |

**12/12 no-throw**, F0 balance clean throughout (the check runs inside every
overlay). Ring counts are unchanged wherever both variants answer. Planar is one
clean rectangle in every cell, before and after.

### Ablation (same sweep, throws out of 12)

| labeller | `rk_compare_along_segment` | throws |
|---|---|---|
| JTS pass 5 | old bound | **5** (`7/1`, `15/1`, `31/1`, `63/1`, `63/2`) |
| F0+F1+F2 | old bound | **3** (`7/1`, `31/1`, `63/1`) |
| F1 only | old bound | 3 |
| F2 only | old bound | 3 |
| JTS pass 5 | F3 | **4** — as `"Ring edge missing"`; F3 alone is *not* sufficient |
| F0 + JTS pass 5 | F3 | 4, now reported at the node |
| F1 only | F3 | **0** |
| F2 only | F3 | **0** |
| F0+F1+F2 | F3 | **0** |

Both halves are load-bearing and together they are sufficient. On this particular
reproducer F1 and F2 are individually enough; they were kept together because the
real reproducer distinguishes them (the diagnosis measured a case that throws with
zero emitted-vs-exact disagreements, which F1 provably cannot reach).

**The real reproducer could not be run.** The brief's watershed WKTs
(`stall_*.wkt`, cached `.jls`, scripts) were to be in a prior session's
scratchpad; that directory exists but is empty, no copy survives anywhere on
disk, and the dataset ("GEOS's own perf dataset, reprojected EPSG:3005→4326") is
not in the `libgeos/geos` tree or findable by GitHub search. The 6/6 fold-order
matrix is therefore unmeasured.

## Test coverage

`test/methods/clipping/overlayng/labeller_robustness.jl`, registered in
`test/runtests.jl`:

* **F1** — planar `_node_kernel_point === node_point` for every node; spherical
  vertex kernel points are `UnitSphericalPoint`s drawn bit-for-bit from the
  input rings, and differ from the lon/lat round trip; crossing kernel points are
  unit.
* **F2** — the invariant directly: at every non-truncated node carrying no
  `gi`-boundary edge, all incident `NotPart`-for-`gi` edges report one `gi`
  location. Checked on the canonical squares (both manifolds × four ops) and on
  the sliver pairs.
* **F0** — balance holds on squares, holed/touching/multipolygon pairs, both
  manifolds, four ops; and a forged single-edge unmark is caught with the node id
  and star in the message.
* **The sweep** — planar exact and single-ringed at every size; spherical
  no-throw in all twelve cells with the area conserved to `rtol = 1e-12`.

`test/methods/relateng/kernel.jl` (appended):

* **F3** — the reproducer pair rebuilt from its four defining segments, with the
  correct order established independently in `Rational{BigInt}` inside the test;
  the comparator's answer in both argument orders and on the reversed segment;
  and the numbers behind the escalation (`ε_a, ε_b > 1`, `|disc| > 64·eps·mag`
  so the retired rule certified it, `|disc| ≤ (ε_a+ε_b)·mag` so the new one does
  not). Plus a well-conditioned control: ordinary transversal crossings keep
  `ε < 1e-13` and stay on the float path, monotone.

Output validity is explicitly not asserted: at these separations the exact union
genuinely contains sliver faces, and how emission renders them is the separate
`_ring_is_collapsed` defect.

## Gates

One file per process, foreground, exit codes checked.

Green (exit 0): `relateng/{kernel, kernel_conformance, node_topology, relate_ng,
fuzz, xml_suite}` — the relate XML suite holding 6537/6537 — and
`overlayng/{noding, overlay_graph, overlay_ng, s2_differential,
labeller_robustness}`. Earlier in the same session, also green:
`overlayng/{overlay_points, faces, api, realdata_identities}` (with
`GO_REQUIRE_DATA=1`).

`overlayng/xml_suite` fails, and did so before this work: the "JTS robust overlay
corpus" ledger now reports six cases NEWLY passing
(`TestOverlay-geos-358.xml` 1, `TestOverlay-misc-1.xml` 2/4/5,
`TestOverlay-misc-2.xml` 2, `TestOverlay-rsf-794.xml` 1). Attributed, not
assumed: reverting BOTH this change's labeller and comparator edits in memory (a
probe, no file touched) reproduces the identical six-case list. It is a
CONCURRENT change to `maximal_edge_ring.jl` / `polygon_builder.jl` /
`collect.jl` / `split.jl` in the same tree, and belongs to whoever owns those.
`overlayng/fuzz`, which failed the same way earlier in the session
(`CENSUS.n_benign <= 8` evaluating `42 <= 8`), is green again as of the final
run. No pin that this change owns has moved.
