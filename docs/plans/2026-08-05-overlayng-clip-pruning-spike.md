# OverlayNG clip pruning — a construct-free `RingClipper`

**Status: spike, complete and green. Recommendation: graduate the pruning; do
NOT stop here — the residual cost is now the self-noding pass, not the split.**

GEOS clips both inputs to the intersection envelope before noding
(`RingClipper` / `LineLimiter`), which is most of why its `intersection` is
several times cheaper than its own `union`. We cannot adopt it as written:
closing a clipped ring along the clip box *constructs* coordinates that then
feed the noder, and the one governing decision of this engine (design §0) is
that no constructed coordinate ever feeds a decision.

This spike takes the other half of the idea. **Nothing is clipped.** A parent
segment whose bounding box misses the clip box simply contributes no
`NodedEdge`, and the surviving sub-chains of a ring are left **open** in the
graph — their ends become ordinary degree-1 graph nodes. No coordinate is
synthesized anywhere, and the exactness of the arrangement is untouched.

---

## 1. What changed

Five files, `+424 / −19`, all uncommitted in the working tree.

| file | ± | what |
|:--|--:|:--|
| `src/methods/clipping/overlayng/overlay_ng.jl` | +85 / −9 | the per-op envelope rule |
| `src/methods/clipping/overlayng/noding/split.jl` | +124 / −2 | the prune itself |
| `src/methods/clipping/overlayng/noding/noded_arrangement.jl` | +27 / −8 | plumbing + one new field |
| `src/methods/clipping/overlayng/overlay_labeller.jl` | +21 / −0 | the one pass that needed to know |
| `test/methods/clipping/overlayng/overlay_ng.jl` | +167 / −0 | six new testsets |

### `overlay_ng.jl` — where the boxes come from

* `_overlay_envelopes(m, op, input)` — the input envelopes *this op* can use, or
  `nothing`s. Planar only; skipped for an empty operand (LibGEOS raises on
  `GI.extent` of an empty geometry, and every op that reads them short circuits
  on emptiness first).
* `_overlay_clip_envelopes(op, ea, eb)` — the per-side clip boxes:
  * `OVERLAY_INTERSECTION` → both sides get `env(A) ∩ env(B)`;
  * `OVERLAY_DIFFERENCE` → **B only**, to `env(A)`; A is never pruned;
  * `OVERLAY_UNION` / `OVERLAY_SYMDIFFERENCE` → nothing (GEOS does not clip
    these either);
  * a non-`_OverlayOpCode` `op` (the labeller accepts any
    `(loc0, loc1) -> Bool`) → nothing, by a separate untyped method.
* `_clip_box` normalizes to a concrete `Extents.Extent{(:X, :Y),
  NTuple{2, NTuple{2, Float64}}}`.
* `_empty_result_short_circuit` now takes the two envelopes instead of
  recomputing them, so `GI.extent` still runs **once per side per call** — the
  disjoint-envelope short circuit and the clip pruning share the traversal.

### `noding/split.jl` — the prune

* `NodedArrangement(...; clip_a, clip_b)` threads two optional boxes to
  `_split_edges!(m, table, seg_nodes, segstrings, na, clip_a, clip_b; exact)`.
  `na` is the A/B split point in `segstrings`, so the side of each string is a
  comparison.
* `_seg_in_clip(clip, p0, p1)` — **closed-interval** bbox overlap (a segment
  merely touching the box is kept). `_seg_in_clip(::Nothing, …) = true`, so the
  unclipped path is byte-for-byte what it was.
* A pruned segment does not emit edges and **does not intern its endpoint vertex
  nodes**. A chain end still gets its endpoint node from the surviving
  neighbour.
* `_clip_dropped_points` / `_mark_dropped!` / `_truncated_bits` compute
  `arr.truncated`, the node ids whose star lost an edge (see §2.3).

### `noding/noded_arrangement.jl` — one new field

`NodedArrangement` gains `truncated :: BitVector`. It is `BitVector()` — empty,
allocation-free — for every construction that passes no clip envelope, which is
every existing call site (`antimeridian_split.jl`, `test/.../noding.jl`,
`overlay_graph.jl`, `faces.jl`). No existing construction was touched.

### `overlay_labeller.jl` — the one pass that needed to know

`_compute_labelling!` skips pass 1 (`_propagate_area_locations!`) at nodes
flagged in `arr.truncated`. Guarded by `!isempty(truncated)`, so the default
path is unchanged. Passes 2–5 are untouched. See §2.3 for why this is both
necessary and sound.

---

## 2. The correctness argument, as implemented

Write **E** for the clip box of the side in question.

### 2.0 The premise the driver owns

Every edge that can be in the result lies inside E.

* *Intersection.* The result is contained in `env(A) ∩ env(B) = E`, and so is
  every piece of either boundary that bounds it.
* *Difference (A − B).* The result is contained in A. Its boundary is pieces of
  ∂A (A is never pruned) and pieces of ∂B lying inside A, hence inside
  `env(A) = E`.
* *Union / symdifference.* No such box exists, so no pruning.

### 2.1 The two engine properties that make open chains legal

These are properties of *this* engine, not of JTS/GEOS, which is why the same
trick would be wrong there:

1. **Winding authority is per segment string, not per edge.**
   `EdgeSourceInfo.depth_delta` is derived once at ingest from the ORIGINAL ring
   by `_ring_material_interior_on_left` (`edge_source.jl`). Deleting edges from
   the graph cannot perturb it — `_edge_source_infos` reads `arr.segstrings`,
   which pruning never touches.
2. **Unknown locations fall back to the ORIGINAL input.**
   `_label_disconnected_edges!` locates through `_OverlayInput`'s lazily-built
   `IndexedPointInAreaLocator`s over the original `a`/`b` (design §3
   amendment 7), never over the surviving linework. This is exactly the property
   that makes "prune the whole of B away and still get the right answer"
   work — and it is what the new degenerate tests exercise.

### 2.2 The invariant the prune exports

> **Every node lying inside E keeps its full star.**

A pruned segment's bbox misses E, so *every* point of it — endpoints
included — is outside E. Contrapositively, every segment incident on a point in
E has a bbox meeting E and survives. Two corollaries used below:

* **Every node shared by both inputs is inside E.** A crossing node lies on an A
  segment (⊆ `env(A)`) and a B segment (⊆ `env(B)`), so it is in
  `env(A) ∩ env(B)`; likewise a touch node, which is always an input vertex of
  one side lying on a segment of the other. For difference the same argument
  runs against `env(A)` alone. Hence **a node outside E carries edges from
  exactly one input**, with or without pruning.
* **A pruned segment carries no crossing node**, since a crossing node would be
  in E and would keep the segment.

Consequences, stage by stage:

* **`_merge_noded_edges`** keys on the unordered node-id pair; coincident
  segments have identical bboxes, so they prune together and never
  half-merge. (Confirmed: the "answer-preserving" testset includes
  boundary-sharing pairs, and the whole XML/fuzz/realdata corpus is green.)
* **Star sorting** (`he_order_star!`) is unaffected — it sorts whatever edges
  exist about an exact symbolic apex.
* **Pass 2 BFS** (`_label_connected_linear_edges!`) propagates "same node ⇒ same
  location", which is point-local and does not read the star as a cycle. The
  only thing that could break it is losing a *blocking* boundary edge, and a
  pruned B boundary edge lies entirely outside `env(A)`, so it was never
  incident on an A node in the first place.
* **Pass 3** (collapse ring role) and **pass 5** (PIP against the originals) are
  independent of the graph's completeness.
* **`_mark_result_area_edges!`** cannot mark an edge outside E: for intersection
  such an edge is a kept A edge whose bbox misses `env(B)`, so pass 5 gives it
  `LOC_EXTERIOR` for B; for difference the pruned side is B, and a B edge outside
  `env(A)` gets `LOC_EXTERIOR` for A.
* **Ring assembly** (`_link_result_area_max_ring_at_node!`, maximal/minimal ring
  linking, `_ring_is_collapsed`) only ever touches marked edges, all of which lie
  in E, all of whose endpoint stars are complete. `_build_lines`
  (`oe_mark_visited_both!`) and `_build_points` likewise only see in-result /
  both-input-incident nodes.

**Nothing asserts closure of the INPUT rings anywhere in the graph half.** That
was the spike's main open risk and it did not materialize.

### 2.3 The one thing that did NOT hold: pass 1 reads a star as a closed cycle

`_propagate_area_locations!` walks a node's star CCW carrying a side location and
checks every boundary edge's Right against it, raising
`_OverlayTopologyError("side location conflict")` on a mismatch. Under pruning a
star can be **missing a wedge**, and then the walk reports a phantom conflict on
a perfectly valid input.

This is real, not theoretical. `test/.../overlay_ng.jl` builds the witness: four
valid `MultiPolygon` components of B meeting at one apex far outside `env(A)`,
where exactly one ray of one component's wedge prunes away. Running pass 1 over
that graph **without** the fix raises
`_OverlayTopologyError("side location conflict: arg 1")`; with it, all four ops
match GEOS.

Two details worth recording, because they explain why this is easy to miss:

* Pass 1 never re-checks its own start edge (JTS's `propagateAreaLocations`
  starts with `curr_loc = Left(e_start)` and walks `onext` until it returns).
  A wedge broken at the star's **wrap-around** point is therefore silently
  tolerated. Only a wedge broken in the *middle* of the CCW order raises. The
  first witness I built was accidentally of the tolerated kind.
* The dangerous configuration needs ≥ 2 surviving edges plus ≥ 1 pruned edge at
  one node, i.e. **three or more rings meeting at a point**. An ordinary ring
  vertex has degree 2, so losing one of its segments leaves degree 1, which pass
  1 already skips.

**The fix.** `arr.truncated` marks every node id that a pruned segment
carried — its endpoints (interned only if a survivor also owns them) and any node
the collect stage put in its interior — and `_compute_labelling!` skips pass 1
at those nodes.

Why that is sound:

* **It is complete.** By §2.2 those are the only stars that can have lost an
  edge.
* **It is conservative.** Skipping pass 1 at a node is information-losing, never
  wrong: the edges it would have labelled stay `LOC_NONE` and are resolved by
  passes 2–5, the last of which locates against the ORIGINAL input and is total.
* **It costs nothing that matters.** A truncated node lies strictly outside E, and
  by §2.2 it carries edges from ONE input only — so the star pass 1 is denied was
  never going to relate the two inputs to each other. For the other input,
  `_find_propagation_start_edge` would have returned `0` anyway, pruned or not.

### 2.4 Degenerate cases (the new tests)

All six are in `test/methods/clipping/overlayng/overlay_ng.jl` §5b, and each one
**asserts that pruning fired** by counting the noded edges each side contributed
(`side_edges(arr, is_a)`). Without that assertion an "A ∩ hugeB == A" test is
satisfied by the unpruned pipeline too and proves nothing.

| test | pruning fired | answer |
|:--|:--|:--|
| A strictly inside a huge B | `side_edges(arr, B) == 0` — **all** of B gone | `A ∩ B == A` (from PIP against the original B); `A ∪ B == B`, and the union arrangement has exactly the unpruned edge count |
| A inside B's hole | all of B gone, both for `∩` and for `∖` | `A ∩ B` empty (EXTERIOR from the original B, hole and all); `A ∖ B == A` |
| difference against a far-away B | `clip_a === nothing`; all of B gone, A's 4 edges intact | `A ∖ B == A` |
| MultiPolygon with a far component | B's edges reduced but non-zero | all four ops equal GEOS **and** equal the unpruned pipeline |
| truncated star (§2.3) | `count(arr.truncated) > 0` | all four ops equal GEOS; raises without the pass-1 skip |
| answer-preserving battery | four pairs × four ops | pruned result GEOS-`equals` the unpruned result, edge for edge |

---

## 3. Perf

Natural Earth 10 m, macOS/M-series, Julia 1.12, `Chairmarks` medians over ≥ 2 s
per cell. "pruned" and "unpruned" are the *same* pipeline (arrangement → graph →
labelling → extract) with `clip_a`/`clip_b` set and unset, so the column
difference is the prune and nothing else. LibGEOS is the exactly-matched
reference (GEOS's default overlay since 3.9 *is* OverlayNG).

| case | op | pruned | unpruned | LibGEOS | prune speedup |
|:--|:--|--:|--:|--:|--:|
| France × Italy (4641/3361 v) | intersection | **3.73 ms** | 5.50 ms | 1.22 ms | **1.47×** |
| | union | 5.97 ms | 5.81 ms | 2.06 ms | 0.97× |
| | difference | 5.83 ms | 5.74 ms | 2.01 ms | 0.98× |
| Brazil × Argentina (11121/4674 v) | intersection | **3.97 ms** | 10.15 ms | 1.46 ms | **2.55×** |
| | union | 10.77 ms | 10.68 ms | 4.07 ms | 0.99× |
| | difference | 9.63 ms | 10.53 ms | 3.86 ms | 1.09× |
| Norway × Sweden (15817/4665 v) | intersection | **10.37 ms** | 16.66 ms | 4.72 ms | **1.61×** |
| | union | 17.55 ms | 17.27 ms | 6.54 ms | 0.98× |
| | difference | 17.65 ms | 17.12 ms | 6.56 ms | 0.97× |
| Brazil × 2° box (11121/5 v) | intersection | **1.79 ms** | 7.04 ms | 259 µs | **3.94×** |
| | union | 7.48 ms | 7.39 ms | 2.28 ms | 0.99× |
| | difference | 7.82 ms | 7.61 ms | 2.62 ms | 0.97× |
| Brazil × longest river (11121/1317 v) | intersection | **3.14 ms** | 6.55 ms | 644 µs | **2.09×** |
| | union | 7.06 ms | 6.92 ms | 2.79 ms | 0.98× |
| | difference | 7.03 ms | 6.91 ms | 2.82 ms | 0.98× |

Union and symdifference are not pruned at all, so their ±3 % wobble is
measurement noise; difference gains only where B genuinely leaves `env(A)`
(Brazil × Argentina, 1.09×). **Intersection gains 1.5×–3.9×**, tracking how much
linework the box removes:

| case | op | noded edges kept | truncated nodes |
|:--|:--|--:|--:|
| France × Italy | intersection | 4081 / 7949 (51.3 %) | 6 |
| Brazil × Argentina | intersection | 2795 / 15744 (17.8 %) | 4 |
| Brazil × Argentina | difference | 12922 / 15744 (82.1 %) | 2 |
| Norway × Sweden | intersection | 8719 / 20321 (42.9 %) | 22 |
| Brazil × 2° box | intersection | **4 / 11082 (0.0 %)** | 0 |
| Brazil × longest river | intersection | 2139 / 12416 (17.2 %) | 8 |

### 3.1 Where the remaining time goes — the real finding

Stage breakdown of the extreme case (Brazil × 2° box, intersection: the prune
removes **all but 4 of 11082** segments, yet we are still 7× slower than GEOS):

```
GI.extent(A) + GI.extent(B)          80 µs
ingest both sides (segment strings)  35 µs
arrangement PRUNED                 1.35 ms      (unpruned: 3.69 ms)
  ├─ collect (A×B + self-noding)     983 µs
  │    ├─ A×B dual tree search         0.4 µs
  │    └─ self-vertex pass, A side   947 µs   ←── 71 % of the whole query
  ├─ merge coincident nodes            0.3 µs
  └─ split PRUNED                    293 µs      (unpruned: 2.54 ms)
edge sources + graph + label + build 203 µs
────────────────────────────────────────────
FULL _overlay_ng (pruned)           1.83 ms      LibGEOS: 259 µs
```

The A×B pass costs *nothing* on a pruned workload — the dual-tree traversal
already prunes by extent, so segments outside the box were never enumerated
against the other side. **The residual is `_collect_self_vertex_nodes!`**: the
side's own vertices noded against the side's own segments, walked with
`_self_pair_search` over the *whole* side's index. It is O(all segments) and
completely blind to the clip box. Same picture on Brazil × Argentina: 948 µs of
a 3.95 ms query.

That is the graduation follow-up, and it is a bigger win than this spike's:
building the per-side segment index over the KEPT segments only would remove the
traversal, not just the classification.

**Do not do it without settling this question first.** Self-noding a *kept*
segment against a *pruned* one can find a real node — one ring's vertex lying
strictly inside another ring's segment, both outside E. Dropping the pruned
segment from the index loses that node, which breaks the arrangement invariant
"no node lies strictly inside a noded edge" outside E. My provisional argument is
that it is harmless (an in-result edge has uniform labels along its whole length,
hence lies wholly inside E, hence cannot contain such a node; and the affected
node is already in `truncated`, so pass 1 skips it) — but that is an argument
about label uniformity that this spike did NOT test, and the corpus would have to
be re-run against it. The strictly safe subset — skipping only pairs where BOTH
segments are pruned — is sound with no new argument, but it saves the
classification and not the traversal, which is where the time is.

Two smaller items, both already measured:

* `_clip_dropped_points` + the `_mark_dropped!` probes are the 293 µs of "split
  PRUNED". The first cut looked pruned endpoints up in the node table afterwards
  and cost 968 µs — `NodeKey` has no custom `Base.hash`, so hashing one is nine
  field hashes. Folding the check into the kept pass against a `Set{P}` of
  dropped vertices, with a cheap in-box prefilter so a point inside every active
  box is never hashed, brought it to 293 µs. A custom `Base.hash(::NodeKey)`
  would speed up `_intern_node!` across the whole engine (relate included) and is
  worth its own look.
* `GI.extent` on both operands is 80–114 µs of every planar intersection. It is
  computed once and shared with the disjoint short circuit, so this is not new
  cost for intersection — but `difference` now pays `GI.extent(A)` where it did
  not before. It is repaid many times over whenever B leaves `env(A)`, and is
  under 2 % otherwise.

---

## 4. Gates

Each in its own Julia process, all green.

| test file | exit | notes |
|:--|:--|:--|
| `overlay_ng.jl` | 0 | + 6 new testsets (60 new assertions) |
| `api.jl` | 0 | |
| `noding.jl` | 0 | constructs `NodedArrangement` directly — untouched |
| `overlay_graph.jl` | 0 | ditto |
| `faces.jl` | 0 | ditto |
| `overlay_points.jl` | 0 | |
| `xml_suite.jl` | 0 | exact per-file (pass, skip) pins hold; 715 assertions |
| `realdata_identities.jl` | 0 | worst planar residual 6.57e-15 (unchanged) |
| `fuzz.jl` | 0 | 1596/1600 exact `equals` + valid, 4 benign, 0 divergent |
| `s2_differential.jl` | 0 | spherical is not pruned; 204/204 agree with s2geography |

One real bug was caught by the gates and fixed: computing the envelopes before
the empty-input short circuit made `GI.extent` run on empty operands, and
LibGEOS raises `GEOSGeom_getXMin_r` on those (5 errors in `overlay_points.jl`).
`_overlay_envelopes` now returns `nothing`s for an empty operand.

---

## 5. Open questions

1. **Prune the segment index, not just the split** (§3.1). The single biggest
   remaining win, gated on the kept-vs-pruned self-noding argument above.
2. **Spherical.** Deliberately excluded: the lon/lat box is unreliable across the
   antimeridian and the poles, which is the same reason the disjoint-envelope
   short circuit is planar-only. The natural spherical form is an S2-cell-cap or
   the PR #434 `Extents.extent(Spherical(), …)` crossing-parity extent, not a
   lon/lat rectangle.
3. **Tighter boxes.** For intersection, `env(A) ∩ env(B)` is the cheapest sound
   box, but per-*component* boxes would be much tighter on multipolygons (an
   island of A only needs `env(island) ∩ env(B)`). That is a per-string clip
   rather than a per-side clip, and the `truncated` machinery already generalizes
   to it unchanged.
4. **`NodeKey` hashing.** Nine field hashes per intern, engine-wide. Orthogonal
   to this spike but it showed up as 675 µs of the first cut.
5. **`arr.truncated` as a public-ish invariant.** It is currently a `BitVector`
   field with a documented meaning; if per-component clipping lands it may want
   to become a small struct alongside the boxes, so that `_compute_labelling!`
   stops reaching through `g.arr`.

---

## 6. Verdict

**Graduate.** The prune is 424 lines including tests, defaults to today's exact
behaviour, buys 1.5×–3.9× on planar `intersection` (the most-used op) and nothing
worse than noise elsewhere, and the one invariant it breaks — pass 1 reading a
star as a closed cycle — is fixed precisely rather than papered over, with a test
that fails without the fix.

The honest caveat is that it does not close the gap to GEOS on its own: after
pruning 11078 of 11082 segments we are still 7× slower than LibGEOS on that
query, because the self-noding vertex pass never learned about the box. Ship the
prune, then go after `_collect_self_vertex_nodes!`.

---
---

# Follow-up (same day): the self-noding residual, and `hash(::NodeKey)`

§3.1 and §5 named two items. Both are now implemented and green. Together they
take the extreme case from **1.79 ms to 574 µs** against LibGEOS's 375 µs, and
they make the *unpruned* arrangement build — the path union, symdifference,
`antimeridian_split` and every spherical overlay take — **1.8× faster** on its
own.

Delta on top of the spike above: `+240 / −28` across three new files.

| file | ± | what |
|:--|--:|:--|
| `noding/collect.jl` | +96 / −17 | the clip-aware self-noding traversal |
| `noding/split.jl` | +91 / −2 (net −33 from the spike's version) | `truncated` recomputed positionally |
| `relateng/kernel.jl` | +23 / −0 | `Base.hash(::NodeKey, ::UInt)` |
| `test/.../overlay_ng.jl` | +89 | three adversarial testsets |

## 7. Self-noding under the clip box

### 7.1 What changed

`_collect_crossings!` now takes `clip_a`/`clip_b` and passes each side's box to
**both** self-noding passes — the linear all-pairs pass
(`_collect_self_crossings!`) and the areal vertex pass
(`_collect_self_vertex_nodes!`). Both run through `_self_pair_search`, which
grew a clip argument:

* a tree node whose extent misses the box is never descended into
  (`_ext_in_clip(clip, node_extent(ca)) || continue`), so a subtree of purely
  pruned segments costs O(1) rather than O(subtree);
* a leaf segment whose extent misses the box is dropped before the pairwise
  loop, so it costs no comparisons;
* the pair predicate becomes `intersects(ea, eb) && in_clip(ea) && in_clip(eb)`,
  which also prunes the cross-child `dual_depth_first_search` descent.

`_ext_in_clip` is closed-interval and is applied to `_segment_extent(::Planar, p,
q)`, which **is** the endpoint bbox — the same quantity, same comparison, as
`_seg_in_clip` in split.jl. So collect and split agree segment for segment;
there is no window where a segment is self-noded but emits no edge, or vice
versa.

The A×B pass is deliberately left alone: a pruned A segment's bbox misses
`env(B)`, so the dual traversal never enumerated it in the first place. Nothing
is gained by testing for it and the completeness argument below depends on that
pass staying whole.

**`_collect_self_crossings!` has no consumer outside `_collect_crossings!`**
(checked), so pruning it changes nothing else. Its scope question (linear
strings only) is unchanged.

### 7.2 The argument, and how it was validated

The claim: dropping self-noding pairs with at least one **pruned** member cannot
change the result.

A pruned segment lies wholly outside the box, so any node such a pair would
create is outside the box too. Three cases:

1. **both pruned** — neither emits an edge, so the node would reference nothing;
2. **a pruned segment's vertex inside a KEPT segment** — the kept segment is left
   unsplit there. *This is the only case that changes the graph;*
3. **a kept segment's vertex inside a PRUNED segment** — the pruned segment emits
   no edge either way, so only the node's existence is lost.

Case 2 cannot change the result because the **A×B pass is untouched**: every kept
segment is still split at every incidence with the *other* input, so between
consecutive nodes a kept edge has a uniform location with respect to the other
input, and its label — the ring's own `depth_delta` plus that uniform location —
is correct along its whole length. An edge marked in-result therefore lies wholly
on the result boundary, hence inside the result's closure, hence inside the box;
an edge spanning a dropped same-side node runs through a point outside the box
and so is never marked. What is lost is only the subdivision of **non-result**
linework at same-side touch points outside the box.

Two invariants survive verbatim, and they are what the rest of the pipeline
rests on:

* **The arrangement restricted to the box is identical to the unpruned one.** Any
  incidence at a point inside the box involves only segments whose bbox meets the
  box, and every one of those is still in the traversal. In particular every node
  inside the box still keeps its full star, so ring assembly — which links only
  in-result edges, all inside the box — reads exactly the graph it would have
  read.
* **Coincident segments stay in lockstep.** They have identical extents, so they
  are pruned together *and* noded together; the edge merger never sees a
  half-split pair.

**What did not survive from the coordinator's sketch.** The sketch also proposed
restricting the *vertex* leg — sweeping only vertices whose point lies in the
box. That is not implemented, on purpose: it drops strictly more nodes than the
segment-level filter for no additional traversal saving (the pair has already
been enumerated by the time the vertex leg runs), so it would enlarge the
correctness surface for nothing. The shipped filter is the segment-level one
only.

**What the argument forced.** The spike's `arr.truncated` was maintained
incrementally — it recorded the endpoints and interior nodes of the segments
*split.jl* dropped. Case 3 breaks that: a kept segment's vertex interior to a
pruned segment loses an edge from its star, and it is neither a pruned segment's
endpoint nor (once the self-noding pass is filtered) a recorded interior node. So
the bookkeeping would have gone stale exactly where it matters.

`truncated` is therefore now derived **positionally**: a node is flagged iff its
key is a VERTEX key whose coordinate lies outside a clip box (`_truncated_bits`).
Crossing keys are never flagged and need not be — a crossing lies on a segment of
each input, hence inside both envelopes, hence inside the box. This is a
superset of the truly-thinned stars, and skipping pass 1 at a node is
information-losing only (§2.3), so the superset is sound. It also:

* cannot go stale, whatever a later stage drops — it is a function of the node
  table and the boxes, not a log of decisions;
* costs O(#nodes) branch-free comparisons instead of a `Set` build plus a hash
  per pruned endpoint (the spike's 293 µs on Brazil × 2° box became **13.5 µs**);
* reads an **input** coordinate, not an emitted one — a vertex key stores the
  original vertex verbatim — so design §0 is untouched.

It over-marks only the fringe: a node outside the box that kept its full star.
That fringe is bounded by the kept segments that straddle the box boundary, not
by the interior, and is small in practice (4–22 nodes on the NE cases).

### 7.3 Adversarial tests

Three new testsets, each constructing a **case-2** incidence — a same-side touch
between a kept and a pruned segment, outside the box — and each asserting the
drop actually happened before checking the answer. `self_nodes(A, B, op)` runs
both self-noding passes on a fresh node table with and without the box and
returns the interior-node counts, so `@test c.a_clipped < c.a_plain` is a direct
measurement that the pruned path fired (and `@test c.a_plain > 0` proves the
touch is a real self-node in the first place, so the test cannot pass vacuously).

| test | the incidence | box drops | checked against |
|:--|:--|:--|:--|
| two A components | A2's apex `(15,0)` strictly inside A1's bottom edge; A1's edge spans into the box, A2's apex segments do not | vertex pass | unpruned pipeline + GEOS, all 4 ops |
| shell–hole touch | hole vertex `(15,0)` on the shell's bottom edge | vertex pass | unpruned pipeline + GEOS, all 4 ops |
| self-crossing LineString | A crosses itself at `(15,5)`; the long horizontal segment survives, the vertical one prunes | **linear all-pairs pass** | GEOS, all 4 ops |

All three pass. The wider corpora are the real adversarial coverage and they are
green: `fuzz.jl` runs 400 generated pairs × 4 ops against LibGEOS including a
dedicated *"hole apex on the shell edge (self-touching input)"* class, and
`xml_suite.jl`'s exact pins hold.

### 7.4 Perf

Same run, same machine, four columns measured back to back so they are directly
comparable (this session's machine is ~1.6× slower than the one in §3 — compare
ratios, not absolutes, across sections).

* **v0** — no pruning (pre-spike).
* **v1** — the spike above: clip at split only.
* **v2** — this follow-up: clip at split *and* both self-noding passes (shipped).

| case | op | v0 unpruned | v1 split | **v2 +self-node** | LibGEOS | v0/v2 |
|:--|:--|--:|--:|--:|--:|--:|
| France × Italy | intersection | 6.84 ms | 4.57 ms | **4.10 ms** | 2.00 ms | 1.67× |
| | union | 7.33 ms | 7.40 ms | 7.35 ms | 3.39 ms | 1.00× |
| | difference | 7.17 ms | 7.29 ms | 7.28 ms | 3.16 ms | 0.98× |
| Brazil × Argentina | intersection | 11.89 ms | 4.85 ms | **2.90 ms** | 2.26 ms | 4.10× |
| | union | 12.67 ms | 12.77 ms | 12.78 ms | 6.33 ms | 0.99× |
| | difference | 12.34 ms | 11.20 ms | 10.93 ms | 6.03 ms | 1.13× |
| Norway × Sweden | intersection | 20.87 ms | 13.41 ms | **10.83 ms** | 7.30 ms | 1.93× |
| | union | 22.04 ms | 21.98 ms | 22.05 ms | 10.12 ms | 1.00× |
| | difference | 21.77 ms | 21.82 ms | 21.83 ms | 10.17 ms | 1.00× |
| Brazil × 2° box | intersection | 7.94 ms | 2.21 ms | **573.8 µs** | 374.9 µs | **13.83×** |
| | union | 8.36 ms | 8.50 ms | 8.48 ms | 3.52 ms | 0.99× |
| | difference | 8.76 ms | 8.87 ms | 8.87 ms | 4.08 ms | 0.99× |
| Brazil × longest river | intersection | 6.76 ms | 3.80 ms | **2.11 ms** | 997.8 µs | 3.20× |
| | union | 7.35 ms | 7.40 ms | 7.47 ms | 4.40 ms | 0.98× |
| | difference | 7.41 ms | 7.51 ms | 7.49 ms | 4.43 ms | 0.99× |

All three variants produce GEOS-`equals` results on every cell (asserted in the
benchmark itself, not just the test suite).

The gap to LibGEOS on `intersection` closes from 7× (§3.1) to **1.53×** on the
extreme case and 1.28× on Brazil × Argentina. Stage breakdown of the extreme case
now:

```
GI.extent(A) + GI.extent(B)          138 µs   ←── now the largest single item
ingest both sides                     54 µs
arrangement PRUNED                    61 µs      (was 1.35 ms; unpruned 2.66 ms)
edge sources + graph + label + build  315 µs
────────────────────────────────────────────
FULL _overlay_ng (pruned)            574 µs      LibGEOS: 375 µs
```

The arrangement is no longer the bottleneck at all — it went 1.35 ms → 61 µs,
and `split` alone went 293 µs → 13.5 µs from the positional `truncated`. What is
left is `GI.extent` over both operands (shared with the disjoint short circuit,
so not new cost, but now 24 % of the query) and the labelling/build tail.

## 8. `Base.hash(::NodeKey, ::UInt)`

### 8.1 Equality first

Checked before writing anything, as instructed. `NodeKey` is `isbits`
(72 bytes) and carries **no custom `==` or `isequal`**, so Julia's defaults
apply: `==(x, y)` falls through to `===`, i.e. bitwise equality, and `isequal`
follows `==`. That is correct here *only* because both constructors canonicalize,
and they do:

* `vertex_node(pt)` writes the signed-zero-normalized coordinate into **all four**
  point fields — there are no "unused" fields carrying junk, so two vertex keys
  for the same point are bit-identical;
* `crossing_node(a0, a1, b0, b1)` orders each segment lexicographically and then
  orders the pair, so any order/orientation of the same segment pair yields
  identical bytes (verified at runtime: `crossing_node(a0,a1,b0,b1) ===
  crossing_node(b0,b1,a1,a0)`).

**No pre-existing equality bug.** Nothing to report there.

### 8.2 The method

The default `hash` for such a type is `hash(objectid(x), h)`, and `objectid` of a
72-byte isbits value memhashes all 72 bytes. Measured: **57.8 ns**, against 6.1 ns
to hash one coordinate tuple. Every `_intern_node!` pays it — twice per segment in
`split.jl` alone.

The first version routed each coordinate through `Base.hash(::Float64)`:

```julia
Base.hash(k::NodeKey, h::UInt) =                       # superseded, see 8.2b
    k.is_crossing ? hash(k.b1, hash(k.b0, hash(k.a1, hash(k.pt, h ⊻ …01)))) :
                    hash(k.pt, h ⊻ …00)
```

It hashes exactly what distinguishes the key: the discriminant plus the one
coordinate for a vertex node, all four for a crossing node. Agreement with
bitwise equality is by construction — equal keys have equal `is_crossing`, and
then equal values in every field the taken branch reads. Generic in `P`, so the
spherical `UnitSphericalPoint` path works unchanged.

### 8.2b Hashing raw bits instead

Because the equality being matched is **bitwise** (8.1), the hash contract only
requires agreement with `===` — not with Julia's cross-type numeric `isequal`.
So the coordinates are hashed by their raw bit patterns, which removes
`Base.hash(::Float64)`'s integer-collapse (integer-valued floats route through
`hash(::Int64)`) along with its branch and its `Float64` arithmetic:

```julia
@inline _nk_word(x::Float64) = reinterpret(UInt64, x)
@inline _nk_word(x) = UInt64(hash(x))            # non-Float64 coords: value-determined

@inline function _nk_mix(h::UInt64, w::UInt64)   # shape of Base.hash_mix
    x = widemul(h ⊻ w, 0x9e3779b97f4a7c15)
    return (x % UInt64) ⊻ ((x >> 64) % UInt64)
end

@inline _nk_mix_point(h::UInt64, p) = _nk_mix_point(booltype(GI.is3d(p)), h, p)
@inline _nk_mix_point(::False, h, p) = _nk_mix(_nk_mix(h, _nk_word(GI.x(p))), _nk_word(GI.y(p)))
@inline _nk_mix_point(::True,  h, p) = _nk_mix(_nk_mix(_nk_mix(h, _nk_word(GI.x(p))),
                                                       _nk_word(GI.y(p))), _nk_word(GI.z(p)))

function Base.hash(k::NodeKey, h::UInt)
    x = _nk_mix(UInt64(h), k.is_crossing ? …01 : …00)
    x = _nk_mix_point(x, k.pt)
    if k.is_crossing
        x = _nk_mix_point(x, k.a1); x = _nk_mix_point(x, k.b0); x = _nk_mix_point(x, k.b1)
    end
    return x % UInt
end
```

Three things had to hold, and were checked rather than assumed:

* **Signed zeros are normalized before the key exists.** `vertex_node` and
  `crossing_node` both route every point through `_node_point`, which applies
  `_pos_zero` to each coordinate (`_node_point(::True, …)` rebuilds the 3D point
  from `_pos_zero`-ed `x`/`y`/`z`, so the spherical path is covered too). So
  `-0.0` never reaches a `NodeKey` and bit-hashing cannot split a node in two.
  Verified at runtime: `vertex_node((-0.0,-0.0)) === vertex_node((0.0,0.0))`,
  hashes equal.
* **Finiteness is *not* separately guaranteed** — see 8.4. It does not need to
  be: bit-hashing agrees with `===` exactly, so it is strictly *more* faithful
  than the previous form, which unified bit patterns that are not `===`.
* **Both kernel point types are covered generically.** Coordinates are read
  through `GI.x`/`GI.y`/`GI.z` under the same `booltype(GI.is3d(p))` dispatch
  `_node_point` uses, so the planar `Tuple{Float64,Float64}` folds 2 words and
  the spherical `UnitSphericalPoint{Float64}` (`<: StaticArrays.FieldVector{3,
  Float64}`, fields `x`/`y`/`z`, `src/utils/UnitSpherical/point.jl:41`) folds
  all 3. Any other coordinate type falls back to `hash`.

**The mix matters more than the words.** The obvious `(h ⊻ w) * K` fold is
*wrong here* and was measured to be so: multiplication propagates carries only
upward, and a small-magnitude Float64 carries all its entropy in the high bits
(`1.0` is `0x3ff0_0000_0000_0000`), so the state ends up with ~20 usable bits
and collides at birthday rates — **2193** collisions over the 300×300 integer
grid, **5499** over the fractional grid, far worse than the `Base.hash` form it
replaced. Folding the 128-bit product's halves together with xor diffuses
downward as well and fixes it completely; a `fmix64` finalizer on top changes
nothing measurable and was dropped.

| fold | int grid | fractional grid | random |
|:--|--:|--:|--:|
| `(h ⊻ w) * K` | 87 807 / 90 000 | 84 501 / 90 000 | 90 000 / 90 000 |
| `bitrotate((h ⊻ w) * K, 31)` | 90 000 | 90 000 | 90 000 |
| **`widemul` hi⊻lo fold (shipped)** | **90 000** | **90 000** | **90 000** |

Dict-slot occupancy (90 000 keys into a 2^17-slot table, random-level ≈ 65 100)
is 65 023 / 65 141 for the shipped fold — random-level in the low bits too, which
is what `Dict` actually indexes on.

### 8.3 Measurements

| | `objectid` | `hash(::Float64)` fold | **bit fold (shipped)** |
|:--|--:|--:|--:|
| `hash` of a vertex key | 57.8 ns | 6.5 ns | **3.44 ns** |
| `hash` of a crossing key | 57.8 ns | 16.8 ns | **4.94 ns** |
| `hash` of a spherical vertex key | — | 8.0 ns | **3.82 ns** |
| `Dict` fill, 10 000 vertex keys | 1468 µs | 522 µs | **438.6 µs** |
| collisions, 300×300 integer grid | 0 | 2 | **0** |
| collisions, 300×300 fractional grid | 0 | 0 | **0** |
| collisions, 90 000 random points | 0 | 0 | **0** |

Engine-wide effect, A/B by disabling the method and rebuilding (**unclipped**
arrangement build over four NE 10 m pairs, so this is the path union,
symdifference, `antimeridian_split` and every spherical overlay take):

| pair | planar: none | `hash` fold | **bit fold** | spherical: none | `hash` fold | **bit fold** |
|:--|--:|--:|--:|--:|--:|--:|
| France × Italy | 3.60 ms | 1.63 ms | **1.57 ms** | 18.16 ms | 16.29 ms | **16.12 ms** |
| Brazil × Argentina | 7.65 ms | 3.65 ms | **3.33 ms** | 38.09 ms | 34.20 ms | **33.29 ms** |
| Norway × Sweden | 11.52 ms | 7.15 ms | **6.28 ms** | 45.71 ms | 45.71 ms | **43.47 ms** |
| France shifted-self | 4.50 ms | 2.55 ms | **2.25 ms** | 18.91 ms | 18.95 ms | **16.91 ms** |
| **total** | **27.26 ms** | 14.98 ms | **13.43 ms** | **123.56 ms** | 115.15 ms | **109.78 ms** |

**Planar arrangement build is 2.03× faster than the stock `objectid` hash**
(1.12× over the intermediate `Base.hash(::Float64)` fold). Spherical gains
1.13× — there the exact kernel predicates, not interning, dominate. This is the
single largest whole-engine win in either half of this spike, and it is unrelated
to clipping: it applies to every op on every input.

The integer-grid collisions that the intermediate form exhibited (2 per 90 000,
from `Base.hash(::Float64)` routing integer-valued floats through
`hash(::Int64)`) are **gone** — grid-aligned coordinates are exactly what test
suites and rasterized data produce, so it was worth removing rather than
tolerating, and the fix is faster besides.

### 8.4 One assumption that did not hold: finiteness

Nothing on the ingest path checks that coordinates are finite.
`grep -rn 'isnan\|isfinite\|NaN'` over `relate_geometry.jl`, `kernel.jl` and
`kernel_planar.jl` returns nothing: a NaN or Inf coordinate in the input will be
carried into a `NodeKey` unremarked.

This is a pre-existing property of the engine, not something the bit fold
introduces or worsens, and it does not threaten the hash contract — quite the
opposite. Bit-hashing agrees with `===` *exactly* (equal bytes ⇒ equal words ⇒
equal hash). The form it replaced is the one that unified values which are not
`===`: `hash(NaN) == hash(NaN)` across different NaN payloads, and
`hash(-0.0) == hash(0.0)`. Both unifications are legal (they only merge, never
split), and the signed-zero one is redundant because `_pos_zero` already ran —
but the point stands that the new form needs *no* assumption about the numeric
character of the coordinates. Whether the engine should reject non-finite input
at ingest is a separate question, out of scope here; it is recorded so it is not
mistaken for a consequence of this change.

## 9. Gates (follow-up)

Each in its own Julia process. All exit 0, zero failures. The relateng suites are
included because `hash(::NodeKey)` lives in the shared kernel.

| suite | exit | | suite | exit |
|:--|:--|:--|:--|:--|
| `overlayng/overlay_ng.jl` | 0 (+3 testsets, 38 assertions) | | `relateng/kernel.jl` | 0 |
| `overlayng/api.jl` | 0 | | `relateng/node_topology.jl` | 0 |
| `overlayng/noding.jl` | 0 | | `relateng/predicates.jl` | 0 |
| `overlayng/overlay_graph.jl` | 0 | | `relateng/relate_ng.jl` | 0 |
| `overlayng/faces.jl` | 0 | | `relateng/de9im.jl` | 0 |
| `overlayng/overlay_points.jl` | 0 | | `relateng/allocations.jl` | 0 |
| `overlayng/xml_suite.jl` | 0 (pins hold) | | `relateng/xml_suite.jl` | 0 |
| `overlayng/realdata_identities.jl` | 0 | | `relateng/fuzz.jl` | 0 |
| `overlayng/fuzz.jl` | 0 | | `relateng/kernel_conformance.jl` | 0 |
| `overlayng/s2_differential.jl` | 0 | | `transformations/correction/crossing_edge_split.jl` | 0 |
| `transformations/antimeridian_split.jl` | 0 | | | |

No pin moved and no comparison loosened.

### 9b. Re-run after the bit fold (§8.2b)

Changing the hash changes `Dict` iteration order, which is the one way a
"behaviour-free" hash change can reach output: node/edge visit order, and
therefore result ring and component order. The pin-heavy suites are re-run for
exactly that reason. Each in its own Julia process; all exit 0.

| suite | exit | pass / total | | suite | exit | pass / total |
|:--|:--|:--|:--|:--|:--|:--|
| `relateng/kernel.jl` | 0 | 403 / 403 | | `overlayng/noding.jl` | 0 | 53 / 53 |
| `relateng/kernel_conformance.jl` | 0 | 76 908 / 76 908 | | `overlayng/overlay_ng.jl` | 0 | 675 / 675 |
| `relateng/node_topology.jl` | 0 | 190 / 190 | | `overlayng/s2_differential.jl` | 0 | 260 / 260 |
| `relateng/relate_ng.jl` | 0 | 1883 / 1883 | | `overlayng/xml_suite.jl` | 0 | 744 / 746 (2 pre-existing `@test_broken`) |
| `relateng/fuzz.jl` | 0 | 500 / 500 | | `overlayng/realdata_identities.jl` | 0 | 105 / 105 |
| `relateng/xml_suite.jl` | 0 | 6537 / 6537 | | `overlayng/faces.jl` | 0 | 221 / 221 |

Zero failures and zero errors anywhere; the only non-pass entries are the two
`@test_broken` in the sub-ULP sliver case, which were already broken before this
change. No pin moved.

## 10. Verdict (follow-up)

**Graduate both.**

* **Self-noding prune** — the argument survived adversarial testing; no consumer
  of same-side nodes outside the box was found, so the strictly-safe fallback
  (skip only both-pruned pairs) was not needed. It turns the spike's 1.5×–3.9×
  on `intersection` into 1.7×–13.8×, and brings the extreme case within 1.5× of
  LibGEOS. It also *simplified* the spike: `truncated` lost its incremental
  bookkeeping and became a positional rule that cannot go stale, which is a
  smaller correctness surface than what §2.3 shipped.
* **`hash(::NodeKey)`** — unconditional 2.03× on planar arrangement build,
  engine-wide, with no behavioural surface at all. Equality was already correct
  and canonical; nothing had to be fixed to make the hash valid. Hashing raw bit
  patterns rather than `Float64` *values* (§8.2b) is both the faster and the more
  faithful choice given that the equality is bitwise, and it removes the
  integer-grid collisions entirely.

Remaining, in order of size: `GI.extent` is now 24 % of the extreme-case query
(one traversal per operand, already shared — the fix is an extent cache on the
input, not a cheaper traversal); the labelling/build tail is untouched by any of
this; and open questions 2, 3 and 5 of §5 (spherical, per-component boxes,
`arr.truncated` as a struct) still stand.
