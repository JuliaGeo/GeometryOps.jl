# OverlayNG: sub-grid sliver collapse (Defect B)

Status: implemented on the working tree, not committed. Two pins in files outside
this change's scope now fail and are diagnosed but deliberately not re-pinned —
see "What did not survive contact".

## The defect

The 1384-polygon Vancouver watershed cascade (STR leaf order, binary merge tree)
produced, on `Spherical()`:

```
1 polygon, 7 holes, LG.isValid = false
  Ring Self-intersection[-124.123260249449 49.1900388738046]
```

against `Planar()`'s 1 polygon, 0 holes, valid. The seven holes are needles
roughly 110 m long and 1 nm wide. Each is a genuine, non-degenerate face of the
exact arrangement: four distinct exact nodes, exact area between 1.6e-21 and
7.7e-21 sr. They are real geometry, because the operands really do disagree along
those boundaries — a cascade re-ingests every intermediate result through Float64
lon/lat, and only 45% of the corpus's 384,306 vertices survive that round trip
bit-exactly (max drift 2 ulps).

The defect is therefore not that the faces exist. It is that they are emitted.
Diagnosed as **L2 (emission policy)**: a correct exact face whose Float64 image
degenerates undetected. `_ring_add!` drops only *consecutive* repeated points, so
an exact ring `A→B→C→D` whose `B` and `D` round to the same coordinate emits as
`A→B→C→B→A` — five points, so the old collapse test (`length(ring_pts) < 4`)
sees nothing wrong, and an invalid ring reaches the result.

## The fix

`_ring_is_collapsed` gains a second half, `_ring_is_subgrid`, and the call sites
in `polygon_builder.jl` pass the builder context so it can reach the exact nodes.

    drop the ring when     2·|A_exact| < 4 · u · P

* `A` and `P` are the ring's exact area and perimeter over its NODE KERNEL
  POINTS — `node_point` output is never measured (design §0). `2A/P` is the
  ring's mean width.
* `u = _ring_grid_step` is the spacing of representable Float64 coordinates where
  the ring sits: `min(eps(|x|), eps(|y|))` on the plane, and on the sphere
  `min(eps(|lon|)·cos φ, eps(|lat|))·π/180` to convert the emitted degrees into
  the radians the width is measured in. Both take the FINER axis, which
  under-fires rather than over-fires.
* `4` is `_RING_GRID_MARGIN`. Emission moves each vertex by up to ½ ulp per
  coordinate, so two facing vertices can move `√2 · u ≈ 1.41 u` relative to one
  another; 4 is the next power of two above that.

### Arithmetic

Filter/escalate, the shape every predicate in the port uses:

1. certified Float64 pass (`_ring_area2_bounded` on the plane, an inline
   translated cross-product accumulation on the sphere) with an explicit error
   bound; if it certifies either verdict, stop;
2. otherwise escalate the AREA to exact `Rational{BigInt}`.

The perimeter enters as a certified Float64 **lower** bound in both branches
(there is no rational form to escalate to — it is a sum of square roots), which
makes the test fire strictly less often than the ideal one. A ring is dropped
only when `2|A|` is *certified* below `4uP`.

Plain Float64 cannot decide these rings and is not trusted to. A first
calibration pass computed the spherical areas in Float64 and reported ratios two
to four orders of magnitude wrong — the areas are ~1e-22 while the summed terms
are ~1e-5, seventeen digits of cancellation, so the answer was pure rounding
noise and showed *no* separation between the spike population and the legitimate
one. The separation is only visible in extended precision.

On the sphere the area used is the CHORDAL one, `½‖Σ (vᵢ−v₁)×(vᵢ₊₁−v₁)‖`, not
the spherical one. The reason is escalation: the chordal area is a polynomial and
goes to `Rational{BigInt}` exactly, while the spherical area is transcendental
(Girard, or Van Oosterom–Strackee) and would need extended-precision floating
point plus a bespoke error analysis. Chord and arc agree to relative `d²/24`, the
test only fires below `4u ≈ 5e-16` rad, and the corpus's chordal diameters are
~1e-5 — a ~4e-12 relative discrepancy against margins of order 1. It also cannot
misfire on a large ring: for any ring spanning an appreciable part of the sphere
both chordal area and chordal perimeter are O(1), so `2A/P` is O(1), sixteen
orders above `4u`.

## Why this is not the tolerance the design forbids

Design §0 says no decision reads a constructed coordinate and the engine carries
no tolerance. This is the engine's first magnitude-relative threshold, so the
claim needs an amendment rather than a silent exception.

**Proposed amendment text, for §0 or §2.6 of the design doc** (not applied — the
tracked design doc is outside this change's scope):

> **Amendment (emission-grid collapse).** The engine carries no *geometric*
> tolerance: no stage decides that two distinct features are the same feature
> because they are close, and no coordinate is snapped, rounded, or merged to
> make a decision come out. That prohibition is unchanged and unqualified.
>
> Result *emission* is not a decision about the geometry; it is a conversion into
> a format. `Vector{Tuple{Float64,Float64}}` has a finite resolution, `eps(m)` at
> coordinate magnitude `m`, and a face of the exact arrangement narrower than
> that has no faithful image in it — every candidate image either merges its two
> sides or pushes them past each other. Such a face is dropped
> (`_ring_is_subgrid`). This is a statement about the output type, not about the
> sphere or the plane, and it is admitted as such:
>
>   * the length scale is read from the *format at that location*
>     (`_ring_grid_step`), never from a user parameter or a global constant, and
>     it vanishes at the coordinate origin — where the format really can resolve
>     arbitrarily fine detail, and where the test correspondingly never fires;
>   * the decision runs on node ids and kernel points, like every other decision;
>     `ring_pts` is read only to locate the ring on the grid;
>   * nothing upstream of emission may consult it. Noding, labelling and graph
>     traversal stay exact, and a ring dropped here was still built and labelled
>     exactly.
>
> The engine is therefore tolerance-free up to emission, and format-limited at
> it. Those are different claims and the second one is unavoidable for any engine
> that returns Float64 coordinates; what is avoidable, and avoided, is letting
> the second leak upstream into the first.

The honest part of the amendment is the last paragraph of the calibration comment
in `maximal_edge_ring.jl`: there is **no empty band** between the spike population
and the legitimate one. They overlap. Choosing 4 is a policy about which failure
matters more, argued from the format — not a separation discovered in the data.
An earlier draft of this document claimed a "measured empty band 3.16 … 807"; that
was an artifact of measuring only the watershed corpus, and the fuzz corpus
refutes it (legitimate slivers from 0.44 grid steps upward). The claim is
withdrawn.

## Results

### Synthetic reproducer

Two rectangles sharing a subdivided vertical edge, B's interior copies of the
shared vertices alternating −1/+1/−1 ulp across it (`shared_edge_pair` in
`test/methods/clipping/overlayng/overlay_ng.jl`):

| subdivisions | manifold | before | after |
|---|---|---|---|
| 3 | spherical | nring 2, valid | nring 1, valid |
| 5 | spherical | nring 3, valid | nring 1, valid |
| 9 | spherical | **nring 5, INVALID** (`Ring Self-intersection[1 49.4]`) | nring 1, valid |
| 3/5/9 | planar | nring 2/3/5, valid | nring 1, valid |

Planar reproduces the sliver but never the invalidity — see below.

### Watershed cascade

| run | before | after |
|---|---|---|
| spherical 4326, STR order | 1 poly, **7 holes, INVALID** | 1 poly, 0 holes, valid |
| spherical 4326, reversed | 1 poly, 1 hole, valid | 1 poly, 0 holes, valid |
| planar 4326, both orders | 1 poly, 0 holes, valid | unchanged |
| planar 3005, both orders | 1 poly, 0 holes, valid | unchanged |

0 throws in every configuration. Area `3.939112918` (deg²) identical across all
three of {spherical before, spherical after, planar} to 10 significant digits, so
the dropped rings carry no measurable area.

### The seven spike holes

`OLD` is `_ring_image_is_degenerate`; `NEW` is `_ring_is_subgrid`; the last two
columns are what an emitted-image test could have seen.

| hole | npts/distinct | OLD | NEW | signed area (deg²) | image is a legal ring |
|---|---|---|---|---|---|
| 2 | 5/4 | no | **yes** | 0.0 | yes |
| 3 | 5/4 | no | **yes** | 4.547e-13 | yes |
| 4 | 5/4 | no | **yes** | 0.0 | no |
| 5 | 4/3 | no | no | 0.0 | yes |
| 6 | 5/4 | no | **yes** | 0.0 | yes |
| 7 | 4/3 | no | no | 0.0 | yes |
| 8 | 5/3 | no | **yes** | 0.0 | no |

Direct: **5 of 7**. Holes 5 and 7 are 30 µm-wide triangles at mean widths of 807
and 12450 grid steps — far too wide for this test and correctly so. They are
*consequences*: they stop forming once the other five are no longer emitted into
the next cascade level, which is why the end-to-end result is 0 of 7.

Hole 3 is the case that justifies deciding on the exact side. Its image is a
four-vertex ring of positive area that `LG.isValid` accepts; no test of the
emitted coordinates rejects it. Only its exact width (0.089 grid steps) does.

### Ring-drop accounting

| corpus | rings built | dropped by OLD | dropped by NEW only |
|---|---|---|---|
| spherical 4326 fwd | 2974 | 0 | 10 |
| spherical 4326 rev | 3001 | 0 | 2 |
| planar 4326 fwd / rev | 2964 / 2999 | 0 | **0** |
| planar 3005 fwd / rev | 2969 / 3023 | 0 | **0** |

The old test fires on nothing at all in this corpus. All 10 spherical drops have
3–4 nodes and bounding boxes 10–300 m long, i.e. they are the spike population;
none is a plausible real feature.

## Planar reachability

The question was whether the plane can produce a sub-grid ring at all through
`_ring_add!`, and it can — the synthetic above does it on both manifolds, at the
same perturbation, and the keep/drop boundary lands on exactly the same `k`
(between 6 and 8 ulps of edge separation) on both. The threshold is a property of
the output grid, and both manifolds write to the same grid.

What the plane does *not* reproduce is the invalidity. In every planar experiment
run here — the synthetic at 3/5/9 subdivisions and 1–100 ulps, an earlier sweep
at scales 1.0 and 1e6 with 3–31 subdivisions — the sub-grid ring emitted as a
*valid* thin ring, while the spherical one at 9 subdivisions emitted a
self-intersection. The asymmetry is in emission: a planar vertex node emits its
ingested coordinate unchanged and a planar crossing node emits a correctly-rounded
exact rational, so the order of two nearby nodes along an edge survives; the
spherical path converts a unit vector back to lon/lat through trigonometry, which
moved 55% of the watershed corpus's vertices and can reorder them.

This is an observation, not a proof — a planar crossing pair can in principle
round into the wrong order too, and nothing here rules it out. So the test stays
on both manifolds. But it does mean the planar half is currently paying for a
failure only the spherical half has been observed to suffer, which is the
strongest argument available for the follow-up in the next section.

## What did not survive contact

Dropping sub-grid faces means no longer matching GEOS, which emits them. Two pins
outside this change's file scope now fail. Both are diagnosed, neither is
re-pinned.

**1. `test/methods/clipping/overlayng/xml_suite.jl:157** — `LG.equals(ours2, geos2)`
on JTS `TestOverlayMisc` case 2 (GEOS ticket 488). We now differ from GEOS by one
triangle:

```
(2.408076000706469, 48.87511031832609)
(2.40943843942901,  48.87510950043716)
(2.4094384394290103, 48.87510950043716)
```

0.00136 long, area 1.8e-22, and **2.7e-19 wide** — 6e-4 of one grid step, so its
three vertices are distinct Float64s but the shape between them is 1600 times
finer than anything the format can express. Our result stays valid and its area
is unchanged at Float64 resolution. No margin saves this pin: the ring is three
orders of magnitude inside the drop zone. Suggested resolution: weaken that one
assertion to "equal to GEOS up to sub-grid faces", i.e. compare
`LG.symmetricDifference` area against the same rounding band the fuzz suite
already uses, rather than asserting exact `equals`.

**2. `test/methods/clipping/overlayng/fuzz.jl:382** — `CENSUS.n_benign <= 8`
(measured control 4). Now 42. The sweep's own classification is unchanged in
substance: `n_op` 1600, `n_divergent` **0**, `n_broken` **0**, every result valid,
every difference inside the rounding band, worst difference 7.2e-16 in area =
0.4% of the band. What moved is `n_equal` 1596 → 1558: 38 ops that used to equal
GEOS exactly now differ by a sub-grid face. Their widths run 0.44 … 3.98 grid
steps, i.e. right through the drop zone, so this is the margin's direct cost.
Suggested resolution: raise the bound with a comment recording the cause, since
the pin exists to catch *silent* result-dropping and this drop is neither silent
nor a result.

The margin sweep quantifies the trade:

| margin | watershed holes / valid | fuzz ops differing from GEOS |
|---|---|---|
| 0 | 7 / **false** | 0 |
| 1 | 2 / true | 8 |
| 2 | 1 / true | 19 |
| 3 | **0** / true | ~38 |
| 4 | **0** / true | 42 |

The contract (never emit invalid output) is met from 1 up; parity with the planar
census needs 3; the derived value is 4.

## Follow-ups, not done here

* **Decide collapse by comparing the exact ring to its image, not by a width
  threshold.** Drop the ring exactly when the emitted image fails to represent
  it — fewer than 3 distinct vertices, or a signed area whose sign disagrees with
  the exactly-computed `is_hole`, or a non-simple image. That has no constant in
  it at all, so the amendment above would not be needed, and it would keep the
  xml-case-2 splinter (whose image is faithful) while still dropping holes 2, 4,
  5, 6, 7, 8. It would NOT catch hole 3, whose image is faithful-looking and
  whose only defect is exact width — so it is a complement to this test, not a
  replacement, and the pair wants designing together rather than bolting on.
* **Exact spherical kernel positions.** `_node_kernel_point` rounds the exact
  crossing direction to a Float64 unit vector, so on the sphere the "exact" width
  is measured on positions carrying ~1e-16 of representation error — the same
  order as `u` itself. `_ring_is_ccw_exact` already has this property, so this
  test is no worse than the orientation decision beside it, but neither is
  certified at the 1-ulp level and both would want an exact spherical predicate.
* **Anisotropy.** `_ring_grid_step` collapses a two-axis grid to its finer axis.
  A needle aligned with the coarser axis is missed. Comparing the width in the
  direction perpendicular to the needle against the grid step in that same
  direction would be tighter, at the cost of a direction computation.

---

# Round 2: the orphaned-hole throw, and what fixed it

Reported after pin arbitration: JTS robust corpus `TestOverlay-misc-4.xml` case 5
started raising `OverlayTopologyError("unable to assign free hole to a shell")`
under `symdifference` where it had previously failed as `:invalid`. A throw is
strictly worse than invalid output — it breaks the never-throw contract.

## Mechanism, confirmed

The hypothesis (a dropped sub-grid shell orphans its holes) is correct, and it is
one maximal ring in the whole case. Instrumenting `_assign_shells_and_holes!`
over the 34 maximal rings of that symdifference:

```
maximal rings where a SHELL was dropped         : 30
... and holes survived with NO shell (orphans)  :  1
     n_min=2  dropped=1 (shells=1, holes=0)  kept=(shells=0, holes=1)
```

`_find_single_shell` then returns 0, the survivor is appended to
`free_hole_list`, and `_place_free_holes!` finds no containing shell among the
result's other (disjoint) faces, so it throws.

Measuring the offending pair settled the design question empirically rather than
by assumption:

| ring | nodes | bbox | area | verdict |
|---|---|---|---|---|
| shell | 7 | 1010 x 834 | 1.53e-7 | dropped, w/u = 0.237 |
| hole | 3 | 1.4e-9 x 9.3e-10 | 1.08e-19 | **kept** |

`LG.contains(shell, hole)` is true, so the nesting is real. But the hole is a
3-ULP triangle of area 1.08e-19 — it is not a legitimate feature that the shell's
drop stranded. **It should have been dropped on its own merits and was not.**

## Root cause: the perimeter's error bound, not the orphaning

`_ring_width_below` compared `2|A|` against `thr·P`, and `P` came from a Float64
perimeter over `node_point` whose error bound is driven by the ½-ulp displacement
of each endpoint — i.e. by the ring's COORDINATE MAGNITUDE, not by its size. For
this ring at magnitude 3.3e6:

```
perimeter          ~4e-9
its error bound    ~9e-9        <-- larger than the quantity
```

so no positive lower bound on `P` was certifiable, and the code took its
"cannot certify, keep the ring" branch. The ring at the very centre of the test's
purpose was the one shape the test could not evaluate.

**Fix**: when the Float64 filter cannot certify, the escalation now rebuilds the
perimeter from the arrangement's EXACT node coordinates, where there is no
endpoint displacement at all. The exact perimeter is irrational, but a certified
rational LOWER bound suffices — `P` is on the larger side of the comparison, so
under-stating it can only make the test fire less often. Each root is taken in
Float64 and backed off two ulps (`_sqrt_lo`).

A root-free alternative was tried and abandoned, and the failed attempt is worth
recording because it looked strictly better on paper. Substituting
`P ≥ 2·max(bbox extent, longest edge)` removes square roots entirely and
escalates as a pure rational comparison. It is exactly tight for an unsubdivided
axis-aligned needle — and loose for everything else: measured at 0.76 on the
spherical synthetic, whose needle is both diagonal in R³ and subdivided at its
midpoint. The consequences were visible immediately: the watershed cascade's
margin requirement moved 3.0 → 3.5, and the synthetic drop boundary split between
the manifolds (planar k=8, spherical k=6). A threshold that is supposed to be a
property of the output grid must not depend on which way the ring points or how
many vertices sit along its side, so the bound was discarded and `P` kept.

With the exact-perimeter escalation the boundary is back to **k=8 on both
manifolds**, and the margin sweep is back to: validity from margin ≥ 1, planar
census parity from margin ≥ 3.

## The orphan safety net (kept anyway)

The root-cause fix makes misc-4 case 5's hole get dropped on its own, so the
orphan path is no longer reached there. The defensive handling is kept, because
the geometric question the coordinator raised does NOT have a clean "yes":

> Is the hole of a sub-grid shell necessarily sub-grid itself?

For convex rings, yes — `2A/P` is the inradius, and a region inside another has
the smaller inradius. For non-convex rings there is no such monotonicity, and a
compact hole inside a long thin shell is a shape the two statistics genuinely
disagree about. So the code does not assume it. `_assign_shells_and_holes!`
records holes orphaned specifically by a sub-grid shell drop in
`ctx.orphan_hole_list`; `_place_free_holes!` re-offers them to every surviving
shell exactly as before, and only if placement genuinely fails concludes they were
the dropped sliver's interior and drops them. Any OTHER unplaceable free hole
still throws, because that is still a broken arrangement.

Measured on the robust corpus, the re-offer never succeeds (result shells are
disjoint faces, so nothing else contains such a hole). It costs one pass over a
list that is almost always empty, and it is the difference between a justified
drop and an assumed one.

## Ledger

`newly_failing` empty, `newly_passing` empty after the update, all 42 remaining
failures of kind `:invalid`, `worst.rel` = 1.385e-11 (against the 1e-9 bar).

Eight entries removed from `ROBUST_KNOWN_DEFECTS` (50 → 42) — the six identified
during arbitration plus two this round's fixes healed:

| entry | healed by |
|---|---|
| `TestOverlay-geos-358` 1 | sub-grid collapse |
| `TestOverlay-misc-1` 2, 4, 5 | sub-grid collapse |
| `TestOverlay-misc-2` 2 | sub-grid collapse |
| `TestOverlay-rsf-794` 1 | sub-grid collapse |
| `TestOverlay-misc-4` 5 | exact-perimeter escalation (was the THROW) |
| `TestOverlay-qgis-29400` 6 | exact-perimeter escalation |

`TestOverlay-misc-4` case 5 disposition: **healed**. All four ops now valid, no
throw, areas unchanged (`difference` 4.20743940282e8, `symdifference`
5.06358123445e8 — identical to the pre-fix invalid results).

## One more pin moved, and it is out of scope

`xml_suite.jl:261`, the reduced jts-798 case-1 reproducer:

```julia
r_int = go_overlay(:intersection, A, B)
@test 0 < lg_area(r_int) < 1e-11     # now 0.0
```

That assertion pins the OLD behaviour explicitly, and its own comment says so:
"the result is a valid sliver polygon of the needle's true (sub-picounit) area,
where GEOS — having snapped — returns a zero-area MULTILINESTRING". The needle
measures **w/u = 0.038**, twenty-six times inside the drop zone, so the sub-grid
test drops it and the intersection is now empty.

Note the direction: GEOS's answer for this intersection has area **0.0**, and so
does ours now. This pin moving brings us CLOSER to GEOS, not further — the
opposite of the two pins arbitrated last round. Suggested resolution: assert
`lg_area(r_int) == 0` with the needle's measured width recorded, and keep the
`lg_valid` leg.
