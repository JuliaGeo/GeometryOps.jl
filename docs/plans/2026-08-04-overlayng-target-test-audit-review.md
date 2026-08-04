# Review of the `target` test audit

Reviewed against branch `ng-port` at `91531d75e`. Every load-bearing claim was
re-verified independently — by re-reading the code and by re-running the decisive
experiments from scratch (`review1.jl` in the session scratchpad, same harness,
own fixtures, own mutations via `Core.eval` with restore; working tree left
clean). I did not reuse the audit's scripts.

---

## 1. Verdict

**Yes, with two exceptions.** The audit's findings are accurate — every
experiment I re-ran reproduced its result exactly, including the three headline
ones (the T3 short-circuit non-test, the `isconcretetype(typeof(gc))` tautology,
and the 12-of-13 sweep-pair distinctness). Its recommended replacements for T3
and T9 were verified to work and to discriminate. The exceptions, both in the
direction of keeping slightly more than the audit would:

1. **T7's op axis is not fully inert** (§3.6 of the audit). The structural
   argument — every `_overlay_ng` return site funnels through
   `_target_result` — is correct *for the driver*, but the op axis also probes
   the four distinct public wrapper methods (`api.jl` 232–239, 274–275), each
   with its own `target` keyword forwarding. The audit's replacement probes only
   `GO.intersection`, so an inference regression confined to `union`,
   `difference`, or `symdifference`'s wrapper would go uncaught (T2's runtime
   `typeof` checks cannot see inference-only regressions). Keep one four-op
   inference row: ~56 assertions, not 52.
2. **The T2 trim (27 → 21) should not be applied as written** — or only after
   exception 1 is applied. Its justification leans on T7 covering these argtypes,
   but the audit's shrunken T7 probes `intersection` only, while T2's two empty
   routes go through `intersection` (envelope short circuit) *and* `difference`
   (pipeline-empty). Cutting one `typeof` per shape leaves one route's container
   type resting on inference coverage that the shrunken T7 no longer provides.
   Six assertions at ~0 s runtime; leave T2 at 27 (the audit itself offers this
   as acceptable).

Everything else — the T1 pair deletion and hoist, the T3 rewrite, the T4 loop
deletion, the T5 line-666 deletion, the T8 deletions, the T9 rewrite, and the
`api.jl` block as written — is safe to act on, with one small improvement to the
T1 hoist noted in §4.

---

## 2. Claim-by-claim

### Claim 1 — T3 does not test the short circuit: **CONFIRMED**

Re-ran independently: replicated all 10 of T3's assertions, then redefined
`GO._target_above_result(t, op::_OverlayOpCode, dim_a, dim_b) = false` via
`Core.eval` and replayed through `Base.invokelatest`:

    baseline T3 failures: 0 / 10
    with _target_above_result DISABLED: T3 failures: 0 / 10
    restored: 0 / 10

Every assertion establishes only emptiness, which the full pipeline reproduces.
I also verified both halves of the audit's proposed replacement:

- **Probe (a), antipodal spherical edge**: baseline, untargeted
  `intersection(OverlayNG(Spherical()), Ant, SA)` throws `ArgumentError` while
  the `PolygonTrait()`-targeted call returns `n=0`; with the short circuit
  disabled the targeted call **throws** — so this probe genuinely discriminates,
  which nothing in T3 today does.
- **Probe (b), allocation gate**: measured 128 B (short circuit) vs 21 504 B
  (full pipeline) warm — the audit's exact numbers; the proposed `< 1024` gate
  passes with a 8× margin against the short-circuit path and 168× against the
  pipeline path.

The audit is also right that the comment at 604–605 picks a non-rejected input:
untargeted `intersection(alg, LineString([(0,0),(1,1)]), A)` returns a
`LineString` without error.

### Claim 2 — `isconcretetype(typeof(gc))` is a tautology: **CONFIRMED**

`typeof(x)` returns a concrete type for every Julia value; the assertion cannot
fail for any behaviour of the code under test. Re-ran the `Vector{Any}`
mutation myself (reverted `_create_result_geometry`'s mixed branch to
`comps = Any[]`):

| assertion | baseline | with `Vector{Any}` |
|---|---|---|
| 758 `trait isa GeometryCollectionTrait` | pass | pass (no catch) |
| 761 `eltype(...) === _ResultComponent` | pass | **fail (catches)** |
| 762 `isconcretetype(typeof(gc))` | pass | pass (no catch) |
| 767 `all(isa _ResultComponent, ...)` | pass | pass (no catch) |
| T7:753 `all(isconcretetype, ms)` | pass | **fail (catches)** |

Identical to the audit's experiment 20. One nit: the audit cites the tautology
as line 763; in the tree it is line 762 (its other citations — 758, 761, 767 —
are correct). The 767-subsumed-by-761 argument is sound: the components come
back as the collection's own `Vector`, whose eltype invariant makes 767
unfailable while 761 passes.

### Claim 3 — 12 of 13 sweep pairs distinct; only `"overlap+boundary"` duplicates: **CONFIRMED**

Recomputed the branch signatures from scratch (driver path × per-op component
counts at each dimension):

    area x area         AR|1,0,0,1,0,0,1,0,0,2,0,0
    overlap+boundary    AR|1,0,0,1,0,0,1,0,0,2,0,0   <- only duplicate
    ... (11 others, all distinct)
    distinct signatures: 12 of 13

Exactly the audit's result, including the byte-identical duplicate signature.
The comment at line 512 ("a mixed-dimension result") is confirmed false: the
intersection is one polygon, zero lines — as the comment at 498–500 in the same
block already says. The signature argument is sound for this testset
specifically because T1 compares targeted output against *untargeted output of
the same engine* — an engine bug affects both sides, so T1 only ever tests the
target filtering machinery, whose branches the signature fully captures. The
deletion is safe.

### Claim 4 — T7's axes are inert; 208 → ~52: **PARTLY RIGHT**

Confirmed, by measurement:

- Op axis, one cell (Planar, Poly×Poly, `PolygonTrait()`): all four wrappers
  infer to one concrete type.
- Untargeted return types: `Poly,Poly` / `Poly,MPoly` / `Line,Poly` are each a
  7-member `Union`; `MPoint,Poly` and `MPoint,MPoint` are `Any`. The 4→3
  argtype collapse and the block-comment fix ("7-member Union" is only true for
  edged inputs) both check out.
- The missing `MultiPoint × MultiPoint` case: confirmed absent today, and
  confirmed that the audit's proposed addition **passes** — all 6 targets infer
  concretely on that argtype, to the same type as `(A, B)`. So the replacement
  is viable as written.

The exception (see Verdict): "the op axis cannot affect the answer" is proven
only about `_overlay_ng`'s internals. The four public wrappers are four separate
methods, each separately inferred; the op axis is the only place in the suite
that asserts *inference* (not just runtime types) through `union`, `difference`
and `symdifference`'s keyword forwarding. The realistic regression class is an
edit to one wrapper (validation branch, target normalisation at the wrapper
level) that widens its inferred return type without changing runtime values —
invisible to T2 and T9. Cost of retaining it: one row of 4 assertions at a
single (target, manifold, argtype) cell. Recommended shape: the audit's 48 + 4
(four-op row) + 4 trailing = **56**.

### Claim 5 — T5 line 666 duplicates line 391: **CONFIRMED**

By reading: line 651's `A` is character-identical to line 380's `A`; lines
665–666 construct the identical
`_OverlayInput(Spherical(), A, A, 2, 2, EX, false, false, nothing, nothing)`
and assert the identical
`@test_throws ArgumentError _resolve_empty_result(Spherical(), OVERLAY_UNION, inp)`
as lines 389–391. Deleting 666 is safe; the untargeted-control role is served by
667–670 two lines below.

### Claim 6 — T9's symdifference rows cannot discriminate: **CONFIRMED**

Measured untargeted traits on T9's fixture: intersection/union/difference →
`PolygonTrait()` (line 776 discriminates for those), symdifference →
`MultiPolygonTrait()` in **all four** forms (alg, bare, `Planar`, `Spherical`) —
so for symdifference both the trait assertion and the area assertion pass even
if `target` were silently dropped, in lines 776–777 (symdiff row) and 780–784
entirely. The area assertions are non-discriminating for all four ops. The
proposed replacement (`PolygonTrait()` → `isa AbstractVector`) discriminates,
since no untargeted call returns a `Vector`.

### Claim 7 — 926 today, ~426 recommended: **CONFIRMED** (arithmetic)

Recounted the inventory per testset from the source: 625 + 27 + 10 + 11 + 12 +
16 + 208 + 6 + 11 = 926 — every row of the audit's §2 table matches. The §6
column sums to 426 as claimed. With my two exceptions the target becomes ~436
(§5 below).

### Claim 8 — interface-vs-internals in good shape: **CONFIRMED**

Spot-checked every row of the §4 table: the wrappers are indeed one-line
pass-throughs (`api.jl` 232–239, 274–275); `TraitTarget` is exported public
surface (GeometryOpsCore `traittarget.jl` line 8); and grep over lines 464–785
confirms the audit's two corrections to its brief — neither `_empty_geom` nor
`_covers_everything` is reached by any new test. The characterisation of the
four deep reaches as elision/guard tests is fair. One small slip in §3.10, see
§4 below.

### Claim 9 — `_FULL_SPHERE_MSG` publicly unreachable: **CONFIRMED**, and strengthenable

I tried my own route rather than accepting the audit's hemispheres experiment: a
three-lune cover (`a` = one 120° lune, `b` = MultiPolygon of the other two, all
meridian boundaries shared pairwise, union = whole sphere with every boundary
cancelling). Result: the pole-vertex lunes are degenerate in GO's spherical
model (`GO.area(Spherical(), lune)` = 0.0) and the union returns a
`MultiLineString` — no full-sphere throw, no areal content at all. Same failure
mode as the audit's great-circle hemispheres: the representable spelling of a
cancelling cover is exactly what the model lacks.

There is also a structural argument the audit did not make, which upgrades
"could not find one" to "cannot exist for valid inputs": under enclosed-region
semantics every non-degenerate ring encloses < 2π, and a polygon is a disk
minus disks. For a boundaryless full-sphere union, the (connected — a k-holed
sphere) complement of `int(a)` must lie inside one component of `b`, hence
inside that component's exterior-ring disk `D` with area(D) < 2π; then the
complementary disk of `D`, with area > 2π, must lie inside a single component
of `a` — forcing both to be exact hemispheres, whose great-circle ring is
degenerate in this model. So only invalid inputs (edge-sharing or overlapping
multipolygon components) could even in principle reach the gate, and both
attempts at those collapse to degeneracy before noding completes. The internal
`_OverlayInput` construction in T5 is the right instrument; keep it.

### The audit's §5 ("do not delete") — accepted

§5.1's logic is airtight without re-running its mutation: the sweep's coordinate
comparison compares empty lists as equal regardless of container trait, so
`GI.trait(mul) === multi` is the only absolute container pin — the hoist must
not become a deletion. §5.3–§5.5 are keep-recommendations (no coverage risk if
the audit were wrong), and their structural claims check out against the pairs
list and `overlay_mixed_points.jl`; I verified the sole-catcher pair claims
structurally (only `(P, P2)` is dim-0×dim-0; only `(P, L)` has a line non-point
side; only `(A, Etouch)` is an areal pair with a point intersection).

---

## 3. Recommendations I would NOT act on as written

1. **T7: 208 → 52** (§3.6). Act on it as 208 → **56**: keep one four-op
   inference row (e.g. `PolygonTrait()`, `Planar`, `Poly×Poly`) alongside the
   audit's intersection-only sweep. Coverage otherwise lost: inference of the
   `union`/`difference`/`symdifference` wrappers' `target` forwarding — the only
   inference probes those three methods have anywhere in the suite.
2. **T2: 27 → 21** (§3.9). Keep at 27, or apply only together with (1). Coverage
   otherwise lost: the absolute container type of one of the two empty routes
   (envelope short circuit vs pipeline-empty), whose backing inference argument
   assumed op coverage that the shrunken T7 no longer has. The audit itself
   flags this trim as optional; take the option.

Both exceptions are additive; no audit deletion needs to be *reversed*.

---

## 4. What the audit missed

Nothing that changes a verdict; three genuine misses and one replacement flaw:

1. **T6 has a small redundancy after all** (§3.10 says it has none, and that the
   four accepted spellings "exercise four distinct `_ov_target` methods").
   Measured: `GO.TraitTarget(GI.PolygonTrait()) === GO.TraitTarget(GI.PolygonTrait)`
   — spellings 3 and 4 construct the identical value, so they exercise the same
   `_ov_target` method (three distinct methods total, not four), and 2 of T6's
   16 assertions are literal duplicates. (Lines 696 vs 697 are *not* duplicates
   — different Unions — that part of the audit holds.) Defensible to keep both
   spellings as user-facing documentation; the audit's justification for
   "KEEP all 16" is nevertheless wrong as stated.
2. **The §3.2 hoist has no non-empty dim-1 fixture.** With `(A, B)` and
   `(A, Etouch)` under intersection only, the lines list is empty in every cell
   (verified from the signature table: `AR|1,0,0…` and `AR|0,0,1…`), so the
   hoisted `GI.trait(mul) === multi` never pins a *non-empty*
   `MultiLineString`, and `vec isa AbstractVector` never sees a non-empty
   line vector. T2's relative `typeof(e) === typeof(nonempty_m)` closes most of
   the resulting gap (another reason not to trim T2), but the fix is one line:
   add `(A, Cedge)` as a third fixture pair — 18 hoisted assertions instead
   of 12.
3. **The wrapper-inference blind spot of §3.6** — covered above as exception 1;
   listing it here because it is a miss in the audit's *argument* (structural
   inertness proven one level too deep), not just in its recommendation.

I looked for further uncaught redundancy (T5's spherical loop vs T9's spherical
row, T4's lines 638–639, the api.jl block) and found none worth acting on.

---

## 5. Revised final shape

Same as the audit's §6 except three rows:

| testset | audit | revised | why |
|---|---:|---:|---|
| T1 sweep | 301 | **307** | hoist gets `(A, Cedge)` as a third fixture (miss 2) |
| T2 return type | 21 | **27** | keep both empty-route `typeof`s (exception 2) |
| T7 inference | 52 | **56** | keep one four-op inference row (exception 1) |
| all other rows | — | unchanged | as audited |
| **total** | **426** | **~436** | |

The audit's three do-it-regardless edits (fix T3, delete the tautology at 762,
fix the line-512 comment) are all confirmed and worth doing first — plus a
fourth of the same kind: T3's replacement probe (a) and gate (b) are verified
working exactly as proposed, so the T3 rewrite carries no implementation risk.
