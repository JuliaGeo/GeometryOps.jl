# Audit: the `target` tests for OverlayNG

Scope: the tests added or modified by

    f70da2956  Add a `target` keyword to the OverlayNG operations
    f9de84a32  Assert that a targeted overlay result is type stable
    91531d75e  Give a mixed-dimension overlay result a typed component vector

i.e. section 7 of `test/methods/clipping/overlayng/overlay_ng.jl` (lines 464–785)
and the `"error paths"` testset in `test/methods/clipping/overlayng/api.jl`
(lines 185–191). Sections 1–6 of `overlay_ng.jl` were read but are not in scope
for removal.

Everything below was verified by running code, not by reading test names. The
scratch scripts live in
`/private/tmp/claude-501/-Users-anshul-temp-GO-jts/36b4dc87-418b-4e8c-8a1f-b2b36a1ca475/scratchpad/`
(`exp1.jl`…`exp6.jl`, `sec7.jl`). No source or test file was modified;
mutations were applied at runtime with `Core.eval(GO, …)` and restored
immediately afterwards.

---

## 1. Summary

**926 assertions today** in the nine new testsets (measured — `sec7.jl` run
under `@testset verbose=true` reports exactly 926/926 in 4.7 s), plus one net
new assertion in `api.jl`. **I would keep about 425.**

The three findings that matter:

1. **`"target above the result dimension short-circuits before noding"` does not
   test the short circuit.** All ten of its assertions still pass with
   `_target_above_result` forced to `false` (experiment 17). Every assertion in
   it establishes only "the result is empty", which is true with or without the
   short circuit, because the full pipeline reaches the same empty answer. There
   *is* a decisive, cheap, public-level probe — an input the noder throws on —
   and it is one line (§3.3).

2. **`@test isconcretetype(typeof(gc))` (line 763) is a tautology, and it is the
   assertion aimed squarely at the regression commit 91531d75e was written to
   prevent.** `typeof(x)` always returns a concrete type in Julia. I re-ran T8
   with `_create_result_geometry` reverted to `Vector{Any}` components: line 763
   still passes, as do lines 758 and 767. Only line 761
   (`eltype(GI.getgeom(gc)) === GO._ResultComponent`) and line 753 in T7
   (`all(isconcretetype, ms)`) catch it (experiment 20).

3. **The two big cross products are asymmetric in how much they earn.** The
   625-assertion sweep is *mostly* justified — 12 of its 13 input pairs reach
   genuinely distinct branch signatures, and a designed 8-mutation battery is
   caught by different subsets of them (experiments 10, 15). Only one pair,
   `"overlap+boundary"`, is a provable duplicate, and its explanatory comment
   (line 512) is factually wrong. The 208-assertion inference testset is *not*
   justified: its op axis (×4) and manifold axis (×2) cannot affect the answer
   and are demonstrably inert, and 3 of its 4 input shapes are the same
   const-folded path (experiments 4, 7, 19). It should shrink by ~75% while
   *gaining* a case it currently misses.

Two secondary findings: line 666 in T5 is a verbatim duplicate of the
pre-existing line 391; and the `symdifference` row of T9 (lines 776–777, 780–784)
cannot distinguish `target` being honoured from `target` being ignored, because
the untargeted symmetric difference of that fixture is already a `MultiPolygon`.

---

## 2. Inventory

Assertion counts are measured, not estimated.

| # | testset (line) | asserts | property it establishes | verdict |
|---|---|---|---|---|
| T1 | `target is exactly the untargeted result, filtered by dimension` (507) | 625 | For 13 input pairs × 4 ops: the targeted result's component coordinates equal the untargeted result filtered to that dimension; and the container shape is `Vector` for a singular target / the right `Multi` for a multi target | **TRIM** 625 → ~300 |
| T2 | `target fixes the return type, empty or not` (549) | 27 | An empty targeted result has the same concrete type as a non-empty one at the same dimension, for both empty routes (pipeline-empty and disjoint-envelope short circuit) | **TRIM** 27 → 21 |
| T3 | `target above the result dimension short-circuits before noding` (583) | 10 | *Claims* the pre-noding short circuit. *Actually* establishes only that an above-dimension target returns the empty container, and (line 609) that the short circuit does not misfire at or below the result dimension | **TRIM + FIX** 10 → 7 |
| T4 | `target elides work it cannot want` (612) | 11 | The mixed-point path skips the locator/point-map for a non-point target and the `tuples` copy for a point target; the elision is visible in allocations end to end; the areal answer is unchanged | **TRIM** 11 → 7 |
| T5 | `target on the sphere, and the full-sphere gate` (650) | 12 | Targeting is manifold-agnostic on `Spherical()`; the full-sphere rejection fires for `nothing` and areal targets and is skipped for line/point targets | **TRIM** 12 → 11 |
| T6 | `target accepts the Foster–Hormann spellings and rejects the rest` (678) | 16 | `_ov_target` normalisation: trait instance / trait type / `TraitTarget` / `TraitTarget{Type}` all agree; seven bad spellings raise `ArgumentError` | **KEEP** all 16 |
| T7 | `target makes the return type concrete and input-independent` (715) | 208 | A targeted call infers to one concrete type, the same for every op/manifold/input shape, equal to the runtime type; untargeted infers to a 7-member `Union` all of whose members are concrete | **TRIM** 208 → ~52 |
| T8 | `a mixed-dimension result is a typed collection, not Vector{Any}` (756) | 6 | The `GeometryCollection` branch builds `Vector{_ResultComponent}`, and the fixture really is mixed-dimension | **TRIM** 6 → 4 |
| T9 | `target on all four public operations, including symdifference` (770) | 11 | Each of the four public wrappers forwards `target`; `symdifference`'s bare and manifold-only forms forward it too | **TRIM + FIX** 11 → 7 |
| — | `api.jl` `"error paths"` (185–191) | +1 net | `target` is accepted by `intersection`/`symdifference` under `OverlayNG`; `T` and `fix_multipoly` are not | **KEEP** |

---

## 3. Redundancy findings

### 3.1 T1: `"overlap+boundary"` is a duplicate of `"area × area"`, and its comment is wrong

**Claim.** The pair at line 520–521, `("overlap+boundary", A, Dmixed)` with
`Dmixed = giwkt("POLYGON ((1 0, 4 0, 4 2, 1 2, 1 0))")`, reaches exactly the
same branches as `("area × area", A, B)`, and its comment at line 512 —
`#-- overlaps AND shares a boundary segment: a mixed-dimension result` — is
false.

**Evidence.** Experiment 1 prints the untargeted result trait for every
(pair, op):

```
area x area         int=Polygon(1)  uni=Polygon(1)  dif=Polygon(1)  sym=MultiPolygon(2)
overlap+boundary    int=Polygon(1)  uni=Polygon(1)  dif=Polygon(1)  sym=MultiPolygon(2)
```

The intersection is a plain `Polygon`, not a mixed-dimension collection — which
is exactly what the *other* comment in the same block already says (lines
498–500: "A single polygon that overlaps and shares boundary yields just the
polygon"). Experiment 10 computes, for every pair, the vector of component
counts the sweep actually compares (`path | count(dim2), count(dim1), count(dim0)`
for each of the four ops):

```
area x area         AR|1,0,0,1,0,0,1,0,0,2,0,0
overlap+boundary    AR|1,0,0,1,0,0,1,0,0,2,0,0
...
  distinct signatures: 12 of 13 pairs
  DUPLICATE signature AR|1,0,0,1,0,0,1,0,0,2,0,0 : ["area x area", "overlap+boundary"]
```

Identical driver path (`AR` = arrangement/`_extract_result`), identical
emptiness pattern at every dimension for every op. The target machinery branches
only on (driver path, op, target dimension, emptiness of the component list), so
identical signature means identical branch coverage. The 8-mutation battery
(experiment 15) confirms it: the two pairs are caught by, and only by, exactly
the same mutations.

**Recommended edit.** Delete `("overlap+boundary", A, Dmixed)` from `pairs`
(line 520), delete the now-unused `Dmixed` binding (line 513) and its wrong
comment (line 512). Saves 48 assertions. (Collinear-shared-boundary *noding* is
already covered by the pre-existing testset at line 90,
`"collinear shared boundary (degenerate intersection, merged union)"`.)

**Confidence: high.**

### 3.2 T1: the two container-shape assertions do not belong inside the pair × op loop

**Claim.** Lines 543–544

```julia
@test GI.trait(mul) === multi
@test vec isa AbstractVector
```

are 312 of T1's 624 loop assertions (156 each). Both establish a property of
`_target_result`'s six trait methods, not of the input pair: the shape is a
function of the target alone. One instance per target — plus one empty and one
non-empty instance for the multi targets, whose method has two branches
(`isempty(v) ? _empty_geom(d) : GI.MultiXxx(v)`, source lines 333–338) — covers
everything the 312 do.

**Evidence.** Experiment 6 shows the returned types are fixed by the target.
Experiment 4 shows inference proves the return type is a single concrete type
per target, independent of every other axis. The remaining question is whether
running them 156 times catches something one run does not; it does not, because
the only inputs to `_target_result` that vary are the three component vectors,
and the shape branch reads only `isempty`.

**Recommended edit.** Remove lines 543–544 from the inner loop and add, after
the loop, a small block:

```julia
#-- the container shape is a function of the target alone: `Vector` for a
#-- singular trait, the right `Multi` for a multi trait, empty or not
for (single, multi, _) in TARGETS
    for (X, Y) in ((A, B), (A, Etouch))   # non-empty and empty at most dims
        @test GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, X, Y;
                             exact = EX, target = single) isa AbstractVector
        @test GI.trait(GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, X, Y;
                                      exact = EX, target = multi)) === multi
    end
end
```

12 assertions in place of 312. Combined with §3.1 this takes T1 from 625 to
`12 pairs × 4 ops × 3 targets × 2 + 12 + 1 guard = 301`.

**Do not simply delete these two assertions** — see §5.1; they are the only
absolute pinning of the container shape for four of the six targets.

**Confidence: high.**

### 3.3 T3 establishes emptiness, not the short circuit

**Claim.** None of T3's ten assertions can fail when the short circuit is
removed. The testset's stated property is untested.

**Evidence.** Experiment 17 replays all ten of T3's assertions, then redefines

```julia
GO._target_above_result(t, op::GO._OverlayOpCode, dim_a, dim_b) = false
```

and replays them again:

```
  baseline T3 failures: 0
  with `_target_above_result` DISABLED, T3 failures: 0
  restored T3 failures: 0
```

Why: with the short circuit gone, the full pipeline runs, produces no components
at the target's dimension, and `_target_is_empty` → `_resolve_empty_result` →
`_empty_target_result` returns the very same empty container. The observable
answer is identical. The comment at lines 604–605 — "the short circuit reads the
input dimensions only, so it survives inputs the pipeline would reject outright"
— names the right idea but picks an input that the pipeline does *not* reject:
experiment 5 shows `GO.intersection(GO.OverlayNG(), GI.LineString([(0.0,0.0),(1.0,1.0)]), A)`
succeeds untargeted and returns a `LineString`.

**Recommended edit.** Two replacements, both at the public level.

(a) An input the noder genuinely rejects. Antipodal spherical vertices throw
inside noding; with an above-dimension target the call must return empty
*instead of throwing*, which is only possible if noding never ran
(experiment 22):

```julia
#-- the short circuit reads the input dimensions only, so it answers inputs the
#-- noder would throw on: an antipodal spherical edge cannot be noded at all
Ant = GI.LineString([(0.0, 0.0), (180.0, 0.0)])
SA  = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
sph = GO.OverlayNG(Spherical())
@test_throws ArgumentError GO.intersection(sph, Ant, SA)
@test isempty(GO.intersection(sph, Ant, SA; target = GI.PolygonTrait()))
```

Measured: untargeted throws `ArgumentError: spherical edge between antipodal
vertices (1.0, 0.0, 0.0) and (-1…`; targeted returns `Vector n=0`.

(b) An allocation gate, in the same spirit as T4's (experiment 21, planar,
warm):

| call | allocations |
|---|---|
| `intersection(alg, L, A; target = GI.PolygonTrait())` (short-circuits) | 128 B |
| `intersection(alg, L, A; target = GI.LineStringTrait())` (full pipeline) | 21 504 B |
| same, with `_target_above_result` disabled | 1 779 456 B |

```julia
GO.intersection(alg, L, A; target = GI.PolygonTrait())   # warm up
@test @allocated(GO.intersection(alg, L, A; target = GI.PolygonTrait())) < 1024
```

Then keep four of the eight `unsat` rows (one per op — the row's op is the only
axis `_result_dimension` reads) and keep the negative control at line 609.
Net: 4 + 2 + 1 = 7 assertions, and for the first time they test the named
property.

**Confidence: high** for the finding; **high** for probe (a), which I ran.

### 3.4 T4's last loop is subsumed by T1's `"area × area"` pair

**Claim.** Lines 641–647,

```julia
for op in OPS
    plain = GO._overlay_ng(Planar(), op, A, B; exact = EX)
    tgt = GO._overlay_ng(Planar(), op, A, B; exact = EX, target = GI.MultiPolygonTrait())
    @test isapprox(GO.area(tgt), GO.area(plain); rtol = 1e-12)
end
```

is a strictly weaker restatement of what T1's `("area × area", A, B)` row already
checks. T1 compares `GI.coordinates` of every component of the targeted result
against the untargeted result filtered to dimension 2, for all four ops and for
both `GI.PolygonTrait()` and `GI.MultiPolygonTrait()` — coordinate equality
implies area equality. The `A` and `B` fixtures are identical geometries in both
testsets (lines 613/642 vs 508/509).

**Recommended edit.** Delete lines 641–647 (and the now-unused `B` at line 642).
Saves 4 assertions.

**Confidence: high.**

### 3.5 T5 line 666 is a verbatim duplicate of the pre-existing line 391

**Claim.** Line 665–666,

```julia
inp = GO._OverlayInput(Spherical(), A, A, 2, 2, EX, false, false, nothing, nothing)
@test_throws ArgumentError GO._resolve_empty_result(Spherical(), GO.OVERLAY_UNION, inp)
```

is character-for-character the same assertion as the pre-existing line 389–391 in
`"spherical empty-vs-full disambiguation (§3 amendment 6)"`, on the same `A`
(both are `GI.Polygon([[(0.0,0.0),(10.0,0.0),(10.0,10.0),(0.0,10.0),(0.0,0.0)]])`,
lines 380 and 651).

**Recommended edit.** Delete line 666. Its role as the "target = nothing still
throws" control is already served by lines 667–670, which assert the same throw
for both areal targets in the same block, immediately above the non-throwing
line/point cases. Saves 1 assertion.

**Confidence: high** on the duplication; **medium** on removing it — one could
argue the untargeted control belongs beside the targeted ones for readability.
If you keep it, add `# (also asserted at line 391)`.

### 3.6 T7: the op axis and the manifold axis are inert; the input-shape axis collapses 4 → 2 and is missing a third case

**Claim.** T7 runs `6 targets × 2 manifolds × 4 ops × 4 input shapes = 192`
`isconcretetype` assertions (plus 16 more) where roughly 36 would give strictly
more coverage.

**Evidence.**

*The result is structurally independent of op and manifold.* Every `return` site
in `_overlay_ng` funnels into `_target_result(target, …)` or
`_empty_target_result(target)` (source lines 113, 116, 121, 123, 133 →
`_dimensional_result`/`_extract_result`). `target` is a singleton, so those
resolve to one method with one concrete return type. `op` is only ever compared
to a constant (`build_points`, `_result_dimension`), never dispatched on;
`Spherical`'s only extra path is a `throw`, which contributes `Union{}`.
Experiment 4 confirms: for each of the six targets, all 32 cells infer to
**exactly one** type.

*The input-shape axis is real, but not four-wide.* Inference **does**
const-fold the driver's dimension branches, so different argument types reach
different return-site sets (experiment 7, untargeted, so the union size is
visible):

```
  untargeted Poly,Poly      -> union of 7 members
  untargeted Poly,MPoly     -> union of 7 members
  untargeted Line,Poly      -> union of 7 members
  untargeted MPoint,Poly    -> Any            (the mixed-point path)
  untargeted MPoint,MPoint  -> Any            (the point×point path)
```

The three argtypes at lines 721–722 that are *not* `MultiPoint`-first are the
same const-folded shape. Experiment 19 makes this concrete with a regression
that only affects the mixed-point path (`_typed`'s fallback method returning
`Any[…]`, which only `_mixed_components` reaches):

```
  regression A (shared assembler untyped): non-concrete cells = 32/32
  regression B (`_typed` fallback untyped): non-concrete cells = 8/32
    caught in: ["MPoint,Poly"]
```

Regression A is caught by every cell — one would do. Regression B is caught by
all 8 `MPoint,Poly` cells (2 manifolds × 4 ops) and by none of the other 24 —
so the argtype axis earns its keep, and the manifold and op axes within it do
not. And `Tuple{MultiPoint, MultiPoint}` — the `_overlay_points` path, a third
distinct const-folded shape with its own `_dimensional_result` return site — is
**not probed at all**.

**Recommended edit.** Replace lines 721–722 and 726–733 with:

```julia
#-- three input shapes, because inference const-folds the driver's dimension
#-- branches: these are the three distinct return-site sets (arrangement,
#-- mixed-point, point×point). The op and manifold axes cannot move the answer —
#-- every return site funnels through `_target_result(target, …)`.
argtypes = [Tuple{typeof(A), typeof(B)}, Tuple{typeof(P), typeof(A)},
            Tuple{typeof(P), typeof(P)}]

for (single, multi, _) in TARGETS, t in (single, multi)
    inferred = Type[]
    for m in (Planar(), Spherical()), Ts in argtypes
        probe = (a, b) -> GO.intersection(GO.OverlayNG(m), a, b; target = t)
        R = Base.return_types(probe, Ts)[1]
        @test isconcretetype(R)
        push!(inferred, R)
    end
    @test length(unique(inferred)) == 1
    @test only(unique(inferred)) === typeof(GO.intersection(GO.OverlayNG(), A, B; target = t))
end
```

`6 × 8 = 48` plus the 4 trailing assertions = 52, down from 208, with the
point×point path newly covered. Delete the now-unused `MB` (line 718) and `L`
(line 719) — or keep `L` if you prefer a third arrangement-shaped case; it adds
nothing measurable.

Also worth fixing while there: the block comment at lines 704–708 says
"Untargeted, the return type is a 7-member `Union`". That is true only for
inputs with edges; for `MultiPoint`-first inputs it is `Any` (experiment 12).
Either narrow the sentence or make it a claim about areal inputs.

**Confidence: high** on op/manifold inertness (structural argument plus
measurement); **high** on the argtype collapse (measured).

### 3.7 T8: one tautology, one subsumed assertion, one duplicate of T1's guard

**Claim and evidence.** Experiment 20 reverts `_create_result_geometry` to build
`Vector{Any}` components — the exact regression commit 91531d75e prevents — and
re-runs each of T8's four distinct assertions:

| line | assertion | baseline | with `Vector{Any}` |
|---|---|---|---|
| 758 | `GI.trait(gc) isa GI.GeometryCollectionTrait` | pass | **pass** (does not catch) |
| 761 | `eltype(GI.getgeom(gc)) === GO._ResultComponent` | pass | **fail** (catches) |
| 763 | `isconcretetype(typeof(gc))` | pass | **pass** (does not catch) |
| 767 | `all(g -> g isa GO._ResultComponent, GI.getgeom(gc))` | pass | **pass** (does not catch) |

- **Line 763 is a tautology.** `typeof(x)` returns a concrete `DataType` for
  every value in Julia. Experiment 3 shows `isconcretetype(typeof(…))` is `true`
  even for `GI.GeometryCollection{false,false,Vector{Any},Nothing,Nothing}`. It
  cannot fail, and in particular cannot fail for the regression it sits next to.
  The property the comment wants — "the wrapper type infers concretely" — is a
  statement about `Base.return_types`, and it *is* tested, at line 753
  (`@test all(isconcretetype, ms)`), which experiment 20 confirms does catch the
  regression (`all members concrete? false`).
- **Line 767 is subsumed by line 761.** If `eltype(GI.getgeom(gc))` is exactly
  `_ResultComponent`, every element of that container is of a type in the union
  by the container's own invariant, so `all(g -> g isa _ResultComponent, …)`
  cannot fail while 761 passes.
- **Line 758 duplicates line 530**, T1's fixture guard, which asserts
  `GI.trait(GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, MIXED_A, MIXED_B; exact = EX))
  isa GI.GeometryCollectionTrait` on the same fixture. (`GO.intersection(GO.OverlayNG(), …)`
  is a one-line pass-through to `_overlay_ng`, api.jl line 232.) Keep one; I'd
  keep line 530, because it guards the sweep's coverage of the flattening branch,
  and keep 758 too only as readable setup — it is one assertion.

**Recommended edit.** Delete lines 763 and 767. Keep 758, 761, 765, 766. 6 → 4.
If you want line 763's *intent* preserved, replace it with a comment pointing at
line 753, since that is where the property lives.

**Confidence: high** (both the tautology and the mutation experiment).

### 3.8 T9's `symdifference` rows cannot distinguish `target` honoured from `target` ignored

**Claim.** For the fixture `A`, `B` used in T9, the *untargeted* symmetric
difference is already a `MultiPolygon`, so
`@test GI.trait(r) isa GI.MultiPolygonTrait` and
`@test isapprox(GO.area(r), GO.area(f(alg, A, B)))` both pass even if the
`target` keyword were silently dropped. This affects lines 776–777 for
`f = GO.symdifference`, and lines 780, 781, 782–784 in their entirety.

**Evidence.** Experiment 24:

```
  intersection   untargeted trait = PolygonTrait()       passes `isa MultiPolygonTrait`? false  areas equal? true
  union          untargeted trait = PolygonTrait()       passes `isa MultiPolygonTrait`? false  areas equal? true
  difference     untargeted trait = PolygonTrait()       passes `isa MultiPolygonTrait`? false  areas equal? true
  symdifference  untargeted trait = MultiPolygonTrait()  passes `isa MultiPolygonTrait`? true   areas equal? true
  symdifference(A, B) bare form untargeted trait      = MultiPolygonTrait()
  symdifference(Planar(), A, B) untargeted trait      = MultiPolygonTrait()
```

The area assertion is non-discriminating for *all four* ops (an ignored target
returns the same geometry, hence the same area) — it is dead weight everywhere,
not just for symdifference.

**Recommended edit.** Switch to the singular trait, whose return shape
(`Vector`) no untargeted call can produce, and drop the area assertions:

```julia
for f in (GO.intersection, GO.union, GO.difference, GO.symdifference)
    v = f(alg, A, B; target = GI.PolygonTrait())
    @test v isa AbstractVector
    @test all(g -> GI.trait(g) isa GI.PolygonTrait, v)
end
#-- symdifference's algorithm-free and manifold forms take it too
@test GO.symdifference(A, B; target = GI.PolygonTrait()) isa AbstractVector
@test GO.symdifference(Planar(), A, B; target = GI.PolygonTrait()) isa AbstractVector
@test GO.symdifference(Spherical(), A, B; target = GI.PolygonTrait()) isa AbstractVector
```

11 → 11 with the same shape, or 7 if you drop the `all(…)` row. Note this
overlaps `api.jl` lines 188–189, which already assert exactly this for
`intersection` and `symdifference` under `OverlayNG()` — a deliberate
cross-file duplicate that I would keep, since `api.jl`'s job is the public
surface.

**Confidence: high.**

### 3.9 T2 partially overlaps T7

**Claim.** T2's `@test typeof(e) === typeof(nonempty)` assertions (lines 568,
578) are implied by T7 for the dimension-2 and dimension-1 rows.

**Evidence.** T7 proves `Base.return_types` is a *concrete* type for
`Tuple{typeof(A), typeof(B)}` and a given target. A concrete inferred return
type means every actual return value has exactly that type. T2's `empties`
(`GO.intersection(alg, A, Far; …)` and `GO.difference(alg, A, A; …)`, lines
563–564) and its dimension-2 / dimension-1 `nonempty` values are all
`Tuple{Polygon, Polygon}` calls with the same concrete argument types T7 probes
— `A`, `B`, `Cedge`, `Far` are all
`GI.Polygon{false,false,Vector{Vector{Tuple{Float64,Float64}}},Nothing,Nothing}`.
So the `typeof` equality follows. The dimension-0 `nonempty`
(`MultiPoint × MultiPoint`, lines 561–562, 573–574) is *not* covered, because
T7 does not probe that argtype today — which §3.6 fixes.

What T2 uniquely establishes and T7 does not is the **value** property: that
these particular calls really do come back empty, and that a non-empty one
exists at each dimension. That is worth keeping — an engine that returned the
right *type* but the wrong *contents* would pass T7 and fail T2.

**Recommended edit.** Keep the `isempty` / `GI.ngeom(e) == 0` /
`!isempty(nonempty)` assertions in full (both empty routes matter — one is the
disjoint-envelope short circuit, the other is the pipeline). Reduce the
`typeof` assertions from two per container shape to one:

per dimension: `!isempty(nonempty)` (1) + `isempty(e)` × 2 (2) + one
`typeof(e) === typeof(nonempty)` (1) + `GI.ngeom(e) == 0` × 2 (2) + one
`typeof(e) === typeof(nonempty_m)` (1) = 7, × 3 = **21**, down from 27.

**Confidence: medium-high.** The inference→runtime implication is sound, but it
is an indirect argument; if you prefer a testset that stands on its own without
depending on T7, leave T2 at 27. It costs 0.0 s.

### 3.10 Non-findings, stated so they are not re-litigated

Three things I checked and found **not** redundant, listed here because they
look redundant:

- **`_target_is_empty`'s areal method vs the untargeted rule.** Rewriting
  `_target_is_empty(::Union{PolygonTrait, MultiPolygonTrait}, …)` to the
  untargeted `isempty(polys) && isempty(lines) && isempty(points)` is caught by
  0/13 sweep pairs (experiment 14) — but that is because the mutation is
  behaviour-preserving on the planar path (both routes end at the same empty
  vector). It is *not* behaviour-preserving on the sphere, where it changes
  whether `_resolve_empty_result`'s full-sphere gate is consulted. T5 is what
  covers that.
- **My `_typed`-drops-a-component mutation caught 0/13** — that is a flaw in the
  mutation, not in the sweep: I redefined only the fallback
  `_typed(::Type{T}, v)`, and the component vectors already hit the
  `_typed(::Type{T}, v::Vector{T})` fast path. Disregard that row of
  experiment 15.
- **T6 has no redundancy.** Its four accepted spellings exercise four distinct
  `_ov_target` methods (source lines 273–280), and its seven rejected values
  exercise both rejection sites (the `isconcretetype` guard and the catch-all).
  Keep all 16.

---

## 4. Interface vs internals

I verified the list of internals reached by the new tests rather than trusting
it. The complete list, with line numbers:

| internal | reached at | is a public-level equivalent available? | recommendation |
|---|---|---|---|
| `GO._overlay_ng` | 530, 534, 537–538, 644–645, 656–657 | Yes (`GO.intersection(GO.OverlayNG(), …)` etc.), but this is the file's established convention — sections 1–6 use it throughout (lines 38, 69–72, 255, 271, …) and the public wrappers are one-line pass-throughs (api.jl 232–239, 274) | **KEEP.** Not a finding. |
| `GO.OVERLAY_INTERSECTION` / `_UNION` / `_DIFFERENCE` / `_SYMDIFFERENCE` | directly at 530, 666, 668, 673; via the pre-existing `OPS` const (line 24) at 533, 643, 655 | Same convention; they are the argument `_overlay_ng` takes | **KEEP.** Not a finding. |
| `GO._mixed_points`, `GO._NO_POINTS` | 624, 625 | Partly. Experiment 18 shows the public allocation assertion at 634 *does* catch a global elision regression on its own (targeted 1 720 096 B vs untargeted 1 953 392 B — ratio 1.14, fails the ÷10 gate). But the assertion at 624 pins *which* structure is elided and that the shared constant is returned by identity, which no allocation count can | **KEEP 624.** It is the only test of the shared-empty-list identity contract (source lines 344–346, "Shared and never mutated"). **Consider dropping 625** — the negative control passes for any freshly-allocated vector and is the weakest of the four. |
| `GO._mixed_components`, `GO._NO_COMPONENTS` | 626, 627, 629 | No. The elision this pins — a *point* target skipping the `tuples` deep copy of the non-point input — has no public probe today: the allocation test at 634 uses an areal target and so exercises the other half | **KEEP 626–627.** See §5.3. **Consider dropping 629** (same negative-control weakness as 625). |
| `GO._resolve_empty_result`, `GO._OverlayInput` | 665–675 (and pre-existing 389–396) | No, as far as I could establish. See §5.2 | **KEEP,** minus the duplicate line 666 (§3.5). |
| `GO._ResultComponent` | 761, 767 | Only partly: `GI.getgeom(gc)`'s `eltype` is observable publicly, but the *name* of the union is not. Line 761 is the sole catcher of the `Vector{Any}` regression (§3.7) | **KEEP 761.** Drop 767 (subsumed). |
| `GO._empty_geom` | **not reached by any new test** — the list in the brief is wrong here (verified by scanning lines 464–785 for `GO._empty_geom`) | — | n/a |
| `GO.TraitTarget` | 684, 696–697 | `TraitTarget` is exported and documented (AGENTS.md, `TraitTarget{GI.PointTrait}()`), so this is public surface | **KEEP.** Not an internal. |
| `GO._covers_everything` | **not reached by any new test** (only by the pre-existing line 390/394) | — | n/a |

**Overall verdict on interface-vs-internals: the new tests are in good shape.**
Nothing new reaches into the graph, the labeller or the builders. The four
genuine deep reaches — `_mixed_points`/`_mixed_components` (T4) and
`_resolve_empty_result`/`_OverlayInput` (T5) — are all testing *elisions* and
*guards*, i.e. things whose defining property is that they change no answer.
That is precisely the category where an internal test is the right instrument.

The one place I would push toward the public level is T3, and there the fix is
not "reach for an internal" but "pick a better public input" (§3.3).

---

## 5. Do not delete

### 5.1 T1's container-shape assertions (lines 543–544) — reduce the count, keep the property

§3.2 recommends cutting these from 312 instances to 12. Do **not** cut them to
zero. `@test GI.trait(mul) === multi` is the sole catcher of a wrong-dimension
empty container: experiment 23 replaces `_empty_geom` with one that always
returns a `MultiLineString`, and

```
  baseline: ["coordlist-equal=true", "trait===MultiPointTrait=true"]
  with `_empty_geom` always returning a MultiLineString:
            ["coordlist-equal=true", "trait===MultiPointTrait=false"]
```

The coordinate-list comparison passes, because both lists are empty. Only the
trait assertion fails. The same argument applies to `vec isa AbstractVector`:
if `_target_result(::GI.PolygonTrait, …)` returned a `MultiPolygon`, the sweep's
`targeted_parts` helper (line 487) would silently fall through to
`collect(GI.getgeom(r))` and the coordinate comparison would still pass. These
two are the only *absolute* pinning of the container shape for the
`LineString` / `MultiLineString` / `Point` / `MultiPoint` targets — T2 and T7
both compare types *relative* to another result, and `api.jl` 188–189 only pins
the polygon targets.

### 5.2 T5's `_resolve_empty_result` / `_OverlayInput` block (lines 665–675)

This is the reverse problem: the property is only testable at the internal
level, and deleting the internal test leaves a real gap.

I tried hard to reach `_FULL_SPHERE_MSG` from the public API. Experiments 8 and
13 build two complementary "hemispheres" as spherical polygons whose ring is a
great circle and union them:

```
  union target=MultiPolygonTrait -> ArgumentError
  union target=MultiLineStringTrait -> ok, ngeom=4
  union target=MultiPointTrait -> ok, ngeom=0
```

which *looks* like exactly the property T5 asserts — but the message is
`rk_point_in_ring: every anchor edge of the ring is degenerate with respect to
the query point`, not `_FULL_SPHERE_MSG`. A great-circle ring has no enclosed
region under `Spherical(; oriented = false)`, so the fixture is degenerate; the
throw comes from the locator inside `_covers_everything`, not from the
full-sphere gate. It would be a bad test: it asserts an unrelated failure.

The source itself explains why (overlay_ng.jl lines 425–429): reaching
`_FULL_SPHERE_MSG` "requires an overlay whose result has *no* boundary at all
and covers everything: in practice a union/symdifference of complementary
hemispheres, or a union whose operands' boundaries cancel exactly" — and the
representable spelling of "complementary hemispheres" is exactly what the model
lacks. Constructing `_OverlayInput` by hand is the only way to put the resolver
in that state. Lines 665, 667–675 stay.

(If someone later finds a genuinely reachable end-to-end full-sphere case, that
public test should be *added*, not substituted — the internal one pins the
`_target_admits_area` branch directly, which is the new behaviour.)

### 5.3 T4's `_mixed_components(A, 2, mpt)` assertion (lines 626–628)

The point-target half of the mixed-point elision — skipping the `tuples` deep
copy of the non-point input — has no public probe in the block. T4's allocation
gate at line 634 uses `target = mpoly`, which exercises the *other* half (the
locator/point-map elision). Deleting lines 626–628 would leave
`_mixed_components`' `_target_needs_dim` guard (overlay_mixed_points.jl lines
83–87) untested. Experiment 18 confirms line 626–627 is one of the three
assertions that catch a global elision regression.

If you would rather have this at the public level, the replacement is a second
allocation gate:

```julia
BigA = ...          # a polygon with many vertices, so the `tuples` copy dominates
@test @allocated(GO.union(alg, Pts, BigA; target = mpt)) <
      @allocated(GO.union(alg, Pts, BigA)) ÷ 2
```

I did not measure this one; treat it as a suggestion, not a verified
replacement.

### 5.4 T4's allocation gate (lines 632–635)

It looks fragile (a `÷ 10` threshold on `@allocated`) and it is not: measured
ratio is **3815×** (512 B targeted vs 1 953 392 B untargeted, experiment 9), and
under the elision regression it drops to 1.14× (experiment 18). The margin is
three orders of magnitude. Allocation counts are also unaffected by
`--code-coverage`, which is why this is the right instrument here rather than
wall clock — the comment at lines 631–632 says so and is correct.

### 5.5 The `"corner touch"`, `"point × point"` and `"point × line"` sweep pairs

The designed mutation battery (experiment 15) found pairs that are the *sole*
catcher of a path-specific regression:

```
  point x point path ignores its target                1/13 : ["point x point"]
  _mixed_components always reports the non-point side as areal
                                                       1/13 : ["point x line"]
  _extract_result never builds points (over-elision)   2/13 : ["corner touch", "line x line"]
```

`("point × point", P, P2)` is the only pair reaching `_overlay_points`;
`("point × line", P, L)` is the only one with a *line* non-point side in
`_mixed_components`; `("corner touch", A, Etouch)` is the only areal pair whose
intersection is a point, i.e. the only areal pair reaching `_build_points`.
None of these three may be dropped.

For the record, a greedy cover of all 8 caught mutations is
`["line × line", "point × area", "point × line", "point × point"]` — four pairs.
I am **not** recommending that cut. The mutation battery is my construction and
is not exhaustive; the branch-signature analysis (12 distinct of 13) is the
better guide, and the whole sweep costs 2 ms warm / 3.9 s cold, of which almost
all is compilation that trimming pairs would not recover (the 13 pairs use only
about ten distinct argument-type combinations, and removing `"overlap+boundary"`
removes none of them).

---

## 6. Recommended final shape

| testset | now | after | what changed |
|---|---:|---:|---|
| `target is exactly the untargeted result, filtered by dimension` | 625 | **301** | drop the `"overlap+boundary"` pair and its wrong comment (§3.1); hoist `vec isa AbstractVector` / `GI.trait(mul) === multi` out of the pair × op loop into a 12-assertion per-target block (§3.2) |
| `target fixes the return type, empty or not` | 27 | **21** | one `typeof` assertion per container shape instead of two; all value assertions kept (§3.9) |
| `target above the result dimension short-circuits before noding` | 10 | **7** | four `unsat` rows (one per op) instead of eight; the non-probe at 606–607 replaced by the antipodal-spherical probe; add an `@allocated` gate; keep the negative control (§3.3) |
| `target elides work it cannot want` | 11 | **7** | delete the four-op areal-area loop at 641–647 (§3.4); drop the two weak negative controls at 625 and 629 (§4) |
| `target on the sphere, and the full-sphere gate` | 12 | **11** | delete the duplicate at 666 (§3.5) |
| `target accepts the Foster–Hormann spellings and rejects the rest` | 16 | **16** | unchanged (§3.10) |
| `target makes the return type concrete and input-independent` | 208 | **52** | drop the op axis (inert) and shrink the argtype axis 4 → 3, adding the missing `MultiPoint × MultiPoint` path; keep both manifolds; fix the block comment's "7-member Union" claim (§3.6) |
| `a mixed-dimension result is a typed collection, not Vector{Any}` | 6 | **4** | delete the tautology at 763 and the subsumed assertion at 767 (§3.7) |
| `target on all four public operations, including symdifference` | 11 | **7** | target `GI.PolygonTrait()` and assert `isa AbstractVector` so the assertions actually discriminate; drop the non-discriminating area comparisons (§3.8) |
| **total** | **926** | **426** | |

`api.jl`'s `"error paths"` edit (lines 185–191) stays exactly as written. I
verified all four of its assertions behave as claimed (experiment 25):
`intersection`/`symdifference` with `target = GI.PolygonTrait()` return an
`AbstractVector` where the untargeted calls return a `Polygon` and a
`MultiPolygon` respectively — so both are discriminating — and both the
`Float32` positional and the `fix_multipoly` keyword raise `MethodError`.

Three edits are worth doing even if none of the trimming is:

1. **Fix T3** (§3.3). It currently asserts a property it cannot observe.
2. **Delete line 763** (§3.7). It is a tautology sitting where a reader will
   assume the `Vector{Any}` regression is guarded.
3. **Fix the `"overlap+boundary"` comment at line 512** (§3.1), whether or not
   the pair is removed. It contradicts the comment 12 lines above it and is
   wrong.
