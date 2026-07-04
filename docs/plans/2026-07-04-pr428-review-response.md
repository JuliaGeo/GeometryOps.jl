# PR #428 review response plan (2026-07-04)

27 inline review comments on the prepared-geometry PR.  They cluster into
seven themes; this plan addresses each theme, maps every comment to an
action, and extends the same standards to the parts of the diff the review
didn't reach.

## Theme A — comment doctrine: restore what was cut, fix the register

The review draws one consistent editorial line, from two directions:

1. **The "condense comments" pass (a107d0a7e) went too far.**  Didactic
   comments that teach the reader what the code is doing — the sorted-query
   explanation, the LoopStateMachine control-flow note, the nested-loop
   equivalence block, the pre-allocation note — were load-bearing.  Restore
   them.
2. **New comments must describe behavior, not tell stories.**  Comments that
   narrate the port's history ("adapted from X with two structural
   changes"), cite profiling sessions ("dominated traversal time in
   profiles"), or congratulate the design are the wrong register.  Rewrite
   to state what the code does and what invariants hold; API-shaped items
   get docstrings, not section comments.

Actions:
- Apply the three verbatim suggestions: accelerator docstring bullet list
  (clipping_processor:14, updated to the final accelerator vocabulary from
  Theme B), iteration-order paragraph (:196), `query_result` pre-allocation
  comment (:255), and the `RTree` constructor `where A` change
  (FlexibleRTrees:122).
- Restore the deleted comments at clipping_processor :298 (sorted to mimic
  nested-loop order), :342 (query-sort + LoopStateMachine control flow),
  :566 (the "equivalent to this nested loop" julia code block in
  `_build_a_list`), and sweep `git diff main...HEAD` for every other
  deleted `#` line in clipping_processor.jl, intersection.jl,
  sutherland_hodgman.jl, union.jl — restoring the explanatory ones
  (the Hao–Sun case tables already survived inside `_hao_sun_edge`; the
  stale `l = GI.Line...` TODO and dead-code comments stay gone).
- Rewrite the FlexibleRTrees module header (:22): describe the tree
  (flat level-vectors + leaf permutation, concrete type at any depth,
  algorithm = leaf ordering via `loadorder`), keep a one-line MIT
  attribution to SortTileRecursiveTree.jl, drop the compare-and-contrast
  narrative.
- RTreeNode section (:215): docstring on `RTreeNode` describing the cursor
  scheme (level/index into `levels`, children of level-l node at l+1,
  leaf slots mapped through `indices`); delete the story framing.
- NaturalIndexing :162: replace the profile anecdote with the plain
  invariant ("children of one node share one extents vector; resolve it
  once per node, then index").

## Theme B — an explicit accelerator vocabulary (comment :210)

The `hasprep` branches inside `AutoAccelerator` are opaque: prep-on-a →
`DoubleNaturalTree`, prep-on-b → `SingleNaturalTree`, both-prepared is not
distinguished, and the dual-traversal machinery is never chosen for the
both-prepared case.  Replace implicit knowledge with types that say what
happens on each side.

Design — one parametric accelerator with a per-side policy:

```julia
struct IterateEdges end                 # no tree; walk this side's edges in order
struct BuildTree{B} backend::B end      # build an ephemeral tree over this side
struct ReuseTree{F} fallback::F end     # use this side's prepared tree; fall back on a miss

struct TreeAccelerator{PA, PB} <: IntersectionAccelerator
    a::PA
    b::PB
end
```

- `TreeAccelerator(IterateEdges(), BuildTree(NaturalIndex))` is today's
  `SingleNaturalTree`; `TreeAccelerator(BuildTree(..), BuildTree(..))` is
  the dual traversal; `NestedLoop` stays as-is (it is already the clearest
  name for iterate×iterate).
- Existing singleton names (`SingleSTRtree`, `SingleNaturalTree`,
  `DoubleSTRtree`, `DoubleNaturalTree`) become thin constructors/aliases so
  the exported vocabulary and tests keep working.
- `AutoAccelerator` becomes a documented decision table that resolves to a
  `TreeAccelerator`:
  - both sides prepared → `ReuseTree` × `ReuseTree` (dual traversal —
    the case the review points out we currently drop);
  - one side prepared → `ReuseTree` on that side; other side by size
    (small → `IterateEdges` when it's `a`, else `BuildTree`);
  - neither prepared → current size heuristic (small×small `NestedLoop`,
    one small → single tree on b, big×big → dual build).
- The callback-order contract (in-`eachedge`-order over `a`) is what forces
  the asymmetry (we can never put "iterate" on `b` alone); state that in
  the docstring.

## Theme C — make the parts machinery legible (:296, :302, :366, :384, :418)

Why it exists (the answer the code failed to give): tree queries return
*random* edge indices, and `eachedge` is a sequential iterator — so the
loops need random access to edge coordinates; and a multi-ring geometry's
prepared trees live per-ring, so global edge numbering needs per-ring
offsets.  The per-ring reuse is a measured win (ring-aware accelerator
benchmarks in the v2 design doc) — keep the capability, rebuild the
presentation:

- **One coordinate accessor instead of two.**  The ephemeral path
  materializes the curve to tuple storage (what `prepare` would store)
  instead of a separate `to_edgelist` vector; then a single accessor
  ("edge j = points j, j+1") serves both reused and ephemeral trees, and
  `_EdgeListCoords` is deleted.  Foreign geometries still get materialized
  exactly once (random `getpoint` on e.g. ArchGDAL would cross FFI per
  access, so in-place reads are only safe on native storage).
- **Inline the loops.**  `_single_parts_loop` folds back into the
  single-tree accelerator method; `_dual_pairs_loop` folds into the dual
  method, with restored narrative comments explaining the shape: collect
  candidate (i, j) pairs from the simultaneous traversal, sort them so the
  callbacks fire in the order the nested loop would, then walk.  The only
  extracted helpers that remain are the function barriers that type
  stability genuinely needs, each labeled as such.
- **Domain names.**  `_EdgePart` → `_RingTree` (tree + offset + extent for
  one ring/curve); `_PartsCoords` → `_EdgeCoords` with the offset scan
  commented (geometries hold few rings).
- Fix `_edge_parts`'s `Any[]` + `map(identity, …)` with a properly typed
  build.
- Theme D deletes the closedness checks in `_edge_tree_and_coords`, which
  is most of its remaining bulk.

## Theme D — normalize at `prepare` time; know the manifold (:101, :224, :290, :305)

`prepare` already materializes; make materialization do the normalization
consumers currently re-check, and thread the manifold through so storage
and edge extents are manifold-appropriate:

- **Close rings during materialization.**  A materialized `LinearRing`
  always repeats its first point.  New invariant, documented on
  `Prepared`: *any preparation retrieved from a `Prepared` node was built
  against materialized storage, and materialized rings are closed.*
  Consequences, all deletions:
  - `_edge_tree_and_coords`'s "reusable only when closed" check dies;
  - `_point_ring_orientation`'s closedness check dies — the helper
    reduces to "prep present → indexed walk, else sequential" and can be
    inlined at its two call sites (the direct answer to :101);
  - the indexed PIP kernel's wrap-around branch (`i == n ? 1 : i + 1`)
    dies (prepared trees index exactly n−1 edges);
  - `_edge_extents`'s wrap logic survives only for ephemeral trees built
    over raw unclosed rings (get-or-create idiom), and says so.
- **Preserve point representations** (:290): if a curve's points are
  already `UnitSphericalPoint`s (or the storage is a reusable isbits
  vector), keep that storage instead of destructuring to lon/lat tuples.
  Concretely `_tuple_points` becomes manifold/storage-aware
  materialization.
- **Thread the manifold** (:224, :305): `prepare(m::Manifold, geom; …)`
  with `prepare(geom; …) = prepare(Planar(), geom; …)` per GO convention;
  `default_preparations(m, trait, geom)` and `build_edge_tree(m, backend,
  curve)` gain the manifold argument.  Planar behavior is unchanged; the
  spherical methods (arc extents via `arc_extent`, cap-driven queries)
  land from the `unit-spherical-indexing` branch later — this creates the
  seam they plug into.
- Adopt the suggested collapse
  `default_preparations(::GI.AbstractCurveTrait, geom) = (EdgeTree,)`
  (one method instead of ring + linestring).

## Theme E — one spec concept, no separate selector (:377, :224)

`EdgeTrees` (a selector) vs `EdgeTree` (a preparation) is two names for one
idea.  Unify (decision: **curried constructors**, not pair syntax):

- Delete `EdgeTrees`.  A `preps` tuple entry is a *spec*: a preparation
  type (`EdgeTree`), a curried constructor (`EdgeTree(HPR())` — the
  `EdgeTree` constructor applied to a backend returns a spec that builds
  `EdgeTree(geom; backend)` when applied to a geometry), or a closure.
  The `EdgeTree{HPR}` type-param spelling is *not* supported: `EdgeTree`'s
  type parameter is the stored tree type, and overloading it to mean
  "backend" in spec position would make one parameter mean two things.
- Each spec declares where it applies:
  `appliesto(::Type{EdgeTree}, trait) = trait isa GI.AbstractCurveTrait`
  (likewise for the curried form).  Specs without a declaration keep
  today's semantics (top node only), so `prepare(poly; preps = (MyPrep,))`
  still means "custom prep on the top node" — but
  `preps = (EdgeTree(HPR()),)` now flows to every curve, replacing
  `preps = EdgeTrees(HPR())`.
- This also answers "how do you indicate what kind of edge tree?": the
  spec carries the backend; the manifold comes from `prepare`'s manifold
  argument (Theme D), not from the spec.

## Theme F — closure and allocation hygiene (:401, :588)

- Replace `j -> push!(query_result, offset + j)` and the dual
  `(i, j) -> push!(candidate_pairs, (off_a + i, off_b + j))` with a small
  callable struct (`_OffsetPush`), per the review's suggestion.  (Done as
  part of Theme C's rework, this round.)
- The indexed PIP kernel (:588) — **deferred to the next iteration**, design
  recorded here: a `let` block would *not* remove the allocation — the
  `Ref` exists because the closure mutates shared state, and any
  boxed/mutable captured counter allocates.  The fix is to drop the
  closure entirely: hand-roll the recursive tree walk as a plain function
  returning `(on_hit, crossings)`, with early return on `on`.  This also
  deletes the `Action(:full_return, …)` unwinding.  Land with an
  `@allocated == 0` test for the indexed path.

## Theme G — prep lookup from the type tuple (:193, :197)

`_first_prep`'s value-level tuple recursion likely constant-folds already,
but it reads as a runtime scan and nothing proves it.  Replace with an
explicit compile-time query of the preps type-tuple: an `@generated`
`getprep(p::Prepared{…,P}, ::Type{Q})` that finds the first
`fieldtype(P, i) <: Q` and splices `p.preps[$i]` (or `nothing`) — zero
runtime search by construction, and `hasprep` follows.  Add `@inferred`
tests for hit, miss, and abstract-type queries.

## Theme H — stop reallocating already-native rings (union :142, :194)

**Deferred to the next iteration**, design recorded here.
`_linearring(tuples(ih))` copies every hole even when it is already native
tuple storage (which prepared inputs always are, post-Theme D).  Decision:
the deep clipping helpers should not each worry about storage at all —
normalization belongs at the **top-level ingestion layer**: the clipping
entry points (`union`/`intersection`/`difference`) convert each input once
into fast-access form, and everything below assumes it.  The helper is
named `fast_access_ring` (`Prepared` ring → its stored ring; a ring
already backed by `Vector{NTuple{2,T}}` → itself; anything else →
materialize).  The per-site `_linearring(tuples(ih))` calls then become
plain `_linearring(ih)` because ingestion already guaranteed the storage.
Verification step before landing: confirm clipping never mutates hole ring
point-vectors in place (it appends rings to `poly.geom` vectors, which is
fine — aliasing the *ring vector* is only safe if nobody `push!`es points
into output rings; audit `_add_union_holes*` and `_add_holes_to_polys!`),
otherwise copy exactly at the mutation site.

## Extending the thrust to the rest of the diff

- **intersection.jl**: restore the "iterate over each pair of maybe
  intersecting edges" framing comment and the LoopStateMachine teaching
  comment deleted by the tightening pass.
- **sutherland_hodgman.jl**: restore the three antipodal-candidate
  comments if our commits deleted them (check `git log -S`).
- **hilbert.jl, test/prepared.jl, test/utils/FlexibleRTrees.jl,
  benchmarks/prepared_natural_earth.jl**: audit for story-register
  comments and closure captures; fix in place.
- **prepared.jl module header**: the "Extensibility" / "Materialization"
  sections are user-facing docs and mostly behavioral — keep, but re-read
  against the doctrine and trim anything narrating design history; update
  for Themes D/E (manifold argument, spec protocol, ring closing).
- **Design docs**: update `2026-07-03-prepared-v2-recursive-materialization.md`
  with the revised accelerator vocabulary and invariants; record this plan's
  results section when done.

## Sequencing (each step lands green)

1. **Comments only** (Theme A + sweep): restores + register rewrites +
   verbatim suggestions.  No behavior change; clean first commit to
   re-review against.
2. **`prepare` normalization** (Theme D): manifold threading, ring
   closing, storage preservation; then delete the consumer-side checks it
   obsoletes (ggp `_point_ring_orientation`, `_edge_tree_and_coords`
   closedness, PIP wrap branch).
3. **Accelerator vocabulary** (Theme B): `TreeAccelerator` + per-side
   policies + `AutoAccelerator` decision table incl. both-prepared dual
   reuse; aliases for old names.
4. **Parts machinery** (Theme C + F): single accessor, inlined loops,
   domain names, `_OffsetPush`.
5. **Static prep lookup** (Theme G) with inference tests.
6. **Spec unification** (Theme E): delete `EdgeTrees`, curried
   `EdgeTree(backend)` specs, `appliesto` protocol, docs.
7. **Validation**: full test suite; re-run
   `benchmarks/prepared_natural_earth.jl` and the clipping benchmark
   against pre-change numbers — the ring-reuse and accelerator speedups
   must not regress.  Update design docs; draft a reply to each of the 27
   comments pointing at the commit that addresses it.

Next iteration (written down, not this round): the indexed-PIP walker
replacing the `Ref` closure (Theme F, :588) and the `fast_access_ring`
top-level ingestion for clipping (Theme H).

Decisions confirmed with review (2026-07-04): Themes A–D as written; old
accelerator names kept as aliases; Theme E via curried constructors
(`EdgeTree(HPR())`), no pair syntax; Theme G executed this round; Themes
F(:588)/H deferred as above.
