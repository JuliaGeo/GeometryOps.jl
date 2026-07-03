# Prepared v2: recursive materialization — decisions (2026-07-03)

Ratified by user, written down immediately (low context). This amends
`2026-07-03-prepared-minimal-design.md`; v1 is implemented on the
`prepared-minimal` branch (8 commits, HEAD a85a1e33a).

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
