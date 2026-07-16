# RelateNG Audit Fixes Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix the findings from the 2026-06-11 audit of `src/methods/geom_relations/relateng/`: one real API-misuse bug (`GI.extent` positional argument), one Base-generic hijack (`Base.merge!`), one kernel-contract fragility (`comp == 1`), a namespace-collision rename sweep, and minor style/documentation items.

**Architecture:** All changes are local to `src/methods/geom_relations/relateng/` and its test directory `test/methods/relateng/`. No public API changes except: `string(im::DE9IM)` keeps its value but `print(io, im)` changes from the `show` form to the bare matrix string. The renames only touch internal (unexported, underscore-convention) functions; the tests that call them via `GO.` are updated in the same commits.

**Tech Stack:** Julia, GeometryOps.jl monorepo (flat module), GeoInterface.jl (`GI`), Extents.jl. Tests run with `julia --project=test` from the `GeometryOps.jl/` repo root (per AGENTS.md). Each relateng test file is self-contained and can be `include`d directly.

**Important context for the executor:**

- Work from `/Users/anshul/temp/GO_jts/GeometryOps.jl` (the git repo root). The parent `GO_jts/` directory is NOT a git repo.
- If `julia --project=test` complains about missing packages, run `julia --project=test -e 'using Pkg; Pkg.instantiate()'` once first.
- Single test files run in ~seconds to ~2 minutes. The FULL relateng suite (`test/methods/relateng/runtests.jl`) takes ~25 minutes (JTS XML conformance + LibGEOS differential fuzz) — only run it once, at the end, in a background shell.
- Commit style (per AGENTS.md): imperative mood, capitalized, no trailing period, **no** `feat:`/`fix:` prefixes, backticks around code identifiers.
- The relateng port convention: files diff against their JTS Java counterparts; method order parallels the Java files. Do not reorder functions while editing.
- **Plan reconciled to commit `2a8ad4eab`** (2026-06-11): five performance commits landed after this plan was written (extent-cache wrapper tree `_rce`/`_relate_cache_extents`, NaturalIndex edge indexes, box-overlap node clustering, `IndexedPointInAreaLocator` in the new `indexed_point_in_area.jl` + its test file). All line references below were re-verified against that commit. The new `indexed_point_in_area.jl` is already trait-dispatch style and adds no rename targets (its `locate` method joins the existing keep-list `locate` generic).

---

### Task 1: Fix `GI.extent(c, true)` silently returning `nothing`

**The bug:** `GeoInterface.extent`'s `fallback` parameter is a *keyword* (`extent(obj; fallback=true)`). Passing it positionally as `GI.extent(c, true)` dispatches to the trait-form fallback `GI.extent(::Any, x) = Extents.extent(x)` — i.e. `c` is treated as a trait and `true` as the geometry — and `Extents.extent(true)` returns `nothing`. Net effect: in `_union_stored_extents`, any non-empty child without a stored extent (e.g. a Point member of a GeometryCollection) is silently dropped from the cached union extent.

**Files:**
- Modify: `src/methods/geom_relations/relateng/relate_geometry.jl:175`
- Test: `test/methods/relateng/relate_geometry.jl` (append a testset)

**Step 1: Write the failing test**

Open `test/methods/relateng/relate_geometry.jl`, look at the top of the file to confirm it imports `GeometryOps as GO` and `GeoInterface as GI` (it is self-contained per house convention). Append this testset at the end of the file:

```julia
@testset "GC extent cache includes point members" begin
    # A GC whose point member lies far outside its polygon member.
    # _union_stored_extents must include the point's coordinates in the
    # cached wrapper extent (audit 2026-06-11: GI.extent(c, true) passed
    # `fallback` positionally and silently returned `nothing`).
    gc = GI.GeometryCollection([
        GI.Point(10.0, 10.0),
        GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]]),
    ])
    cached = GO._relate_cache_extents(GO.Planar(), gc)
    ext = GI.extent(cached)
    @test ext.X[2] == 10.0
    @test ext.Y[2] == 10.0
end
```

**Step 2: Run the test to verify it fails**

Run from the repo root:
```bash
julia --project=test -e 'include("test/methods/relateng/relate_geometry.jl")'
```
Expected: the new testset FAILS with `ext.X[2] == 1.0` (the polygon-only extent). All pre-existing testsets in the file must PASS. If anything else fails, stop and investigate before proceeding.

**Step 3: Fix the call**

In `src/methods/geom_relations/relateng/relate_geometry.jl`, inside `_union_stored_extents` (defined at line 167; the call is at line 175), change:

```julia
        else
            GI.extent(c, true)
        end
```
to:
```julia
        else
            GI.extent(c; fallback = true)
        end
```

**Step 4: Check there are no other positional-Bool `GI.extent` calls**

```bash
grep -rn "GI.extent(" src/ | grep -v "fallback"
```
Expected: every hit either takes a single argument or is unrelated. If any other call passes a second positional argument, fix it the same way.

**Step 5: Run the test to verify it passes**

```bash
julia --project=test -e 'include("test/methods/relateng/relate_geometry.jl")'
```
Expected: ALL testsets PASS, including the new one.

**Step 6: Commit**

```bash
git add src/methods/geom_relations/relateng/relate_geometry.jl test/methods/relateng/relate_geometry.jl
git commit -m 'Fix GC extent caching dropping point members (positional `GI.extent` fallback arg)'
```

---

### Task 2: Rename the `Base.merge!` method on `RelateEdge` to `merge_edge!`

**The problem:** `Base.merge!(e::RelateEdge, is_a, dir_pt, dim, is_forward)` extends a Base generic whose contract is collection merging with unrelated semantics (merging edge labels). Not type piracy (we own `RelateEdge`), but it pollutes a Base generic. The Java method is named `merge`; the port-parallel convention only needs a recognizably parallel name, so `merge_edge!` is fine.

**Files:**
- Modify: `src/methods/geom_relations/relateng/relate_node.jl:209` (definition) and `:555` (callsite inside `add_edge!`)
- Test: `test/methods/relateng/node_topology.jl` (append a testset)

**Step 1: Write the failing test**

Append to `test/methods/relateng/node_topology.jl`:

```julia
@testset "no GO methods on Base.merge!" begin
    # RelateEdge label merging must not extend Base.merge! (collection
    # semantics). It lives on the internal function `merge_edge!` instead.
    @test !any(m -> m.module == GO, methods(Base.merge!))
end
```

(Confirm the file imports `GeometryOps as GO` at the top; it does, per house convention.)

**Step 2: Run the test to verify it fails**

```bash
julia --project=test -e 'include("test/methods/relateng/node_topology.jl")'
```
Expected: the new testset FAILS (one `Base.merge!` method is defined in `GeometryOps`). Everything else PASSES.

**Step 3: Rename definition and callsite**

In `src/methods/geom_relations/relateng/relate_node.jl`:

1. Around line 209, change the definition line
   `function Base.merge!(e::RelateEdge, is_a::Bool, dir_pt, dim::Integer, is_forward::Bool)`
   to
   `function merge_edge!(e::RelateEdge, is_a::Bool, dir_pt, dim::Integer, is_forward::Bool)`
2. Around line 555 (inside `add_edge!`), change
   `merge!(e, is_a, dir_pt, dim, is_forward)`
   to
   `merge_edge!(e, is_a, dir_pt, dim, is_forward)`
3. In the `#= Port of RelateEdge.merge(...) =#` comment block above the definition, no change needed — it already cites the Java name.

**Step 4: Verify no remaining callers**

```bash
grep -rn "merge!(" src/methods/geom_relations/relateng/ test/methods/relateng/
```
Expected: only `merge_edge!` hits (plus the new test's `methods(Base.merge!)` introspection line). No bare `merge!(e, ...)` calls remain.

**Step 5: Run the test to verify it passes**

```bash
julia --project=test -e 'include("test/methods/relateng/node_topology.jl")'
```
Expected: ALL testsets PASS.

**Step 6: Commit**

```bash
git add src/methods/geom_relations/relateng/relate_node.jl test/methods/relateng/node_topology.jl
git commit -m 'Rename `Base.merge!` method on `RelateEdge` to internal `merge_edge!`'
```

---

### Task 3: Harden `add_edge!` against the kernel sign contract

**The problem:** `rk_compare_edge_dir`'s documented contract (kernel.jl) is "negative / zero / positive", but `add_edge!` tests `comp == 1`. The planar kernel happens to return exactly ±1/0; a future kernel returning any other positive value would silently corrupt the edge-wheel ordering.

**Note on TDD:** No failing test is possible here — every current kernel returns exactly ±1/0, so the change has no observable behavior under the existing implementations. This is a contract-hardening change verified by the existing suite. Do NOT build a mock manifold to force it (YAGNI).

**Files:**
- Modify: `src/methods/geom_relations/relateng/relate_node.jl:558`

**Step 1: Make the change**

In `add_edge!` (around line 558), change:
```julia
        if comp == 1
```
to:
```julia
        if comp > 0
```

**Step 2: Check for other exact-sign comparisons against kernel results**

```bash
grep -rn "== 1\|== -1" src/methods/geom_relations/relateng/
```
Expected: any hits are NOT comparisons of `rk_compare_edge_dir`/`compare_to`-style results (e.g. counts are fine). If another kernel-comparison result is tested with `== 1`/`== -1`, change it to `> 0`/`< 0` the same way.

**Step 3: Run the node-topology tests**

```bash
julia --project=test -e 'include("test/methods/relateng/node_topology.jl")'
```
Expected: ALL PASS (no behavior change).

**Step 4: Commit**

```bash
git add src/methods/geom_relations/relateng/relate_node.jl
git commit -m 'Compare kernel edge-angle results by sign in `add_edge!`'
```

---

### Task 4: Replace `Base.string(::DE9IM)` with `Base.print`

**The problem:** Defining `Base.string` directly breaks the invariant `string(x) == sprint(print, x)`: `string(im)` gives `"212101212"` while `print(io, im)` falls through to `show` and gives `DE9IM("212101212")`. The Julian fix is to define `Base.print` and let `string` fall out of it.

**Files:**
- Modify: `src/methods/geom_relations/relateng/de9im.jl:102`
- Test: `test/methods/relateng/de9im.jl` (around line 41, where `string(im)` is already tested)

**Step 1: Write the failing test**

In `test/methods/relateng/de9im.jl`, directly after the existing line
`@test string(im) == "212101212"` (line ~41), add:

```julia
    @test sprint(print, im) == "212101212"   # print and string must agree
    @test sprint(show, im) == "DE9IM(\"212101212\")"
```

**Step 2: Run the test to verify it fails**

```bash
julia --project=test -e 'include("test/methods/relateng/de9im.jl")'
```
Expected: the `sprint(print, im)` assertion FAILS (currently yields `DE9IM("212101212")` via `show`). The `sprint(show, ...)` and `string(...)` assertions PASS.

**Step 3: Swap the definition**

In `src/methods/geom_relations/relateng/de9im.jl`, change line 102 from:
```julia
Base.string(im::DE9IM) = join(dim_char(d) for d in im.entries)
```
to:
```julia
# `string(im)` and `"$im"` yield the standard 9-character matrix form via
# this `print`; `show` keeps the constructor form for the REPL.
Base.print(io::IO, im::DE9IM) = join(io, (dim_char(d) for d in im.entries))
```
Leave the `Base.show` method on the next line unchanged. Note `string(im)` continues to work — Julia's generic `string(x)` is `sprint(print, x)`.

**Step 4: Check existing `string(im)` callers still resolve**

```bash
grep -rn "string(p.im)\|string(im)" src/methods/geom_relations/relateng/
```
Expected: hits in `topology_predicate.jl` (`Base.show(::IMPredicate)`) and `de9im.jl` — all fine, they go through generic `string`.

**Step 5: Run the test to verify it passes**

```bash
julia --project=test -e 'include("test/methods/relateng/de9im.jl")'
```
Expected: ALL PASS.

**Step 6: Commit**

```bash
git add src/methods/geom_relations/relateng/de9im.jl test/methods/relateng/de9im.jl
git commit -m 'Define `Base.print` for `DE9IM` instead of `Base.string`'
```

---

### Task 5: Throw `ArgumentError` instead of bare `error(...)` in `TopologyComputer`

**The problem:** Three unreachable-state guards use `error("Unknown target dimension: ...")` (→ `ErrorException`); the rest of the folder throws typed `ArgumentError`s.

**Files:**
- Modify: `src/methods/geom_relations/relateng/topology_computer.jl:293,326,400`
- Test: `test/methods/relateng/topology_computer.jl` (append a testset)

**Step 1: Write the failing test**

Open `test/methods/relateng/topology_computer.jl` and find how existing testsets construct a `TopologyComputer` (they build two `GO.RelateGeometry`s and pass a predicate — mirror that construction exactly, including the `exact = ...` keyword they use). Append:

```julia
@testset "unknown target dimension throws ArgumentError" begin
    a = GI.Point(0.0, 0.0)
    b = GI.Point(1.0, 1.0)
    rga = GO.RelateGeometry(GO.Planar(), a; exact = GO.True())
    rgb = GO.RelateGeometry(GO.Planar(), b; exact = GO.True())
    tc = GO.TopologyComputer(GO.pred_intersects(), rga, rgb)
    @test_throws ArgumentError GO.add_point_on_geometry!(
        tc, true, GO.LOC_INTERIOR, Int8(9), (0.0, 0.0))
end
```

(If the file's existing tests spell the exactness flag differently — e.g. `GO.GeometryOpsCore.True()` — match them.)

**Step 2: Run the test to verify it fails**

```bash
julia --project=test -e 'include("test/methods/relateng/topology_computer.jl")'
```
Expected: the new testset FAILS — an `ErrorException` is thrown where `ArgumentError` is expected.

**Step 3: Replace the three `error` calls**

In `src/methods/geom_relations/relateng/topology_computer.jl`, at the ends of `add_point_on_geometry!` (~line 293), `add_line_end_on_geometry!` (~line 326), and `add_area_vertex!` (~line 400), change each:
```julia
    error("Unknown target dimension: $dim_target")
```
to:
```julia
    throw(ArgumentError("unknown target dimension: $dim_target"))
```

**Step 4: Run the test to verify it passes**

```bash
julia --project=test -e 'include("test/methods/relateng/topology_computer.jl")'
```
Expected: ALL PASS.

**Step 5: Commit**

```bash
git add src/methods/geom_relations/relateng/topology_computer.jl test/methods/relateng/topology_computer.jl
git commit -m 'Throw `ArgumentError` for unknown target dimensions in `TopologyComputer`'
```

---

### Task 6: Rename collision-prone generic internal names

**The problem:** The flat GO module gains unexported internal functions with extremely generic names (`id`, `dimension`, `location`, `matches`, `is_empty`). No collision exists today (verified by grep), but any future GO file defining the same name silently merges methods. We rename the five riskiest; the rest (`locate`, `finish!`, `compare_to`, `get_*`, `is_known`, …) are deliberately KEPT — `locate` is a legitimate multi-locator generic, and the others are subsystem-scoped port-parity names whose rename would cost more diffability than it buys.

**Rename map (definition → all callsites, src AND test):**

| Old | New | Definition | Known callsites |
|---|---|---|---|
| `id(ns::NodeSection)` | `section_id` | `node_sections.jl:94` | `node_sections.jl:126` (`is_same_polygon`, two calls), `polygon_node_converter.jl:107`, `test/.../node_topology.jl:32` |
| `dimension(ns::NodeSection)` | `section_dim` | `node_sections.jl:91` | `node_sections.jl:80` (`is_area_area`, two calls), `relate_node.jl:464`, `test/.../node_topology.jl:31` |
| `dimension(e::RelateEdge, is_a)` | `edge_dim` | `relate_node.jl:351` | `relate_node.jl:243` |
| `location(e::RelateEdge, is_a, pos)` | `edge_location` | `relate_node.jl:337` | `relate_node.jl:254,360,364,407–410,609,615,640–641`; `topology_computer.jl:658–664`; `test/.../node_topology.jl:188–189,340–341,381,383,446–447` |
| `matches(im::DE9IM, pattern)` | `im_matches` | `de9im.jl:112` | `relate_predicates.jl:283` (`value_im(::IMPatternMatcher)`); `test/.../de9im.jl:48–55` |
| `is_empty(rg::RelateGeometry)` | `is_geom_empty` | `relate_geometry.jl:397` | `topology_computer.jl:271,314` |

When renaming, keep the adjacent "Port of `<JavaClass.method>`" comments exactly as they are — they record the Java name, which is the point.

**Step 1: For each row, find every occurrence**

For each old name, run a word-boundary grep BEFORE editing, e.g.:
```bash
grep -rnw "dimension" src/methods/geom_relations/relateng/ test/methods/relateng/
```
Cross-check the hits against the table. If you find a callsite NOT in the table, update it too — the table is the audit's snapshot, the grep is the truth.

**Step 2: Apply the renames**

Edit each definition and callsite. Careful with `dimension`: it has TWO methods becoming TWO differently-named functions (`section_dim` for `NodeSection`, `edge_dim` for `RelateEdge`) — match by argument type at each callsite. Careful with `is_empty` vs. the *field* `rg.is_geom_empty` (the accessor becomes `is_geom_empty(rg) = rg.is_geom_empty`, which is fine — function and field namespaces don't clash). Do NOT touch `GI.isempty` or Julia's `isempty`.

**Step 3: Verify zero stale references**

```bash
grep -rnw -e "id" -e "dimension" -e "matches" src/methods/geom_relations/relateng/ | grep -v "_id\|section_id\|section_dim\|edge_dim\|im_matches\|matches_entry\|pred_matches\|ring_id\|element_id\|node_id\|#"
grep -rnw "is_empty" src/methods/geom_relations/relateng/
grep -rnw -e "GO.id" -e "GO.dimension" -e "GO.location" -e "GO.matches" -e "GO.is_empty" test/methods/relateng/
```
Expected: no hits that are calls/definitions of the old names (comments citing Java method names like `NodeSection.id` are fine and should remain).

**Step 4: Run the affected unit test files**

```bash
julia --project=test -e 'include("test/methods/relateng/de9im.jl")'
julia --project=test -e 'include("test/methods/relateng/node_topology.jl")'
julia --project=test -e 'include("test/methods/relateng/topology_computer.jl")'
julia --project=test -e 'include("test/methods/relateng/predicates.jl")'
```
Expected: ALL PASS. (`MethodError`/`UndefVarError` here means a missed callsite — go back to Step 3.)

**Step 5: Commit**

```bash
git add src/methods/geom_relations/relateng/ test/methods/relateng/
git commit -m 'Rename collision-prone internal relateng helpers (`id`, `dimension`, `location`, `matches`, `is_empty`)'
```

---

### Task 7: Convert trait `isa` ladders to trait-dispatched methods

**The problem:** Some relateng walks classify geometry by `trait = GI.trait(geom)` followed by an `if trait isa X ... elseif trait isa Y ...` ladder, while neighboring functions in the same folder (`_rce`, `_iig_add_geom!`, `_add_rings!`, `_locate_point_in_polygonal`) use the house style: trait-dispatched methods. Unify on dispatch.

**Scope — convert ONLY pure classification ladders** (functions whose body is an if/elseif over the trait covering several geometry classes). Explicitly KEEP as-is:

- Single-trait guards inside larger logic (`is_polygonal`, `locate_with_dim`'s `isa Union{...}` check, `locate_node`/`is_node_in_area` in topology_computer.jl, the inline checks in `_extract_segment_strings!`).
- The sequential walks `_compute_line_ends_walk!` / `_compute_area_vertex_walk!` (relate_ng.jl) and `_extract_segment_strings_from_atomic!` — they mix trait checks with threaded state/side effects and are guard-style, not ladders.

**Functions to convert:**

| Function | File |
|---|---|
| `_extract_elements!` | `point_locator.jl:305` (trait-form method; the 4-arg entry at `:302` already exists) |
| `_relate_is_empty` | `relate_geometry.jl:85` |
| `_relate_extent` | `relate_geometry.jl:99` |
| `_geom_dimension` | `relate_geometry.jl:187` |
| `_analyze_dimensions` + element classification in `_analyze_collection_dimensions` | `relate_geometry.jl:211,227` |
| `_is_zero_length` | `relate_geometry.jl:259` |
| `_add_component_coordinates!` | `relate_geometry.jl:419` |
| `_extract_point_elements!` | `relate_geometry.jl:457` |
| `_segment_string_eltype` | `relate_geometry.jl:497` |

**⚠️ Correctness note for the executor:** in the ladders, *branch order* encodes precedence. In GeoInterface, the `Multi*` traits are subtypes of `GI.AbstractGeometryCollectionTrait`, and the ladders check the point/curve/polygon unions BEFORE the GC branch. The converted methods rely on dispatch specificity instead: a `Union` containing `GI.AbstractMultiPointTrait` (etc.) is strictly more specific than `GI.AbstractGeometryCollectionTrait`, so dispatch picks the same branch. The characterization testset in Step 1 pins exactly these cases — do not skip it.

**Files:**
- Modify: `src/methods/geom_relations/relateng/point_locator.jl`, `src/methods/geom_relations/relateng/relate_geometry.jl`
- Test: `test/methods/relateng/relate_geometry.jl` (append characterization testset)

**Step 1: Write the characterization test (pins behavior BEFORE the refactor)**

Append to `test/methods/relateng/relate_geometry.jl`. This is a refactor, so the test must PASS both before and after — it exists to catch dispatch-specificity mistakes (especially `Multi*` traits vs. the GC trait):

```julia
@testset "trait classification characterization" begin
    mp    = GI.MultiPoint([(0.0, 0.0), (1.0, 1.0)])
    ml    = GI.MultiLineString([[(0.0, 0.0), (1.0, 0.0)], [(0.0, 1.0), (1.0, 1.0)]])
    mpoly = GI.MultiPolygon([[[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]]])
    zline = GI.LineString([(0.0, 0.0), (0.0, 0.0)])   # zero-length
    gc    = GI.GeometryCollection([GI.Point(5.0, 5.0), mpoly])
    empty_gc = GI.GeometryCollection([GI.LineString(Tuple{Float64, Float64}[])])

    #-- _geom_dimension: Multi* must classify by content, not as GC
    @test GO._geom_dimension(mp) == GO.DIM_P
    @test GO._geom_dimension(ml) == GO.DIM_L
    @test GO._geom_dimension(mpoly) == GO.DIM_A
    @test GO._geom_dimension(gc) == GO.DIM_A

    #-- _relate_is_empty: recursive emptiness through collections
    @test GO._relate_is_empty(empty_gc)
    @test !GO._relate_is_empty(gc)

    #-- _is_zero_length: MultiLineString takes the collection branch
    @test GO._is_zero_length(zline)
    @test !GO._is_zero_length(ml)

    #-- _relate_extent: union over GC elements including the bare point
    @test GO._relate_extent(GO.Planar(), gc).X == (0.0, 5.0)

    #-- _analyze_dimensions via the constructor: mixed GC
    rg = GO.RelateGeometry(GO.Planar(), gc; exact = GO.True())
    @test GO.get_dimension(rg) == GO.DIM_A
    @test rg.has_points && rg.has_areas && !rg.has_lines

    #-- _add_component_coordinates!: point coords + first line coords
    pts = Set{Tuple{Float64, Float64}}()
    GO._add_component_coordinates!(pts, GI.GeometryCollection([mp, ml]))
    @test (0.0, 0.0) in pts

    #-- _extract_point_elements!: only the Point element of the GC
    lst = Any[]
    GO._extract_point_elements!(lst, gc)
    @test length(lst) == 1

    #-- _extract_elements!: classification into points/lines/polygons
    points = Set{Tuple{Float64, Float64}}(); lines = Any[]; polygons = Any[]
    GO._extract_elements!(points, lines, polygons,
        GI.GeometryCollection([GI.Point(0.0, 0.0), ml, mpoly]))
    @test length(points) == 1 && length(lines) == 2 && length(polygons) == 1
end
```

(If the existing testsets construct `RelateGeometry` with a different exactness spelling than `GO.True()`, match them.)

**Step 2: Run the test — it must PASS against the CURRENT code**

```bash
julia --project=test -e 'include("test/methods/relateng/relate_geometry.jl")'
```
Expected: ALL PASS. If the new testset fails here, the test is wrong (it mischaracterizes current behavior) — fix the test, not the source.

**Step 3: Convert `_extract_elements!` (point_locator.jl)**

Replace the trait-form method `_extract_elements!(points, lines, polygons, trait::GI.AbstractTrait, geom)` (lines ~305–323) with dispatch methods. Keep the 4-argument entry point above it (`point_locator.jl:302–303`) and its "Port of RelatePointLocator.extractElements" comment unchanged. The original ladder checks the polygonal union before the GC branch — dispatch specificity reproduces this (see the correctness note):

```julia
function _extract_elements!(points, lines, polygons, ::GI.PointTrait, geom)
    GI.isempty(geom) && return nothing
    #-- addPoint: normalized coordinate tuples, as in LinearBoundary
    push!(points, _node_point(geom))
    return nothing
end
function _extract_elements!(points, lines, polygons, ::GI.AbstractCurveTrait, geom)
    GI.isempty(geom) && return nothing
    #-- addLine (Java LinearRing extends LineString, hence AbstractCurve)
    push!(lines, geom)
    return nothing
end
function _extract_elements!(points, lines, polygons,
        ::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, geom)
    GI.isempty(geom) && return nothing
    #-- addPolygonal: whole polygonal geometry kept as one element
    push!(polygons, geom)
    return nothing
end
function _extract_elements!(points, lines, polygons,
        ::GI.AbstractGeometryCollectionTrait, geom)
    GI.isempty(geom) && return nothing
    #-- covers GeometryCollection, MultiPoint, MultiLineString
    for g in GI.getgeom(geom)
        _extract_elements!(points, lines, polygons, g)
    end
    return nothing
end
_extract_elements!(points, lines, polygons, ::GI.AbstractTrait, geom) = nothing
```

**Step 4: Convert the relate_geometry.jl ladders**

Replace each function body with dispatch methods, preserving the existing comments by moving them onto the matching method. Trait argument goes before the geometry, matching `_rce`/`_iig_add_geom!`.

`_relate_is_empty`:
```julia
_relate_is_empty(geom) = _relate_is_empty(GI.trait(geom), geom)
function _relate_is_empty(::GI.AbstractGeometryCollectionTrait, geom)
    for g in GI.getgeom(geom)
        _relate_is_empty(g) || return false
    end
    return true
end
_relate_is_empty(::GI.AbstractTrait, geom) = GI.isempty(geom)
```

`_relate_extent`:
```julia
_relate_extent(m::Manifold, geom) = _relate_extent(m, GI.trait(geom), geom)
function _relate_extent(m::Manifold, ::GI.AbstractGeometryCollectionTrait, geom)
    ext = nothing
    for g in GI.getgeom(geom)
        e = _relate_extent(m, g)
        e === nothing && continue
        ext = ext === nothing ? e : Extents.union(ext, e)
    end
    return ext
end
function _relate_extent(m::Manifold, ::GI.AbstractTrait, geom)
    GI.isempty(geom) && return nothing
    return rk_interaction_bounds(m, geom)
end
```

`_geom_dimension`:
```julia
_geom_dimension(geom) = _geom_dimension(GI.trait(geom), geom)
_geom_dimension(::Union{GI.AbstractPointTrait, GI.AbstractMultiPointTrait}, geom) = DIM_P
_geom_dimension(::Union{GI.AbstractCurveTrait, GI.AbstractMultiCurveTrait}, geom) = DIM_L
_geom_dimension(::Union{GI.AbstractPolygonTrait, GI.AbstractMultiPolygonTrait}, geom) = DIM_A
function _geom_dimension(::GI.AbstractGeometryCollectionTrait, geom)
    dim = DIM_FALSE
    for g in GI.getgeom(geom)
        d = _geom_dimension(g)
        d > dim && (dim = d)
    end
    return dim
end
_geom_dimension(::GI.AbstractTrait, geom) = DIM_FALSE
```

`_analyze_dimensions` (keep the 3-arg entry; the GC/other fallback goes to the collection walk as before):
```julia
function _analyze_dimensions(geom, dim0::Int8, is_geom_empty::Bool)
    is_geom_empty && return (dim0, false, false, false)
    return _analyze_dimensions(GI.trait(geom), geom, dim0)
end
_analyze_dimensions(::Union{GI.AbstractPointTrait, GI.AbstractMultiPointTrait}, geom, dim0) =
    (DIM_P, true, false, false)
_analyze_dimensions(::Union{GI.AbstractCurveTrait, GI.AbstractMultiCurveTrait}, geom, dim0) =
    (DIM_L, false, true, false)
_analyze_dimensions(::Union{GI.AbstractPolygonTrait, GI.AbstractMultiPolygonTrait}, geom, dim0) =
    (DIM_A, false, false, true)
#-- analyze a (possibly mixed type) collection
_analyze_dimensions(::GI.AbstractTrait, geom, dim0) =
    _analyze_collection_dimensions(geom, dim0, false, false, false)
```

`_analyze_collection_dimensions` — extract the per-element ladder into a dispatch helper:
```julia
# The recursive element walk of analyzeDimensions (Java uses a
# GeometryCollectionIterator; only atomic elements match the checks).
function _analyze_collection_dimensions(geom, dim, has_points, has_lines, has_areas)
    for g in GI.getgeom(geom)
        dim, has_points, has_lines, has_areas = _analyze_element_dimensions(
            GI.trait(g), g, dim, has_points, has_lines, has_areas)
    end
    return (dim, has_points, has_lines, has_areas)
end
_analyze_element_dimensions(::GI.AbstractGeometryCollectionTrait, g, dim, hp, hl, ha) =
    _analyze_collection_dimensions(g, dim, hp, hl, ha)
function _analyze_element_dimensions(::GI.AbstractPointTrait, g, dim, hp, hl, ha)
    GI.isempty(g) && return (dim, hp, hl, ha)
    return (max(dim, DIM_P), true, hl, ha)
end
function _analyze_element_dimensions(::GI.AbstractCurveTrait, g, dim, hp, hl, ha)
    GI.isempty(g) && return (dim, hp, hl, ha)
    return (max(dim, DIM_L), hp, true, ha)
end
function _analyze_element_dimensions(::GI.AbstractPolygonTrait, g, dim, hp, hl, ha)
    GI.isempty(g) && return (dim, hp, hl, ha)
    return (max(dim, DIM_A), hp, hl, true)
end
_analyze_element_dimensions(::GI.AbstractTrait, g, dim, hp, hl, ha) = (dim, hp, hl, ha)
```

`_is_zero_length`:
```julia
_is_zero_length(geom) = _is_zero_length(GI.trait(geom), geom)
function _is_zero_length(::GI.AbstractGeometryCollectionTrait, geom)
    for g in GI.getgeom(geom)
        _is_zero_length(g) || return false
    end
    return true
end
_is_zero_length(::GI.AbstractCurveTrait, geom) = _is_zero_length_linestring(geom)
_is_zero_length(::GI.AbstractTrait, geom) = true
```

`_add_component_coordinates!` (the point/curve ternary splits into two methods):
```julia
_add_component_coordinates!(set, geom) =
    _add_component_coordinates!(set, GI.trait(geom), geom)
function _add_component_coordinates!(set, ::GI.AbstractGeometryCollectionTrait, geom)
    for g in GI.getgeom(geom)
        _add_component_coordinates!(set, g)
    end
    return nothing
end
function _add_component_coordinates!(set, ::GI.AbstractPointTrait, geom)
    GI.isempty(geom) && return nothing
    push!(set, _node_point(geom))
    return nothing
end
function _add_component_coordinates!(set, ::GI.AbstractCurveTrait, geom)
    GI.isempty(geom) && return nothing
    push!(set, _node_point(GI.getpoint(geom, 1)))
    return nothing
end
_add_component_coordinates!(set, ::GI.AbstractTrait, geom) = nothing
```

`_extract_point_elements!`:
```julia
_extract_point_elements!(list, geom) =
    _extract_point_elements!(list, GI.trait(geom), geom)
function _extract_point_elements!(list, ::GI.AbstractPointTrait, geom)
    push!(list, geom)
    return nothing
end
function _extract_point_elements!(list, ::GI.AbstractGeometryCollectionTrait, geom)
    for g in GI.getgeom(geom)
        _extract_point_elements!(list, g)
    end
    return nothing
end
_extract_point_elements!(list, ::GI.AbstractTrait, geom) = nothing
```

`_segment_string_eltype` (the separate Polygon/MultiPolygon branches had identical bodies — merge into one `Union` method; keep the explanatory comment block above the entry point):
```julia
_segment_string_eltype(rg::RelateGeometry, geom) =
    _segment_string_eltype(rg, GI.trait(geom), geom)
_segment_string_eltype(rg::RG, ::GI.AbstractCurveTrait, geom) where {RG <: RelateGeometry} =
    RelateSegmentString{Tuple{Float64, Float64}, Nothing, RG}
#-- rings of MultiPolygon elements carry the MultiPolygon as parent
_segment_string_eltype(rg::RG, ::Union{GI.AbstractPolygonTrait, GI.AbstractMultiPolygonTrait},
        geom) where {RG <: RelateGeometry} =
    RelateSegmentString{Tuple{Float64, Float64}, typeof(geom), RG}
function _segment_string_eltype(rg::RelateGeometry, ::GI.AbstractGeometryCollectionTrait, geom)
    T = Union{}
    for g in GI.getgeom(geom)
        T = Union{T, _segment_string_eltype(rg, g)}
    end
    return T
end
#-- point elements produce no segment strings
_segment_string_eltype(::RelateGeometry, ::GI.AbstractTrait, geom) = Union{}
```

**Step 5: Check for ladders you missed and stale callers**

```bash
grep -rn "trait isa" src/methods/geom_relations/relateng/
```
Expected: remaining hits are ONLY the keep-list items (single guards in `is_polygonal`, `locate_with_dim`, `locate_node`/`is_node_in_area`, `_extract_segment_strings!`/`_extract_segment_strings_from_atomic!`, and the relate_ng.jl walks). All converted functions must have no `trait isa` left.

**Step 6: Run the characterization test and neighbors — must still pass**

```bash
julia --project=test -e 'include("test/methods/relateng/relate_geometry.jl")'
julia --project=test -e 'include("test/methods/relateng/point_locator.jl")'
julia --project=test -e 'include("test/methods/relateng/relate_ng.jl")'
```
Expected: ALL PASS. A `MethodError: ... is ambiguous` here means a `Union` method conflicts with the GC-trait method for some `Multi*` trait — resolve by adding an explicit method for that concrete trait that forwards to the intended branch.

**Step 7: Commit**

```bash
git add src/methods/geom_relations/relateng/ test/methods/relateng/relate_geometry.jl
git commit -m 'Convert trait `isa` ladders in relateng to trait-dispatched methods'
```

---

### Task 8: Documentation and comment fixes

No code behavior changes; no tests. Four edits:

**Files:**
- Modify: `src/methods/geom_relations/relateng/de9im.jl` (DE9IM docstring, ~line 79)
- Modify: `src/methods/geom_relations/relateng/relate_ng.jl` (RelateNG docstring ~line 103, `prepare` docstring ~line 650, definition at `:657`)
- Modify: `src/methods/geom_relations/relateng/topology_predicate.jl` (`is_known_entry`, ~line 214)

**Step 1: DE9IM docstring — explain the index convention**

In the `DE9IM` docstring, replace the sentence `Index with `im[locA, locB]`.` with:

```
Index with `im[locA, locB]`, where the indices are the JTS location *codes*
(`0` = Interior, `1` = Boundary, `2` = Exterior — the internal `LOC_*`
constants), **not** 1-based array positions: `im[0, 0]` is the
Interior/Interior entry.
```

**Step 2: RelateNG docstring — document Float64 evaluation**

In the `RelateNG` docstring (after the numbered capability list, before the "Keyword arguments" paragraph), add:

```
All coordinates are evaluated as `Float64`: input coordinates are converted
on extraction, and the exact-predicate machinery (adaptive orientation
predicates, rational-arithmetic node coincidence) assumes `Float64` inputs.
Non-`Float64` geometries are accepted but evaluated at `Float64` precision.
```

**Step 3: `prepare` docstring — state the genericity intent**

At the top of the `prepare(alg::RelateNG, a)` docstring, after the signature line, add:

```
`prepare` is the generic entry point for prepared-geometry optimizations in
GeometryOps; `RelateNG` is currently the only algorithm implementing it.
```

**Step 4: `is_known_entry` — mark the parity dead code**

Above `is_known_entry` in `topology_predicate.jl` (~line 214), add:

```julia
# NOTE: unused; kept for JTS IMPredicate API parity. As ported it can never
# return `false`: matrix entries are initialized to DIM_FALSE and only ever
# increase, so they never hold DIM_UNKNOWN (-3). Do not use as a real check.
```

**Step 5: Sanity-check the doc edits compile**

```bash
julia --project=test -e 'include("test/methods/relateng/de9im.jl")'
```
Expected: PASS (docstrings/comments only; this just catches syntax slips).

**Step 6: Commit**

```bash
git add src/methods/geom_relations/relateng/
git commit -m 'Document DE9IM index codes, Float64 evaluation, `prepare` genericity, and `is_known_entry` parity status'
```

---

### Task 9: Full-suite verification

**Step 1: Run the full relateng suite in the background**

This takes ~25 minutes (XML conformance + LibGEOS differential fuzz). Run it in a background shell and poll:

```bash
cd /Users/anshul/temp/GO_jts/GeometryOps.jl && julia --project=test -e 'include("test/methods/relateng/runtests.jl")' > /tmp/relateng_suite.log 2>&1
```

**Step 2: Check the result**

```bash
tail -20 /tmp/relateng_suite.log
```
Expected: the standard `Test Summary` block with zero `Fail`/`Error`. If anything fails, identify which task introduced it (`git log --oneline` + the failing testset name), fix in place, and re-run only the affected unit file before re-running the full suite.

**Step 3: Review the final diff**

```bash
git log --oneline main..HEAD 2>/dev/null || git log --oneline -9
git diff HEAD~8 --stat
```
Expected: 8 commits, all changes confined to `src/methods/geom_relations/relateng/`, `test/methods/relateng/`, and this plan file's directory.

---

## Explicitly out of scope (decided during the audit)

- **Submodule restructuring** (`module RelateNG ... end`): rejected in favor of targeted renames — the package convention is a flat module (clipping is flat too), and wrapping 14 files would hurt the file-by-file Java diffability that the port is built around.
- **Renaming `get_`/`set_` accessors, `compare_to`, `finish!`, `locate`, `is_known`, `get_geometry`, …**: deliberate port-parity names; `locate` is a proper multi-type generic. Leave them.
- **Converting the single-trait guards and sequential walks to dispatch** (`is_polygonal`, `locate_with_dim`, `locate_node`/`is_node_in_area`, `_compute_line_ends_walk!`, `_compute_area_vertex_walk!`, `_extract_segment_strings_from_atomic!`): these mix trait checks with threaded state or are single predicates inside larger logic — Task 7 deliberately converts only the pure classification ladders.
- **`Vector{Any}` element collections in `RelatePointLocator` / abstract `Vector{NodeSection}`**: documented, deliberate, not on hot paths. No change.
- **Exporting `LOC_*` constants or symbol-based `DE9IM` indexing**: YAGNI for now; the docstring fix (Task 8) covers usability.
- **`IM_PATTERN_*` constants**: NOT dead code — they are used in `test/methods/relateng/relate_ng.jl`. Keep as-is.
