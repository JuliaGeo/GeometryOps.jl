# RelateNG Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement the RelateNG DE-9IM predicate engine in GeometryOps.jl per the validated design in `docs/plans/2026-06-10-relateng-design.md` — exact predicates, no constructed intersection coordinates, manifold-parameterized kernel, Planar first.

**Architecture:** Two layers. The *topology layer* (predicate framework, IM bookkeeping, point location, node-topology analysis, evaluation phases) is a faithful, file-by-file port of JTS `org.locationtech.jts.operation.relateng` so it stays diffable against the Java. The *geometry layer* is a small `rk_*` kernel API (orientation, symbolic segment-intersection classification, point-in-ring, edge ordering around symbolic nodes, bounds) implemented for `Planar` via AdaptivePredicates; it is the only place coordinates are touched, and it never constructs intersection points.

**Tech Stack:** Julia; GeometryOpsCore (`Algorithm{M}`, `Manifold`, `True`/`False`/`booltype`); `GO.Predicates` (AdaptivePredicates exact orient); SpatialTreeInterface (`dual_depth_first_search`, STRtree/FlatNoTree); Extents.jl; XML.jl + WellKnownGeometry (test harness); LibGEOS (differential oracle, GEOS ≥ 3.13 runs RelateNG natively).

---

## Conventions for every task (read once, apply always)

- **JTS source root:** `/Users/anshul/temp/GO_jts/jts/modules/core/src/main/java/org/locationtech/jts/operation/relateng/` (abbreviated `JTS:` below). JTS tests: `/Users/anshul/temp/GO_jts/jts/modules/core/src/test/java/org/locationtech/jts/operation/relateng/`. JTS XML: `/Users/anshul/temp/GO_jts/jts/modules/tests/src/test/resources/testxml/`.
- **Port protocol for "port from Java" steps:** Read the named Java file *in full* first. Port the JUnit tests for it first (RED), then the implementation (GREEN). Keep function-level structure parallel to the Java (same method names in snake_case, same order in file) so the Julia file diffs against its Java counterpart. Where this plan's code conflicts with the Java semantics, the Java wins — note the discrepancy in the commit message.
- **Run command:** `julia --project=docs <file>` or `julia --project=docs -e '...'` (per AGENTS.md the docs environment has GeometryOps devved with all dependencies). Verify once in Task 0.
- **Commit style (AGENTS.md, overrides generic templates):** imperative, capitalized, no `feat:`-style prefixes, no trailing period, backticks for code identifiers. Trailer: `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- **TDD:** every task is test-first. Write the failing test, run it, watch it fail for the right reason, implement, run again, commit. One commit per task minimum; commit more often if a task has natural checkpoints.
- **Naming:** all kernel functions are prefixed `rk_` (RelateKernel). Topology-layer types keep their JTS names (`NodeSection`, `TopologyComputer`, …). Internal constants are prefixed (`LOC_`, `DIM_`, `DL_`, `TRI_`). Nothing internal is exported.
- **Points are coordinate tuples** (`Tuple{Float64,Float64}` typically) obtained via `GO._tuple_point`. The `exact` flag is a keyword taking `True()`/`False()` (GeometryOpsCore BoolsAsTypes), threaded exactly like `Predicates.orient(a, b, c; exact)` in `src/methods/clipping/predicates.jl:13`.
- **File include discipline:** each new `src/.../relateng/*.jl` file gets its `include(...)` line added to `src/GeometryOps.jl` (after line 89, `include("methods/geom_relations/common.jl")`) in the same task that creates the file. Same for test files in `test/methods/relateng/runtests.jl`.

---

## Task 0: Workspace setup and test-command verification

**Files:**
- Create: `test/methods/relateng/runtests.jl` (empty shell)
- Modify: `test/runtests.jl`

**Step 1: Create the working branch/worktree**

The design is committed on branch `relateng-design`. Create the implementation branch from it:

```bash
cd /Users/anshul/temp/GO_jts/GeometryOps.jl
git checkout relateng-design
git checkout -b relateng
```

(If executing via the worktrees skill, create the worktree off `relateng-design` instead.)

**Step 2: Verify the test run command works**

```bash
julia --project=docs -e 'import GeometryOps, Test, Extents; println("ok")'
```

Expected: `ok`. If GeometryOps is not devved in the docs env, fix with `julia --project=docs -e 'using Pkg; Pkg.develop(path=".")'` and record the working command — it is used by every task below.

**Step 3: Create the relateng test shell and register it**

Create `test/methods/relateng/runtests.jl`:

```julia
using SafeTestsets

@safetestset "DE9IM" begin include("de9im.jl") end
# Further files appended here as tasks land:
# @safetestset "Predicates" begin include("predicates.jl") end
# @safetestset "Kernel" begin include("kernel.jl") end
# ...
```

(Comment lines are uncommented as each test file is created.)

In `test/runtests.jl`, after the line `@safetestset "DE-9IM Geom Relations" begin include("methods/geom_relations.jl") end`, add:

```julia
@safetestset "RelateNG" begin include("methods/relateng/runtests.jl") end
```

**Step 4: Commit**

```bash
git add test/methods/relateng/runtests.jl test/runtests.jl
git commit -m "Add RelateNG test scaffolding"
```

---

# Stage 1 — Foundations (combinatorial, no geometry)

## Task 1: Location/dimension codes, `DimensionLocation`, `DE9IM` type

**Files:**
- Create: `src/methods/geom_relations/relateng/de9im.jl`
- Modify: `src/GeometryOps.jl` (add include)
- Test: `test/methods/relateng/de9im.jl`

**Java references:** `org/locationtech/jts/geom/Location.java`, `Dimension.java`, `IntersectionMatrix.java` (matches semantics), `JTS:DimensionLocation.java` (89 lines, constants verbatim).

**Step 1: Write the failing tests** — `test/methods/relateng/de9im.jl`:

```julia
using Test
import GeometryOps as GO
import GeometryOps: DE9IM

@testset "codes" begin
    @test GO.LOC_INTERIOR == 0 && GO.LOC_BOUNDARY == 1 && GO.LOC_EXTERIOR == 2
    @test GO.DIM_FALSE == -1 && GO.DIM_P == 0 && GO.DIM_L == 1 && GO.DIM_A == 2
    @test GO.dim_char(GO.DIM_FALSE) == 'F'
    @test GO.dim_code('T') == GO.DIM_TRUE && GO.dim_code('*') == GO.DIM_DONTCARE
    @test_throws ArgumentError GO.dim_code('X')
end

@testset "DimensionLocation" begin
    # Constants verbatim from JTS DimensionLocation.java
    @test GO.DL_POINT_INTERIOR == 103
    @test GO.DL_LINE_INTERIOR == 110 && GO.DL_LINE_BOUNDARY == 111
    @test GO.DL_AREA_INTERIOR == 120 && GO.DL_AREA_BOUNDARY == 121
    @test GO.dimloc_location(GO.DL_AREA_BOUNDARY) == GO.LOC_BOUNDARY
    @test GO.dimloc_dimension(GO.DL_LINE_INTERIOR) == GO.DIM_L
    @test GO.dimloc_location(GO.DL_EXTERIOR) == GO.LOC_EXTERIOR
    @test GO.dimloc_area(GO.LOC_INTERIOR) == GO.DL_AREA_INTERIOR
    @test GO.dimloc_line(GO.LOC_BOUNDARY) == GO.DL_LINE_BOUNDARY
    @test GO.dimloc_point(GO.LOC_EXTERIOR) == GO.DL_EXTERIOR
end

@testset "DE9IM" begin
    im = DE9IM("212101212")
    @test string(im) == "212101212"
    @test im[GO.LOC_INTERIOR, GO.LOC_INTERIOR] == GO.DIM_A
    @test im[GO.LOC_BOUNDARY, GO.LOC_BOUNDARY] == GO.DIM_P
    @test DE9IM() == DE9IM("FFFFFFFFF")
    im2 = GO.with_entry(DE9IM(), GO.LOC_INTERIOR, GO.LOC_BOUNDARY, GO.DIM_L)
    @test string(im2) == "F1FFFFFFF"
    # pattern matching (JTS IntersectionMatrix.matches semantics)
    @test GO.matches(DE9IM("212101212"), "T*F**FFF*") == false
    @test GO.matches(DE9IM("2FF1FF212"), "T*F**FFF*") == false
    @test GO.matches(DE9IM("2FF1FFFF2"), "T*F**FFF*") == true
    @test GO.matches(DE9IM("0FFFFFFF2"), "0FFFFFFF2") == true
    @test_throws ArgumentError DE9IM("212")        # wrong length
    @test_throws ArgumentError GO.matches(DE9IM(), "T*F**FF")  # wrong length
end
```

**Step 2: Run to verify failure**

```bash
julia --project=docs test/methods/relateng/de9im.jl
```

Expected: FAIL / error — `LOC_INTERIOR not defined`.

**Step 3: Implement** — `src/methods/geom_relations/relateng/de9im.jl`:

```julia
# # DE-9IM matrix, location and dimension codes
#=
Port of code-level concepts from JTS: `Location`, `Dimension`,
`IntersectionMatrix` (pattern matching only), and
`operation/relateng/DimensionLocation.java`.
The `DE9IM` struct is immutable and isbits (`NTuple{9, Int8}`).
=#

# Locations (JTS Location.java)
const LOC_INTERIOR = Int8(0)
const LOC_BOUNDARY = Int8(1)
const LOC_EXTERIOR = Int8(2)
const LOC_NONE     = Int8(-1)

# Dimensions (JTS Dimension.java)
const DIM_FALSE    = Int8(-1)   # 'F'
const DIM_TRUE     = Int8(-2)   # 'T'  (patterns only)
const DIM_DONTCARE = Int8(-3)   # '*'  (patterns only)
const DIM_P = Int8(0)
const DIM_L = Int8(1)
const DIM_A = Int8(2)

function dim_char(d::Integer)
    d == DIM_FALSE    && return 'F'
    d == DIM_TRUE     && return 'T'
    d == DIM_DONTCARE && return '*'
    (0 <= d <= 2)     && return Char('0' + d)
    throw(ArgumentError("invalid dimension code $d"))
end

function dim_code(c::AbstractChar)
    c == 'F' && return DIM_FALSE
    c == 'T' && return DIM_TRUE
    c == '*' && return DIM_DONTCARE
    ('0' <= c <= '2') && return Int8(c - '0')
    throw(ArgumentError("invalid DE-9IM character '$c'"))
end

# DimensionLocation codes (JTS DimensionLocation.java, verbatim values)
const DL_EXTERIOR       = LOC_EXTERIOR  # 2
const DL_POINT_INTERIOR = Int8(103)
const DL_LINE_INTERIOR  = Int8(110)
const DL_LINE_BOUNDARY  = Int8(111)
const DL_AREA_INTERIOR  = Int8(120)
const DL_AREA_BOUNDARY  = Int8(121)

dimloc_point(loc::Integer) = loc == LOC_INTERIOR ? DL_POINT_INTERIOR : DL_EXTERIOR
function dimloc_line(loc::Integer)
    loc == LOC_INTERIOR && return DL_LINE_INTERIOR
    loc == LOC_BOUNDARY && return DL_LINE_BOUNDARY
    return DL_EXTERIOR
end
function dimloc_area(loc::Integer)
    loc == LOC_INTERIOR && return DL_AREA_INTERIOR
    loc == LOC_BOUNDARY && return DL_AREA_BOUNDARY
    return DL_EXTERIOR
end
dimloc_location(dimloc::Integer) = dimloc < 100 ? Int8(dimloc) : Int8(dimloc % 10)
function dimloc_dimension(dimloc::Integer, exterior_dim::Integer = DIM_FALSE)
    dimloc < 100 && return Int8(exterior_dim)
    return Int8(dimloc ÷ 10 - 10)
end

"""
    DE9IM

An immutable DE-9IM intersection matrix. Entries are dimension codes
(`DIM_FALSE`, `DIM_P`, `DIM_L`, `DIM_A`) stored row-major over
(Interior, Boundary, Exterior) of A × B, matching the standard string
form `"212101212"`. Construct from a 9-character string or empty
(all-`F`) via `DE9IM()`. Index with `im[locA, locB]`.
"""
struct DE9IM
    entries::NTuple{9, Int8}
end
DE9IM() = DE9IM(ntuple(_ -> DIM_FALSE, 9))
function DE9IM(s::AbstractString)
    length(s) == 9 || throw(ArgumentError("DE-9IM string must have 9 characters, got $(repr(s))"))
    return DE9IM(ntuple(i -> dim_code(s[i]), 9))
end

@inline im_index(locA::Integer, locB::Integer) = 3 * Int(locA) + Int(locB) + 1
Base.getindex(im::DE9IM, locA::Integer, locB::Integer) = im.entries[im_index(locA, locB)]
with_entry(im::DE9IM, locA::Integer, locB::Integer, dim::Integer) =
    DE9IM(Base.setindex(im.entries, Int8(dim), im_index(locA, locB)))

Base.string(im::DE9IM) = join(dim_char(d) for d in im.entries)
Base.show(io::IO, im::DE9IM) = print(io, "DE9IM(\"", string(im), "\")")

"Match a single matrix entry against a pattern code (JTS `IntersectionMatrix.matches`)."
function matches_entry(dim::Int8, pat::Int8)
    pat == DIM_DONTCARE && return true
    pat == DIM_TRUE     && return dim >= 0
    return dim == pat
end

function matches(im::DE9IM, pattern::AbstractString)
    length(pattern) == 9 || throw(ArgumentError("DE-9IM pattern must have 9 characters, got $(repr(pattern))"))
    return all(matches_entry(im.entries[i], dim_code(pattern[i])) for i in 1:9)
end
```

Add to `src/GeometryOps.jl` after the `common.jl` include:

```julia
include("methods/geom_relations/relateng/de9im.jl")
```

**Step 4: Run tests** — same command. Expected: all PASS.

**Step 5: Commit**

```bash
git add src/methods/geom_relations/relateng/de9im.jl src/GeometryOps.jl test/methods/relateng/de9im.jl
git commit -m "Add \`DE9IM\` matrix type and location/dimension codes for RelateNG"
```

## Task 2: Boundary node rules

**Files:**
- Modify: `src/methods/geom_relations/relateng/de9im.jl` (append; small enough not to need its own file)
- Test: `test/methods/relateng/de9im.jl` (append)

**Java reference:** `org/locationtech/jts/algorithm/BoundaryNodeRule.java`.

**Step 1: Failing tests** (append to `test/methods/relateng/de9im.jl`):

```julia
@testset "BoundaryNodeRule" begin
    @test GO.is_in_boundary(GO.Mod2Boundary(), 1) == true
    @test GO.is_in_boundary(GO.Mod2Boundary(), 2) == false
    @test GO.is_in_boundary(GO.Mod2Boundary(), 3) == true
    @test GO.is_in_boundary(GO.EndpointBoundary(), 1) == true
    @test GO.is_in_boundary(GO.EndpointBoundary(), 2) == true
    @test GO.is_in_boundary(GO.MultivalentEndpointBoundary(), 1) == false
    @test GO.is_in_boundary(GO.MultivalentEndpointBoundary(), 2) == true
    @test GO.is_in_boundary(GO.MonovalentEndpointBoundary(), 1) == true
    @test GO.is_in_boundary(GO.MonovalentEndpointBoundary(), 2) == false
end
```

**Step 2: Run** — expected FAIL (`Mod2Boundary not defined`).

**Step 3: Implement** (append to `de9im.jl`):

```julia
# Boundary node rules (JTS BoundaryNodeRule.java). Zero-size structs.
abstract type BoundaryNodeRule end
"OGC SFS standard rule: a vertex is on the boundary iff an odd number of line ends meet it (Mod-2 rule)."
struct Mod2Boundary <: BoundaryNodeRule end
struct EndpointBoundary <: BoundaryNodeRule end
struct MultivalentEndpointBoundary <: BoundaryNodeRule end
struct MonovalentEndpointBoundary <: BoundaryNodeRule end

is_in_boundary(::Mod2Boundary, boundary_count::Integer) = isodd(boundary_count)
is_in_boundary(::EndpointBoundary, boundary_count::Integer) = boundary_count > 0
is_in_boundary(::MultivalentEndpointBoundary, boundary_count::Integer) = boundary_count > 1
is_in_boundary(::MonovalentEndpointBoundary, boundary_count::Integer) = boundary_count == 1
```

**Step 4: Run** — PASS.

**Step 5: Commit** — `git commit -m "Add boundary node rule types for RelateNG"`

## Task 3: Predicate framework (`TopologyPredicate` API, tri-state, `BasicPredicate`, `IMPredicate` core)

**Files:**
- Create: `src/methods/geom_relations/relateng/topology_predicate.jl`
- Modify: `src/GeometryOps.jl` (include after `de9im.jl`)
- Test: `test/methods/relateng/predicates.jl` (+ register in `test/methods/relateng/runtests.jl`)

**Java references:** `JTS:TopologyPredicate.java` (166 lines — the full contract, read it first), `JTS:BasicPredicate.java` (108), `JTS:IMPredicate.java` (134).

**Design mapping (decided, do not re-derive):** Julia has no field inheritance, so JTS's class triangle becomes two kind-parameterized mutable structs. Per-kind behavior (`is_determined`, `value_im`, requirement flags, init overrides) dispatches on the kind singleton; requirement flags are pure functions of the *type*, so evaluation specializes per predicate (design D1 performance note).

**Step 1: Failing tests** — `test/methods/relateng/predicates.jl`. Exercise the framework through the two `BasicPredicate` kinds (implemented in this task) before any IM kinds exist:

```julia
using Test
import GeometryOps as GO
import Extents

@testset "tri-state and intersects predicate" begin
    p = GO.pred_intersects()
    @test GO.predicate_name(p) == "intersects"
    @test !GO.is_known(p)
    @test GO.require_self_noding(typeof(p)) == false
    @test GO.require_interaction(typeof(p)) == true
    @test GO.require_exterior_check(typeof(p), true) == false
    # interior/interior intersection determines intersects=true immediately
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_P)
    @test GO.is_known(p) && GO.predicate_value(p) == true
    # exterior-only updates never determine it; finish! defaults to false
    q = GO.pred_intersects()
    GO.update_dim!(q, GO.LOC_INTERIOR, GO.LOC_EXTERIOR, GO.DIM_L)
    @test !GO.is_known(q)
    GO.finish!(q)
    @test GO.is_known(q) && GO.predicate_value(q) == false
    # disjoint envelopes determine intersects=false at init_bounds!
    r = GO.pred_intersects()
    GO.init_bounds!(r, Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
                       Extents.Extent(X=(5.0, 6.0), Y=(5.0, 6.0)))
    @test GO.is_known(r) && GO.predicate_value(r) == false
end

@testset "disjoint predicate" begin
    p = GO.pred_disjoint()
    @test GO.require_interaction(typeof(p)) == false
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_P)
    @test GO.is_known(p) && GO.predicate_value(p) == false
    q = GO.pred_disjoint()
    GO.finish!(q)
    @test GO.is_known(q) && GO.predicate_value(q) == true
end
```

Register in `test/methods/relateng/runtests.jl` (uncomment/add the `Predicates` line).

**Step 2: Run** — FAIL (`pred_intersects not defined`).

**Step 3: Implement** — `src/methods/geom_relations/relateng/topology_predicate.jl`. Port `TopologyPredicate.java` + `BasicPredicate.java` + `IMPredicate.java` faithfully. Framework skeleton (complete — the `intersects`/`disjoint` kinds live here because they are `BasicPredicate`s; all IM kinds come in Task 4):

```julia
# # Topology predicate framework
#=
Port of JTS `TopologyPredicate`, `BasicPredicate`, `IMPredicate`
(operation/relateng). Each named predicate is a kind singleton
parameterizing one of two mutable state structs, so requirement flags
and determination checks are compile-time-specialized per predicate.
=#

abstract type TopologyPredicate end

# Tri-state (BasicPredicate.java)
const TRI_UNKNOWN = Int8(-1)
const TRI_FALSE   = Int8(0)
const TRI_TRUE    = Int8(1)

# --- API with JTS defaults (TopologyPredicate.java) ---
require_self_noding(::Type{<:TopologyPredicate}) = true
require_interaction(::Type{<:TopologyPredicate}) = true
require_covers(::Type{<:TopologyPredicate}, is_source_a::Bool) = false
require_exterior_check(::Type{<:TopologyPredicate}, is_source_a::Bool) = true
init_dims!(p::TopologyPredicate, dimA::Integer, dimB::Integer) = nothing
init_bounds!(p::TopologyPredicate, extA, extB) = nothing
# update_dim!(p, locA, locB, dim), finish!(p), is_known(p), predicate_value(p),
# predicate_name(p) are implemented per struct below.

is_intersection(locA::Integer, locB::Integer) =
    locA != LOC_EXTERIOR && locB != LOC_EXTERIOR

ext_intersects(extA, extB) = Extents.intersects(extA, extB)

# --- BasicPredicate kinds ---
struct IntersectsPred end
struct DisjointPred end

mutable struct BasicPredicate{K} <: TopologyPredicate
    kind::K
    value::Int8
end
BasicPredicate(kind) = BasicPredicate(kind, TRI_UNKNOWN)

is_known(p::BasicPredicate) = p.value != TRI_UNKNOWN
predicate_value(p::BasicPredicate) = p.value == TRI_TRUE
set_value!(p, v::Bool) = is_known(p) ? nothing : (p.value = v ? TRI_TRUE : TRI_FALSE; nothing)
set_value_if!(p, v::Bool, cond::Bool) = cond ? set_value!(p, v) : nothing
require!(p, cond::Bool) = cond ? nothing : set_value!(p, false)

# intersects (RelatePredicate.java intersects())
pred_intersects() = BasicPredicate(IntersectsPred())
predicate_name(::BasicPredicate{IntersectsPred}) = "intersects"
require_self_noding(::Type{BasicPredicate{IntersectsPred}}) = false
require_exterior_check(::Type{BasicPredicate{IntersectsPred}}, is_source_a::Bool) = false
init_bounds!(p::BasicPredicate{IntersectsPred}, extA, extB) =
    require!(p, ext_intersects(extA, extB))
update_dim!(p::BasicPredicate{IntersectsPred}, locA, locB, dim) =
    set_value_if!(p, true, is_intersection(locA, locB))
finish!(p::BasicPredicate{IntersectsPred}) = set_value!(p, false)

# disjoint (RelatePredicate.java disjoint())
pred_disjoint() = BasicPredicate(DisjointPred())
predicate_name(::BasicPredicate{DisjointPred}) = "disjoint"
require_self_noding(::Type{BasicPredicate{DisjointPred}}) = false
require_interaction(::Type{BasicPredicate{DisjointPred}}) = false
require_exterior_check(::Type{BasicPredicate{DisjointPred}}, is_source_a::Bool) = false
init_bounds!(p::BasicPredicate{DisjointPred}, extA, extB) =
    set_value_if!(p, true, !ext_intersects(extA, extB))
update_dim!(p::BasicPredicate{DisjointPred}, locA, locB, dim) =
    set_value_if!(p, false, is_intersection(locA, locB))
finish!(p::BasicPredicate{DisjointPred}) = set_value!(p, true)

# --- IMPredicate core (IMPredicate.java) ---
const DIM_UNKNOWN = DIM_DONTCARE   # JTS IMPredicate.DIM_UNKNOWN

mutable struct IMPredicate{K} <: TopologyPredicate
    kind::K
    dimA::Int8
    dimB::Int8
    im::DE9IM
    value::Int8
end
IMPredicate(kind) = IMPredicate(kind, DIM_UNKNOWN, DIM_UNKNOWN, DE9IM(), TRI_UNKNOWN)

is_known(p::IMPredicate) = p.value != TRI_UNKNOWN
predicate_value(p::IMPredicate) = p.value == TRI_TRUE

function init_dims!(p::IMPredicate, dimA::Integer, dimB::Integer)
    p.dimA = dimA; p.dimB = dimB
    init_dims_kind!(p)   # per-kind hook, default no-op
end
init_dims_kind!(p::IMPredicate) = nothing

is_dim_changed(p::IMPredicate, locA, locB, dim) = dim > p.im[locA, locB]

function update_dim!(p::IMPredicate, locA, locB, dim)
    if is_dim_changed(p, locA, locB, dim)
        p.im = with_entry(p.im, locA, locB, dim)
        if is_determined(p)   # per-kind
            p.value = value_im(p) ? TRI_TRUE : TRI_FALSE
        end
    end
    nothing
end

function finish!(p::IMPredicate)
    is_known(p) && return nothing
    p.value = value_im(p) ? TRI_TRUE : TRI_FALSE
    nothing
end

# Shared IM queries (IMPredicate.java helpers)
intersects_exterior_of(p::IMPredicate, is_a::Bool) = is_a ?
    (is_intersects_entry(p, LOC_EXTERIOR, LOC_INTERIOR) || is_intersects_entry(p, LOC_EXTERIOR, LOC_BOUNDARY)) :
    (is_intersects_entry(p, LOC_INTERIOR, LOC_EXTERIOR) || is_intersects_entry(p, LOC_BOUNDARY, LOC_EXTERIOR))
is_intersects_entry(p::IMPredicate, locA, locB) = p.im[locA, locB] >= DIM_P
is_known_entry(p::IMPredicate, locA, locB) = p.im[locA, locB] != DIM_FALSE  # verify vs Java
is_dimension_entry(p::IMPredicate, locA, locB, dim) = p.im[locA, locB] == dim
get_dimension(p::IMPredicate, locA, locB) = p.im[locA, locB]
```

**Port note:** `is_known_entry` semantics must be verified against `IMPredicate.java` (JTS tracks "known" differently from `FALSE` — it initializes entries to `DIM_UNKNOWN`, not `FALSE`; if so, initialize `p.im` entries to `DIM_UNKNOWN` and adjust `finish!` to map unknown→`F`. Follow the Java exactly.)

**Step 4: Run** — PASS (fix until green).

**Step 5: Commit** — `Add topology predicate framework for RelateNG`

## Task 4: Named IM predicates, pattern matcher, matrix predicate

**Files:**
- Create: `src/methods/geom_relations/relateng/relate_predicates.jl`
- Modify: `src/GeometryOps.jl`
- Test: `test/methods/relateng/predicates.jl` (append)

**Java references:** `JTS:RelatePredicate.java` (631 lines — contains all 8 IM kinds with their flag overrides, `init` overrides, `isDetermined`, `valueIM`), `JTS:IMPatternMatcher.java` (111), `JTS:IntersectionMatrixPattern.java` (63), `JTS:RelateMatrixPredicate.java` (55). Test source: `RelatePredicateTest.java` (92 lines).

**Step 1: Port `RelatePredicateTest.java` to Julia** — append to `test/methods/relateng/predicates.jl`. Port every test method; they drive predicates via `update_dim!`/`finish!` sequences and check values. Also add flag-table assertions for each kind:

| factory | self_noding | interaction | covers(A)/covers(B) | ext_check(A)/ext_check(B) |
|---|---|---|---|---|
| `pred_contains()` | per Java | T | per Java | per Java |
| `pred_within()` | per Java | T | per Java | per Java |
| `pred_covers()` | per Java | T | per Java | per Java |
| `pred_coveredby()` | per Java | T | per Java | per Java |
| `pred_crosses()` | per Java | T | F/F | per Java |
| `pred_equalstopo()` | per Java | T | per Java | per Java |
| `pred_overlaps()` | per Java | T | F/F | per Java |
| `pred_touches()` | per Java | T | F/F | per Java |

Fill the "per Java" cells while reading `RelatePredicate.java` — write the assertion for the value the Java declares (each inner class's `requireX` overrides). Do not guess.

Add pattern-matcher tests:

```julia
@testset "IMPatternMatcher" begin
    p = GO.IMPatternMatcher("T*F**FFF*")
    @test GO.predicate_name(p) == "IMPattern"
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_A)
    @test !GO.is_known(p)
    # an entry violating the pattern mask determines false immediately
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_EXTERIOR, GO.DIM_L)  # pattern 'F' at I/E
    @test GO.is_known(p) && GO.predicate_value(p) == false
end

@testset "RelateMatrixPredicate" begin
    p = GO.RelateMatrixPredicate()
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_A)
    GO.finish!(p)
    @test !isnothing(GO.result_im(p))
    @test GO.result_im(p)[GO.LOC_INTERIOR, GO.LOC_INTERIOR] == GO.DIM_A
end
```

**Step 2: Run** — FAIL.

**Step 3: Implement** — `relate_predicates.jl`. One kind singleton + `IMPredicate{K}` method set per JTS inner class, in the same order as `RelatePredicate.java`: `ContainsPred`, `WithinPred`, `CoversPred`, `CoveredByPred`, `CrossesPred`, `EqualsTopoPred`, `OverlapsPred`, `TouchesPred`. For each, port: factory function, `predicate_name`, requirement-flag overrides, `init_dims_kind!` (dimension-compatibility early exits, e.g. contains is false if `dimA < dimB`), `init_bounds!` (`require_covers` envelope checks via `ext_covers(extA, extB)` — implement `ext_covers` with explicit interval comparisons if `Extents` lacks a `contains`), `is_determined`, `value_im`. Then port `IMPatternMatcher` (standalone mutable struct holding `pattern::String`, `pattern_matrix::DE9IM`, plus IM state — including its `require_interaction` computed from the pattern) and `RelateMatrixPredicate` (never determined early; `result_im(p)` accessor returns the accumulated `DE9IM`).

**Step 4: Run** — PASS.

**Step 5: Commit** — `Add named DE-9IM predicates and pattern matcher for RelateNG`

---

# Stage 2 — The exact geometry kernel

## Task 5: Kernel API + planar orientation/on-segment/point-in-ring/bounds

**Files:**
- Create: `src/methods/geom_relations/relateng/kernel.jl` (API contract, docstrings, generic fallbacks)
- Create: `src/methods/geom_relations/relateng/kernel_planar.jl`
- Modify: `src/GeometryOps.jl` (include both, *before* the topology-layer files that will use them)
- Test: `test/methods/relateng/kernel.jl`

**Step 1: Failing tests** — `test/methods/relateng/kernel.jl`:

```julia
using Test
import GeometryOps as GO
import GeometryOps: Planar, True, False
import GeoInterface as GI

const PT = Tuple{Float64, Float64}
m = Planar()

@testset "rk_orient" begin
    @test GO.rk_orient(m, (0.0,0.0), (1.0,0.0), (0.0,1.0); exact = True()) > 0
    @test GO.rk_orient(m, (0.0,0.0), (1.0,0.0), (0.0,-1.0); exact = True()) < 0
    @test GO.rk_orient(m, (0.0,0.0), (1.0,0.0), (2.0,0.0); exact = True()) == 0
    # adversarial near-collinear: exact must get the sign right
    a, b = (0.0, 0.0), (1.0, 1.0)
    c = (0.5, 0.5 + 1e-17)   # above the line by less than eps
    @test GO.rk_orient(m, a, b, c; exact = True()) == GO.rk_orient(m, a, b, (0.5, 0.6); exact = True())
end

@testset "rk_point_on_segment" begin
    @test GO.rk_point_on_segment(m, (0.5,0.5), (0.0,0.0), (1.0,1.0); exact = True()) == true
    @test GO.rk_point_on_segment(m, (2.0,2.0), (0.0,0.0), (1.0,1.0); exact = True()) == false   # collinear, outside
    @test GO.rk_point_on_segment(m, (0.5,0.6), (0.0,0.0), (1.0,1.0); exact = True()) == false
    @test GO.rk_point_on_segment(m, (0.0,0.0), (0.0,0.0), (1.0,1.0); exact = True()) == true    # endpoint inclusive
end

@testset "rk_point_in_ring" begin
    ring = GI.LinearRing([(0.0,0.0), (10.0,0.0), (10.0,10.0), (0.0,10.0), (0.0,0.0)])
    @test GO.rk_point_in_ring(m, (5.0,5.0), ring; exact = True()) == GO.LOC_INTERIOR
    @test GO.rk_point_in_ring(m, (5.0,0.0), ring; exact = True()) == GO.LOC_BOUNDARY
    @test GO.rk_point_in_ring(m, (15.0,5.0), ring; exact = True()) == GO.LOC_EXTERIOR
end

@testset "bounds" begin
    pa = GI.Polygon([[(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,0.0)]])
    ea = GO.rk_interaction_bounds(m, pa)
    @test !GO.rk_bounds_disjoint(ea, ea)
    @test GO.rk_bounds_covers(ea, ea)
end
```

Register `kernel.jl` testset in `test/methods/relateng/runtests.jl`.

**Step 2: Run** — FAIL.

**Step 3: Implement.** `kernel.jl` declares the contract: a docstring block listing every kernel function (this is the spherical-implementation spec) and any manifold-generic helpers. `kernel_planar.jl`:

```julia
# # Planar RelateKernel
#=
The geometry layer of RelateNG (design doc D1/D2): every coordinate-level
question the topology layer may ask, answered with exact predicates and
no constructed coordinates. Planar implementation; a future Spherical
kernel implements the same functions and must pass the same
conformance testset.
=#

rk_orient(::Planar, a, b, c; exact) = Predicates.orient(a, b, c; exact)

# valid only when p is known collinear with (q0, q1)
@inline function _collinear_between(p, q0, q1)
    (min(GI.x(q0), GI.x(q1)) <= GI.x(p) <= max(GI.x(q0), GI.x(q1))) &&
    (min(GI.y(q0), GI.y(q1)) <= GI.y(p) <= max(GI.y(q0), GI.y(q1)))
end

function rk_point_on_segment(m::Planar, p, q0, q1; exact)
    rk_orient(m, q0, q1, p; exact) == 0 || return false
    return _collinear_between(p, q0, q1)
end

function rk_point_in_ring(m::Planar, p, ring; exact)
    o = _point_filled_curve_orientation(m, p, ring; exact)
    o == point_in  && return LOC_INTERIOR
    o == point_on  && return LOC_BOUNDARY
    return LOC_EXTERIOR
end

rk_interaction_bounds(::Planar, geom) = GI.extent(geom, fallback = true)
rk_bounds_disjoint(extA, extB) = !Extents.intersects(extA, extB)
function rk_bounds_covers(extA, extB)
    (extA.X[1] <= extB.X[1] && extB.X[2] <= extA.X[2]) &&
    (extA.Y[1] <= extB.Y[1] && extB.Y[2] <= extA.Y[2])
end
```

Notes: `_point_filled_curve_orientation` is the existing Hao–Sun machinery at `src/methods/geom_relations/geom_geom_processors.jl:495` (signature `(::Planar, point, curve; in, on, out, exact)`); `GI.extent(geom, fallback=true)` — check actual extent-fetch idiom used elsewhere in GO (e.g. clipping) and match it.

**Step 4: Run** — PASS.

**Step 5: Commit** — `Add planar RelateKernel orientation, on-segment, point-in-ring and bounds`

## Task 6: Symbolic segment-intersection classification

**Files:**
- Modify: `src/methods/geom_relations/relateng/kernel.jl` (the `SegSegClass` type is manifold-generic), `kernel_planar.jl`
- Test: `test/methods/relateng/kernel.jl` (append)

This replaces JTS's `RobustLineIntersector`: classification only, no intersection coordinates (design D2).

**Step 1: Failing tests** (append):

```julia
@testset "rk_classify_intersection" begin
    cl(a0,a1,b0,b1) = GO.rk_classify_intersection(m, a0, a1, b0, b1; exact = True())
    # disjoint
    @test cl((0.,0.),(1.,0.),(0.,1.),(1.,1.)).kind == GO.SS_DISJOINT
    # proper crossing
    r = cl((0.,0.),(2.,2.),(0.,2.),(2.,0.))
    @test r.kind == GO.SS_PROPER
    @test !(r.a0_on_b || r.a1_on_b || r.b0_on_a || r.b1_on_a)
    # touch: b0 on interior of a
    r = cl((0.,0.),(2.,0.),(1.,0.),(1.,1.))
    @test r.kind == GO.SS_TOUCH && r.b0_on_a && !r.b1_on_a && !r.a0_on_b && !r.a1_on_b
    # touch: shared endpoint
    r = cl((0.,0.),(1.,0.),(1.,0.),(1.,1.))
    @test r.kind == GO.SS_TOUCH && r.a1_on_b && r.b0_on_a
    # collinear overlap
    r = cl((0.,0.),(2.,0.),(1.,0.),(3.,0.))
    @test r.kind == GO.SS_COLLINEAR && r.b0_on_a && r.a1_on_b
    # collinear disjoint
    @test cl((0.,0.),(1.,0.),(2.,0.),(3.,0.)).kind == GO.SS_DISJOINT
    # collinear, touching only at one shared endpoint -> SS_TOUCH
    r = cl((0.,0.),(1.,0.),(1.,0.),(2.,0.))
    @test r.kind == GO.SS_TOUCH && r.a1_on_b && r.b0_on_a
    # containment: b inside a (collinear)
    r = cl((0.,0.),(3.,0.),(1.,0.),(2.,0.))
    @test r.kind == GO.SS_COLLINEAR && r.b0_on_a && r.b1_on_a
    # degenerate: zero-length b on a
    r = cl((0.,0.),(2.,0.),(1.,0.),(1.,0.))
    @test r.kind == GO.SS_TOUCH && r.b0_on_a && r.b1_on_a
end
```

**Step 2: Run** — FAIL.

**Step 3: Implement.** In `kernel.jl`:

```julia
# Symbolic segment-pair intersection classification (replaces RobustLineIntersector).
@enum SegSegKind::Int8 SS_DISJOINT SS_PROPER SS_TOUCH SS_COLLINEAR

"""
    SegSegClass

Combinatorial classification of the intersection of closed segments
(a0,a1) × (b0,b1). `kind` is `SS_PROPER` only for a crossing in both
segments' interiors (the node is *symbolic*: no coordinate exists for
it anywhere in the engine). All vertex incidences are reported via the
`*_on_*` flags, whose coordinates are exact input vertices.
"""
struct SegSegClass
    kind::SegSegKind
    a0_on_b::Bool
    a1_on_b::Bool
    b0_on_a::Bool
    b1_on_a::Bool
end
```

In `kernel_planar.jl`:

```julia
function rk_classify_intersection(m::Planar, a0, a1, b0, b1; exact)
    oa0 = rk_orient(m, b0, b1, a0; exact)
    oa1 = rk_orient(m, b0, b1, a1; exact)
    ob0 = rk_orient(m, a0, a1, b0; exact)
    ob1 = rk_orient(m, a0, a1, b1; exact)
    # fully collinear configuration (handles zero-length segments too)
    if oa0 == 0 && oa1 == 0 && ob0 == 0 && ob1 == 0
        a0_on_b = _collinear_between(a0, b0, b1)
        a1_on_b = _collinear_between(a1, b0, b1)
        b0_on_a = _collinear_between(b0, a0, a1)
        b1_on_a = _collinear_between(b1, a0, a1)
        n_inc = a0_on_b + a1_on_b + b0_on_a + b1_on_a
        n_inc == 0 && return SegSegClass(SS_DISJOINT, false, false, false, false)
        # single shared endpoint counts twice (one endpoint of each on the other)
        shared_endpoint_only = n_inc == 2 &&
            ((a0_on_b || a1_on_b) && (b0_on_a || b1_on_a)) &&
            (_equals2(a0, b0) || _equals2(a0, b1) || _equals2(a1, b0) || _equals2(a1, b1))
        kind = shared_endpoint_only ? SS_TOUCH : SS_COLLINEAR
        # zero-length degenerate: a point on a segment is a touch, not an overlap
        if _equals2(a0, a1) || _equals2(b0, b1)
            kind = SS_TOUCH
        end
        return SegSegClass(kind, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    end
    a0_on_b = oa0 == 0 && _collinear_between(a0, b0, b1)
    a1_on_b = oa1 == 0 && _collinear_between(a1, b0, b1)
    b0_on_a = ob0 == 0 && _collinear_between(b0, a0, a1)
    b1_on_a = ob1 == 0 && _collinear_between(b1, a0, a1)
    if a0_on_b || a1_on_b || b0_on_a || b1_on_a
        return SegSegClass(SS_TOUCH, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    end
    if (oa0 > 0) != (oa1 > 0) && oa0 != 0 && oa1 != 0 &&
       (ob0 > 0) != (ob1 > 0) && ob0 != 0 && ob1 != 0
        return SegSegClass(SS_PROPER, false, false, false, false)
    end
    return SegSegClass(SS_DISJOINT, false, false, false, false)
end

_equals2(p, q) = GI.x(p) == GI.x(q) && GI.y(p) == GI.y(q)
```

**Step 4: Run** — PASS. Add any edge case that failed during implementation as a permanent test.

**Step 5: Commit** — `Add exact symbolic segment intersection classification to RelateKernel`

## Task 7: Symbolic node identity and the rational coincidence slow path

**Files:**
- Modify: `kernel.jl` (NodeKey type), `kernel_planar.jl` (rational path)
- Test: `test/methods/relateng/kernel.jl` (append)

Design D2/D3: vertex nodes key exactly by coordinate; proper crossings key by canonicalized segment pair; cross-kind coincidence uses exact rational arithmetic and is only invoked on self-noding paths.

**Step 1: Failing tests** (append):

```julia
@testset "NodeKey" begin
    v = GO.vertex_node((1.0, 2.0))
    v2 = GO.vertex_node((1.0, 2.0))
    @test v == v2 && hash(v) == hash(v2)
    c1 = GO.crossing_node((0.,0.), (2.,2.), (0.,2.), (2.,0.))
    c2 = GO.crossing_node((0.,2.), (2.,0.), (2.,2.), (0.,0.))  # same pair, swapped & reversed
    @test c1 == c2 && hash(c1) == hash(c2)
    @test v != c1
end

@testset "exact crossing coincidence (rational slow path)" begin
    # X crossing at exactly (1,1); a vertex node placed there must coincide
    c = GO.crossing_node((0.,0.), (2.,2.), (0.,2.), (2.,0.))
    @test GO.rk_nodes_coincide(m, c, GO.vertex_node((1.0, 1.0)); exact = True()) == true
    @test GO.rk_nodes_coincide(m, c, GO.vertex_node((1.0, 1.0 + eps(1.0))); exact = True()) == false
    # two crossings meeting at the same point
    c2 = GO.crossing_node((1.,0.), (1.,2.), (0.,1.), (2.,1.))
    @test GO.rk_nodes_coincide(m, c, c2; exact = True()) == true
    # crossing point with non-representable rational coordinates
    c3 = GO.crossing_node((0.,0.), (3.,1.), (0.,1.), (3.,0.))  # crosses at (1.5, 0.5)
    @test GO.rk_nodes_coincide(m, c3, GO.vertex_node((1.5, 0.5)); exact = True()) == true
end
```

**Step 2: Run** — FAIL.

**Step 3: Implement.** In `kernel.jl`:

```julia
# Symbolic node identity (design D2). One concrete isbits key type for both
# node kinds so Dict{NodeKey{P}, ...} is type-stable.
struct NodeKey{P}
    is_crossing::Bool
    pt::P          # vertex nodes: the coordinate. crossing nodes: canonical a0.
    a1::P
    b0::P
    b1::P
end

vertex_node(pt) = NodeKey(false, pt, pt, pt, pt)

"Canonicalize: each segment ordered lexicographically by (x, y); segments ordered by their first point."
function crossing_node(a0, a1, b0, b1)
    a0, a1 = _seg_canon(a0, a1)
    b0, b1 = _seg_canon(b0, b1)
    if (GI.x(b0), GI.y(b0), GI.x(b1), GI.y(b1)) < (GI.x(a0), GI.y(a0), GI.x(a1), GI.y(a1))
        a0, a1, b0, b1 = b0, b1, a0, a1
    end
    return NodeKey(true, a0, a1, b0, b1)
end
_seg_canon(p, q) = (GI.x(p), GI.y(p)) <= (GI.x(q), GI.y(q)) ? (p, q) : (q, p)
```

(`==`/`hash` come free from the default struct definitions since fields are tuples — add explicit ones only if `P` is not bitstype-comparable.)

In `kernel_planar.jl`, the rational slow path (Float64 values are exact rationals, so `Rational{BigInt}` arithmetic is exact):

```julia
"Exact intersection point of two properly crossing segments, as rationals."
function _exact_crossing_point(a0, a1, b0, b1)
    R = Rational{BigInt}
    ax0, ay0 = R(GI.x(a0)), R(GI.y(a0)); ax1, ay1 = R(GI.x(a1)), R(GI.y(a1))
    bx0, by0 = R(GI.x(b0)), R(GI.y(b0)); bx1, by1 = R(GI.x(b1)), R(GI.y(b1))
    dax, day = ax1 - ax0, ay1 - ay0
    dbx, dby = bx1 - bx0, by1 - by0
    denom = dax * dby - day * dbx          # nonzero for a proper crossing
    t = ((bx0 - ax0) * dby - (by0 - ay0) * dbx) // denom
    return (ax0 + t * dax, ay0 + t * day)
end

_exact_node_point(k::NodeKey) = k.is_crossing ?
    _exact_crossing_point(k.pt, k.a1, k.b0, k.b1) :
    (Rational{BigInt}(GI.x(k.pt)), Rational{BigInt}(GI.y(k.pt)))

function rk_nodes_coincide(::Planar, k1::NodeKey, k2::NodeKey; exact)
    k1 == k2 && return true
    # Slow path (design D3, follow-up F1): exact rational comparison.
    return _exact_node_point(k1) == _exact_node_point(k2)
end
```

**Step 4: Run** — PASS.

**Step 5: Commit** — `Add symbolic node identity with exact rational coincidence slow path`

## Task 8: Edge ordering around nodes (`PolygonNodeTopology` port with symbolic apex)

**Files:**
- Modify: `kernel.jl`, `kernel_planar.jl`
- Test: `test/methods/relateng/kernel.jl` (append)

**Java reference:** `/Users/anshul/temp/GO_jts/jts/modules/core/src/main/java/org/locationtech/jts/algorithm/PolygonNodeTopology.java` — methods `compareAngle`, `isAngleGreater`, `isBetween`, `compareBetween`, `quadrant`, `isCrossing`, `isInteriorSegment`. Read in full.

**What changes vs Java:** the apex argument becomes a `NodeKey`. For vertex nodes the apex coordinate is exact and the port is direct. For crossing nodes the four incident directions are the segment endpoints themselves, and their CCW cyclic order is derived from orientation signs of the *original endpoints* — never a constructed apex.

**Step 1: Failing tests** (append):

```julia
@testset "edge ordering around a vertex node" begin
    origin = GO.vertex_node((0.0, 0.0))
    east, north, west, south = (1.,0.), (0.,1.), (-1.,0.), (0.,-1.)
    cmp(p, q) = GO.rk_compare_edge_dir(m, origin, p, q; exact = True())
    @test cmp(east, north) < 0   # JTS compareAngle: CCW from positive x-axis
    @test cmp(north, west) < 0
    @test cmp(west, south) < 0
    @test cmp(east, east) == 0
    @test cmp(north, east) > 0
    # same quadrant resolved by orientation
    @test cmp((2.0, 1.0), (1.0, 2.0)) < 0
end

@testset "crossing-node incident edge order" begin
    # a: (0,0)->(2,2), b: (0,2)->(2,0); crossing at symbolic (1,1)
    dirs = GO.rk_crossing_dirs_ccw(m, (0.,0.), (2.,2.), (0.,2.), (2.,0.); exact = True())
    # CCW order starting from direction toward a1=(2,2):
    @test dirs == ((2.,2.), (0.,2.), (0.,0.), (2.,0.))
end

@testset "isCrossing / isInteriorSegment" begin
    n = (1.0, 1.0)
    @test GO.rk_is_crossing(m, GO.vertex_node(n), (0.,0.), (2.,2.), (0.,2.), (2.,0.); exact = True())
    @test !GO.rk_is_crossing(m, GO.vertex_node(n), (0.,0.), (2.,2.), (2.,0.), (2.,2.); exact = True())
end
```

(Adjust the `compareAngle` sign expectations to whatever the Java actually returns — write the test by reading `PolygonNodeTopology.compareAngle`'s contract first, then assert that.)

**Step 2: Run** — FAIL.

**Step 3: Implement.**

- `rk_quadrant(m, origin_pt, p)` — port `quadrant` (uses coordinate comparisons against the origin point).
- `rk_compare_edge_dir(m, node::NodeKey, p, q; exact)` — for `!node.is_crossing`, port `compareAngle(origin, p, q)` verbatim with `origin = node.pt`, routing the orientation test through `rk_orient`.
- For crossing nodes, `rk_compare_edge_dir` is only ever needed among the 4 incident endpoints; implement via `rk_crossing_dirs_ccw`:

```julia
"""
CCW cyclic order of the four half-edge directions incident to the
proper crossing of (a0,a1) × (b0,b1), starting from a1. Since the
crossing is proper, b0/b1 are strictly on opposite sides of line(a0,a1):
if b1 is to the left, CCW order is (a1, b1, a0, b0), else (a1, b0, a0, b1).
"""
function rk_crossing_dirs_ccw(m::Planar, a0, a1, b0, b1; exact)
    if rk_orient(m, a0, a1, b1; exact) > 0
        return (a1, b1, a0, b0)
    else
        return (a1, b0, a0, b1)
    end
end
```

- `rk_is_crossing(m, node, a0, a1, b0, b1; exact)` — port `PolygonNodeTopology.isCrossing` (apex = `node` vertex; used by `TopologyComputer.updateAreaAreaCross`).
- `rk_is_interior_segment(m, node, a0, a1, b; exact)` — port `isInteriorSegment`.

**Step 4: Run** — PASS.

**Step 5: Commit** — `Add edge ordering around symbolic nodes to RelateKernel`

## Task 9: Kernel conformance testset

**Files:**
- Create: `test/methods/relateng/kernel_conformance.jl`
- Test registration: `test/methods/relateng/runtests.jl`

This is the spec the future `Spherical` kernel must pass (design layer contract). Structure it as a function over a manifold:

**Step 1: Write the testset** — a `function kernel_conformance_suite(m; exact)` containing property-style checks, instantiated for `Planar()`:

- `rk_orient` antisymmetry (`orient(a,b,c) == -orient(b,a,c)`), cyclic invariance, degeneracy (`orient(a,a,b) == 0`).
- `rk_classify_intersection` symmetry: swapping A and B swaps the flag pairs; reversing a segment's endpoints swaps its two flags; classification of shared-endpoint configurations is `SS_TOUCH`.
- Consistency: `SS_PROPER` implies all four incidence flags false and all four orientations nonzero; every flagged endpoint passes `rk_point_on_segment`.
- `rk_compare_edge_dir` is a strict weak order on a fan of 16 directions around a vertex node.
- `rk_nodes_coincide` is reflexive/symmetric and agrees with `==` on equal keys.
- `rk_point_in_ring` agrees with `rk_point_on_segment` for points on ring edges.

Use randomized inputs from a fixed-seed RNG (`StableRNGs` if already a test dep, else `Random.seed!`) plus the explicit corner cases above.

**Step 2: Run, fix any kernel bugs it finds, re-run** — PASS.

**Step 3: Commit** — `Add RelateKernel conformance testset`

---

# Stage 3 — Point location

## Task 10: `LinearBoundary`

**Files:**
- Create: `src/methods/geom_relations/relateng/point_locator.jl` (holds LinearBoundary + AdjacentEdgeLocator + RelatePointLocator, in that order — they are small and tightly coupled; JTS file boundaries preserved as clearly-marked sections)
- Modify: `src/GeometryOps.jl`
- Test: `test/methods/relateng/point_locator.jl`

**Java references:** `JTS:LinearBoundary.java` (83 lines); tests from `LinearBoundaryTest.java` (97 lines — port every method; cases cover Mod2 vs Endpoint rules on lines/multilines).

**Step 1: Port `LinearBoundaryTest.java`** to `test/methods/relateng/point_locator.jl`, using `GO.tuples` on WKT strings via the existing test idiom (parse with `jts_wkt_to_geom`-style helper or construct GI geometries directly to avoid the dependency — prefer direct GI construction for unit tests).

**Step 2: Run** — FAIL.

**Step 3: Implement:**

```julia
struct LinearBoundary{BR <: BoundaryNodeRule, P}
    vertex_degree::Dict{P, Int}
    has_boundary::Bool
    rule::BR
end
```

Constructor takes an iterable of linestrings (each a coordinate vector) + rule; port `computeBoundaryPoints` (count endpoint degree per coordinate; closed lines contribute nothing per JTS — verify in Java), `hasBoundary`, `isBoundary(lb, pt)`.

**Step 4: Run** — PASS. **Step 5: Commit** — `Add LinearBoundary for RelateNG point location`

## Task 11: `AdjacentEdgeLocator`

**Files:** same source file; test file appended.

**Java references:** `JTS:AdjacentEdgeLocator.java` (117); `AdjacentEdgeLocatorTest.java` (85 — port all).

Depends on `NodeSection`/`NodeSections` only lightly in Java (it builds sections to test adjacency); port the dependency-minimal version: the Java's `addSections`/`createSection` usage can be expressed with a local struct or by moving `NodeSection` forward — **decision: create `node_sections.jl` with the plain `NodeSection` struct in this task** (fields only, no `NodeSections` logic yet), include it before `point_locator.jl`. Port `locate(ael, p)` returning `LOC_BOUNDARY`/`LOC_INTERIOR` with the section-based exterior-edge check, routing ring orientation through existing GO ring-orientation utilities and all geometry through the kernel.

TDD steps as usual. **Commit** — `Add AdjacentEdgeLocator for polygon union semantics`

## Task 12: `RelatePointLocator`

**Files:** same source file; test appended.

**Java references:** `JTS:RelatePointLocator.java` (347); `RelatePointLocatorTest.java` (101 — port all; it exercises mixed GeometryCollections).

Port: constructor extracts points/lines/polygons from the (possibly GC) input via `apply`/`GI.getgeom` traversal mirroring `extractElements`; methods `has_boundary`, `locate(p)`, `locate_with_dim(p)`, `locate_line_end_with_dim(p)`, `locate_node(p, parent_polygonal)`, `locate_node_with_dim`, and the private precedence chain `compute_dim_location` → `locate_on_points/lines/line/polygons/polygonal` returning `DL_*` codes. Point-on-line checks use `rk_point_on_segment` looped over segments; polygon location uses `rk_point_in_ring` over shell/holes (exterior ring then holes, standard evenodd composition — port `locateOnPolygonal`'s use of JTS `PolygonUtil`/locators as a loop over rings). `isNode` handling and the `AdjacentEdgeLocator` delegation ported exactly.

TDD steps as usual. **Commit** — `Add RelatePointLocator`

## Task 13: `RelateGeometry`

**Files:**
- Create: `src/methods/geom_relations/relateng/relate_geometry.jl`
- Modify: `src/GeometryOps.jl`
- Test: `test/methods/relateng/relate_geometry.jl`

**Java references:** `JTS:RelateGeometry.java` (420); `RelateGeometryTest.java` (72 — port all). Also `JTS:RelateSegmentString.java` (158) — port here or in Stage 5 Task 18; **decision: port `RelateSegmentString` in this task** since `extract_segment_strings` is `RelateGeometry` API.

Struct sketch:

```julia
mutable struct RelateGeometry{G, M, E, BR}
    geom::G
    manifold::M
    exact::E
    boundary_rule::BR
    is_prepared::Bool
    extent::Any            # concrete Extents.Extent in practice — type it
    dim::Int8
    has_points::Bool; has_lines::Bool; has_areas::Bool
    is_line_zero_len::Bool
    is_geom_empty::Bool
    # lazy caches
    unique_points::Union{Nothing, Set{...}}
    locator::Union{Nothing, RelatePointLocator{...}}
end
```

Port: `analyzeDimensions` (GI-trait traversal), `isZeroLength` (per-line all-segments-zero check using coordinate equality — exact), `getDimensionReal`, `hasDimension`, `hasAreaAndLine`, `hasEdges`, `hasBoundary`, `isPolygonal`, `isSelfNodingRequired` (trait-based: anything that can self-cross — port the Java condition), `getUniquePoints`, `getEffectivePoints`, `extract_segment_strings(rg, is_a, ext_filter)` building `RelateSegmentString`s (one per line / per ring, with `is_a`, `dim`, element id, ring id, parent polygonal ref, coordinate vector), `locate_*` delegations to the lazy locator.

`RelateSegmentString` port includes `prev_vertex`/`next_vertex`/ring wraparound and `create_node_section` — but with the symbolic twist: `create_node_section(ss, seg_index, node::NodeKey, ...)` takes a `NodeKey` instead of a constructed `Coordinate`. For `SS_TOUCH`/vertex incidences the key is `vertex_node(pt)`; for `SS_PROPER` it is the `crossing_node(...)` of the two segments.

TDD steps as usual. **Commit** — `Add RelateGeometry input facade and segment strings`

## Task 14: Generalize the JTS XML harness and vendor relate test files

**Files:**
- Modify: `test/external/jts/jts_testset_reader.jl`
- Create: `test/external/jts/relate_runner.jl`, `test/data/jts/general/` + `test/data/jts/validate/` (vendored XML), `test/data/jts/LICENSE-NOTICE.md`
- Test: `test/methods/relateng/xml_harness.jl` (parser smoke test; the full runner activates in Stage 5)

**Step 1: Vendor the XML files**

```bash
mkdir -p test/data/jts/general test/data/jts/validate
cp /Users/anshul/temp/GO_jts/jts/modules/tests/src/test/resources/testxml/general/TestRelate{PP,PL,PA,LL,LA,AA}.xml test/data/jts/general/
cp /Users/anshul/temp/GO_jts/jts/modules/tests/src/test/resources/testxml/general/TestBoundary.xml test/data/jts/general/
cp /Users/anshul/temp/GO_jts/jts/modules/tests/src/test/resources/testxml/validate/TestRelate*.xml test/data/jts/validate/
cp /Users/anshul/temp/GO_jts/jts/modules/tests/src/test/resources/testxml/misc/TestRelate{Empty,GC}.xml test/data/jts/general/ 2>/dev/null || true
cp /Users/anshul/temp/GO_jts/jts/modules/tests/src/test/resources/testxml/robust/TestRobustRelate*.xml test/data/jts/validate/ 2>/dev/null || true
```

(Adjust to actual filenames; the misc/robust names were reported without extensions verified — `ls` first.) Add `LICENSE-NOTICE.md` stating the files are from JTS (EPL/EDL dual license) with the upstream URL and commit.

**Step 2: Failing parser test** — `test/methods/relateng/xml_harness.jl`:

```julia
using Test
include(joinpath(@__DIR__, "..", "..", "external", "jts", "jts_testset_reader.jl"))

@testset "relate XML parsing" begin
    cases = load_test_cases(joinpath(@__DIR__, "..", "..", "data", "jts", "general", "TestRelatePP.xml"))
    @test !isempty(cases)
    item = first(first(cases).items)
    @test item.operation == "relate"
    @test item.expected_result isa Bool       # boolean ops parse as Bool now
    @test item.pattern isa String && length(item.pattern) == 9
end
```

**Step 3: Run** — FAIL (current reader parses every expected result as geometry, `test/external/jts/jts_testset_reader.jl:65`, and `TestItem` has no `pattern` field).

**Step 4: Generalize the reader** (`jts_testset_reader.jl`):

- Add `pattern::Union{Nothing, String}` to `TestItem`; parse by op kind:

```julia
const BOOLEAN_OPS = Set(["relate", "intersects", "disjoint", "contains", "within",
    "covers", "coveredby", "crosses", "touches", "overlaps", "equalstopo", "equals"])

function parse_expected(operation, raw::String)
    lowercase(operation) in BOOLEAN_OPS && return parse(Bool, lowercase(strip(raw)))
    return jts_wkt_to_geom(raw)
end
```

- In `parse_case`, read `arg3` from `op_attrs` when present (`get(op_attrs, "arg3", nothing)`) as the relate pattern; handle ops whose `arg1` isn't `"A"`/`"B"` by skipping with a counter.
- Delete the hardcoded `testfile` path and the executable overlay loop at lines 73–123 — move that overlay-running code to `test/external/jts/overlay_runner.jl` unchanged (it is not currently wired into CI, so this is a pure file move; do not silently delete it).
- Read `<precisionModel>` from the run header into a `Run` metadata field (used later to skip FIXED-precision cases).

**Step 5: Create `test/external/jts/relate_runner.jl`** — the case runner used in Stage 5, parameterized so it can be smoke-tested now:

```julia
"""
    run_relate_cases(relate_fn, pattern_fn, predicate_fns, files; skiplist)

For each XML case: `relate_fn(a, b)::DE9IM`, `pattern_fn(a, b, pattern)::Bool`,
`predicate_fns[opname](a, b)::Bool`. Cases whose (file, description, op) is in
`skiplist` are recorded as skipped, never silently dropped.
"""
```

Iterates cases, `@test`s expected values, collects a summary. Skiplist lives at `test/external/jts/relate_skiplist.jl` as a `Set{Tuple{String,String,String}}` with a mandatory comment per entry explaining the divergence.

**Step 6: Run parser test** — PASS. Register `xml_harness.jl` in the relateng test runner.

**Step 7: Commit** — `Generalize JTS XML test reader for relate operations and vendor test files`

---

# Stage 4 — Node topology

## Task 15: `NodeSection` + `NodeSections`

**Files:**
- Modify/Create: `src/methods/geom_relations/relateng/node_sections.jl` (struct exists since Task 11; add full API + `NodeSections`)
- Test: `test/methods/relateng/node_topology.jl`

**Java references:** `JTS:NodeSection.java` (201), `JTS:NodeSections.java` (122). No dedicated JUnit file — write unit tests for: `EdgeAngleComparator` ordering (via `rk_compare_edge_dir`), `isProper`/`isNodeAtVertex`, `NodeSections.prepareSections` ordering invariant, and `createNode` on a simple two-area touch (assert resulting edge count; full label assertions come with Task 17).

`NodeSection` final form (replacing the Task-11 minimal version):

```julia
struct NodeSection{P, G}
    is_a::Bool
    dim::Int8
    id::Int32
    ring_id::Int32
    polygonal::G                  # parent polygonal geometry or `nothing`
    is_node_at_vertex::Bool
    v0::Union{P, Nothing}         # vertex before node (nothing at line start)
    node::NodeKey{P}              # symbolic — JTS stores a Coordinate here
    v1::Union{P, Nothing}
end
```

`NodeSections` is `mutable struct` holding `node::NodeKey{P}` and `sections::Vector{NodeSection}`; port `addNodeSection`, `hasInteractionAB`, `getPolygonal(is_a)`, `createNode` (sort + `PolygonNodeConverter` delegation — stub the converter call until Task 16; in the meantime the test covers only single-polygon nodes), `prepareSections` (the comparator chain: lines before areas, grouped by polygon).

TDD steps as usual. **Commit** — `Add NodeSection and NodeSections with symbolic node keys`

## Task 16: `PolygonNodeConverter`

**Files:**
- Create: `src/methods/geom_relations/relateng/polygon_node_converter.jl`
- Test: `test/methods/relateng/node_topology.jl` (append)

**Java references:** `JTS:PolygonNodeConverter.java` (148); **`PolygonNodeConverterTest.java` (154 — port every test method first**; they construct section lists and assert converted output, which ports directly).

Port `convert`, `extractUnique`, `findShell`, `convertHoles`, `convertShellAndHoles`, `createSection` exactly; geometric comparisons go through the kernel comparator. Wire the real call into `NodeSections.createNode` (remove Task-15 stub) and re-run Task 15's tests.

TDD steps as usual. **Commit** — `Add PolygonNodeConverter for shell-hole node rewriting`

## Task 17: `RelateEdge` + `RelateNode`

**Files:**
- Create: `src/methods/geom_relations/relateng/relate_node.jl`
- Test: `test/methods/relateng/node_topology.jl` (append)

**Java references:** `JTS:RelateEdge.java` (362), `JTS:RelateNode.java` (230). No dedicated JUnit tests — write unit tests asserting full edge-wheel state for hand-built configurations (each verified by hand against the Java semantics before writing):

1. Two crossing lines at a vertex node → 4 edges, all dims L, no area labels.
2. Area corner (two edges of one polygon) → 2 edges with interior on the correct side (left of forward edge for CCW shell).
3. Area corner of A + line end of B at the same node → line edge labeled interior/boundary w.r.t. A correctly.
4. Two area corners (A and B) overlapping → collinear-edge merge case: coincident edges merged, area-over-line dim override.

`RelateEdge`:

```julia
mutable struct RelateEdge{P}
    node::NodeKey{P}
    dir_pt::P
    # per-geometry: dimension and left/right/on locations
    a_dim::Int8; a_loc_left::Int8; a_loc_right::Int8; a_loc_line::Int8
    b_dim::Int8; b_loc_left::Int8; b_loc_right::Int8; b_loc_line::Int8
end
```

Port the factory `relate_edge(node, dir_pt, is_a, dim, is_forward)`, `compare_to_edge` (via `rk_compare_edge_dir` with the symbolic node as apex), `merge!`, `set_area_interior!`, `is_known`, `is_interior`, `set_location!`, `get_location`, plus statics `find_known_edge_index` and `set_all_area_interior!`. `RelateNode` ports `add_edges!` (both arities), `add_line_edge`/`add_area_edge`/`add_edge` insertion-or-merge into the sorted wheel, and the `update_edges_in_area`/`update_if_area_prev/next` label propagation with circular `next_index`/`prev_index`.

**Additional step (added after Task 11 review):** Task 11's `AdjacentEdgeLocator` shipped with a private sequential slice of the node-wheel pipeline (`_AelEdge`, `_create_node_edges`, etc. in `point_locator.jl`, marked `TODO(Task 17)`). Once the real `NodeSections`/`RelateNode` exist, rewire `locate(::AdjacentEdgeLocator, p)` onto them (build `NodeSections`, push sections, `create_node`, `has_exterior_edge(node, true)`) and delete the slice. The ported AdjacentEdgeLocatorTest cases are the regression gate.

TDD steps as usual. **Commit** — `Add RelateNode and RelateEdge wheel with label propagation`

## Task 18: `TopologyComputer`

**Files:**
- Create: `src/methods/geom_relations/relateng/topology_computer.jl`
- Test: `test/methods/relateng/topology_computer.jl`

**Java reference:** `JTS:TopologyComputer.java` (520 — the heart of the topology layer; read fully, port in Java method order). No dedicated JUnit file; unit-test through the public entry points with a `RelateMatrixPredicate` attached and assert resulting IM strings:

```julia
# Example shape (values verified by hand/LibGEOS before writing):
@testset "addPointOnGeometry P vs A" begin
    # point in area interior: II entry must become 0... etc.
end
```

Cover: `initExteriorDims` for all dim pairs (P/P, P/L, P/A, L/L, L/A, A/A — assert the a-priori exterior entries match `TopologyComputer.java:44-102`), empty-geometry init, `add_point_on_point_interior!`/`_exterior!`, `add_point_on_geometry!`, `add_line_end_on_geometry!`, `add_area_vertex!`, `add_intersection!` + `evaluate_nodes!` (node-section grouping in `Dict{NodeKey{P}, NodeSections}`), `updateAreaAreaCross` via `rk_is_crossing`, short-circuit (`is_result_known` true → entry points become no-ops).

The struct:

```julia
mutable struct TopologyComputer{TP <: TopologyPredicate, RA, RB, P}
    predicate::TP
    geom_a::RA
    geom_b::RB
    node_sections::Dict{NodeKey{P}, NodeSectionsCollector}  # name per Task 15
end
```

Where the predicate's `require_self_noding(typeof(predicate))` is true **and** the geometry's `isSelfNodingRequired`, run the D3 coincidence-merge pass before `evaluate_nodes!`: group keys via `rk_nodes_coincide` (O(k²) over crossing keys — acceptable slow path; reference follow-up F1 in a comment).

TDD steps as usual. **Commit** — `Add TopologyComputer with symbolic node grouping`

---

# Stage 5 — The engine

## Task 19: Edge intersection (`EdgeSegmentIntersector` semantics over `SegSegClass`)

**Files:**
- Create: `src/methods/geom_relations/relateng/edge_intersector.jl`
- Test: `test/methods/relateng/edge_intersector.jl`

**Java references:** `JTS:EdgeSegmentIntersector.java` (89), the canonicality logic in `JTS:RelateSegmentString.java#isContainingSegment`.

Port `add_intersections!(computer, ssA, iA, ssB, iB; m, exact)`:

1. `rk_classify_intersection` on the two segments.
2. `SS_DISJOINT` → return.
3. `SS_PROPER` → one `crossing_node` key; `create_node_section` on both strings; `add_intersection!(computer, nsA, nsB)`.
4. `SS_TOUCH`/`SS_COLLINEAR` → for each incident vertex, create vertex-node sections. **Once-only rule:** JTS ensures a vertex shared by two adjacent segments of the same string is processed once via `isContainingSegment` (start-inclusive, end-exclusive containment). Our equivalent: attribute a vertex incidence to a segment only if the vertex is *not* that segment's end vertex (i.e. `a1_on_b` where the touch point equals `a1` is attributed to the *next* segment's `a0` — except for the final segment of an open line). Port the Java rule exactly; encode it as a predicate `_is_canonical_incidence(ss, seg_index, which_endpoint)` with its own unit tests for: mid-string vertex (two segments), string start, string end, ring wraparound.

Unit tests: hand-built two-string configurations asserting the exact set of `(NodeKey, section count)` recorded by a mock/`RelateMatrixPredicate`-backed computer — crossing, T-touch, shared endpoint, collinear overlap, adjacent-segment shared vertex (must produce sections once, not twice).

**Commit** — `Add edge segment intersector over symbolic classification`

## Task 20: Edge set enumeration via `SpatialTreeInterface`

**Files:**
- Modify: `edge_intersector.jl`
- Test: `test/methods/relateng/edge_intersector.jl` (append)

Replaces `JTS:EdgeSetIntersector.java` (HPRtree + monotone chains).

Implement `process_edge_intersections!(computer, ssa_list, ssb_list, accelerator; m, exact)`:

- Build per-segment extent lists for each side (`Extents.Extent` per segment).
- Accelerator selection mirrors the clipping pattern (`src/methods/clipping/clipping_processor.jl:8-62`): `NestedLoop` below a size threshold or on non-planar manifolds; STRtree-backed otherwise; `AutoAccelerator` picks by total segment count (use the same threshold constant as clipping; read it from the clipping source rather than inventing a new one).
- Tree path: `STRtree` over segment extents per side, then `dual_depth_first_search((i, j) -> ..., Extents.intersects, treeA, treeB)` (`src/utils/SpatialTreeInterface/dual_depth_first_search.jl:25`) mapping flat segment indices back to `(segment_string, seg_index)` pairs via an offset table.
- **Early exit:** after each pair, check `is_result_known(computer)`; since `dual_depth_first_search` has no built-in termination, throw a private `struct _RelateDone <: Exception end` from the callback and catch it at the call site. (Check first whether the traversal respects a `Break` return from `LoopStateMachine` — if so use that instead; note which in the commit.)

Tests: same fixtures as Task 19 run through both `NestedLoop` and the tree path must produce identical computer state; an early-exit test (intersects predicate over two heavily-overlapping rings) asserting the traversal stops (count pairs processed via a counting wrapper).

**Commit** — `Add accelerated edge set enumeration for RelateNG`

## Task 21: `RelateNG` algorithm type and evaluation phases

**Files:**
- Create: `src/methods/geom_relations/relateng/relate_ng.jl`
- Modify: `src/GeometryOps.jl`
- Test: `test/methods/relateng/relate_ng.jl`

**Java reference:** `JTS:RelateNG.java` (549) — port `evaluate` (lines 224–268) and its helpers `computePP` (305–331), `computeAtPoints` (333–361), `computeLineEnds` (392–457), `computeAreaVertex`-related (459–506), `computeAtEdges` (508–546) in order.

The algorithm type (mirrors `FosterHormannClipping` at `src/methods/clipping/clipping_processor.jl:51`):

```julia
struct RelateNG{M <: Manifold, A <: IntersectionAccelerator, E, BR <: BoundaryNodeRule} <: GeometryOpsCore.Algorithm{M}
    manifold::M
    accelerator::A
    exact::E
    boundary_rule::BR
end
RelateNG(; manifold::Manifold = Planar(), accelerator = AutoAccelerator(),
           exact = True(), boundary_rule = Mod2Boundary()) =
    RelateNG(manifold, accelerator, exact, boundary_rule)
RelateNG(m::Manifold; kw...) = RelateNG(; manifold = m, kw...)
GeometryOpsCore.manifold(alg::RelateNG) = alg.manifold
```

Entry points:

```julia
relate(alg::RelateNG, a, b) = ...          # RelateMatrixPredicate → DE9IM
relate(alg::RelateNG, a, b, pattern::String) = ...  # IMPatternMatcher → Bool
relate(a, b) = relate(RelateNG(), a, b)
relate(a, b, pattern::String) = relate(RelateNG(), a, b, pattern)
relate_predicate(alg::RelateNG, pred::TopologyPredicate, a, b)::Bool  # core
```

`relate_predicate` is the ported `evaluate`: build `RelateGeometry` for both sides, then phases 1–7 from the design doc (envelope screen using `require_interaction`/`require_covers` flags → `init_dims!` → exit-if-known → `init_bounds!` → exit-if-known → P/P fast path → points-vs-geometry (`computeAtPoints` both directions, line ends, area vertices, honoring `require_exterior_check`) → edge phase (extract segment strings filtered by interaction envelope, Task 20 enumeration, `evaluate_nodes!`) → `finish!` → `predicate_value`).

**Step 1 (RED): port `RelateNGTest.java` (686 lines).** This is the big one — port every test method into `test/methods/relateng/relate_ng.jl`, using the `RelateNGTestCase.java` helper shape:

```julia
function check_relate(awkt, bwkt, expected_im::String)
    a, b = from_wkt(awkt), from_wkt(bwkt)
    @test string(GO.relate(GO.RelateNG(), a, b)) == expected_im
end
function check_predicate(pred_factory, awkt, bwkt, expected::Bool) ... end
```

(`from_wkt` = the sanitizing WKT parser from the harness; move `jts_wkt_to_geom` into a shared test util include.) Port the whole file even though it takes a while — it is the primary correctness net for the engine and every case is a one-liner through these helpers. Skip prepared-mode checks until Task 22 (mark `@test_broken` or guard with a flag, then flip in Task 22).

**Step 2: Run** — FAIL (engine missing).

**Step 3: Implement** the engine; iterate per phase, running the test file after each helper lands. Debugging discipline: when a case disagrees, diff phase-by-phase against the Java (same WKT through JTS `RelateNG` semantics via LibGEOS `relate` if needed).

**Step 4: Run** — all ported tests PASS (GC-dependent cases may remain `@test_broken` if GC traversal gaps surface; record them in the task notes, fix before Task 23 declares the XML suite green).

**Step 5: Commit** — `Add RelateNG engine with phased evaluation` (commit earlier per-phase if the diff grows beyond ~500 lines).

## Task 22: Prepared mode

**Files:**
- Modify: `relate_ng.jl`
- Test: `test/methods/relateng/relate_ng.jl` (enable prepared checks)

```julia
struct PreparedRelate{RG, T, ALG <: RelateNG}
    alg::ALG
    geom_a::RG          # RelateGeometry with locator/unique-points caches forced
    edge_tree::T        # prebuilt segment tree for A (or nothing below threshold)
    segs_a::...         # extracted segment strings for A
end
prepare(alg::RelateNG, a) = ...
relate(p::PreparedRelate, b), relate(p::PreparedRelate, b, pattern), relate_predicate(p, pred, b)
```

Port the JTS caveat: predicates requiring self-noding bypass the cached interaction-envelope filtering (JTS `RelateNG.java` prepared branch). Tests: every `checkPrepared` case from `RelateNGTest.java` now enabled — prepared results must equal unprepared results across the whole ported suite (add a loop asserting that wholesale); plus a cache-reuse smoke test (two evaluations against one `PreparedRelate`).

**Commit** — `Add prepared mode for RelateNG`

## Task 23: Full JTS XML suite green

**Files:**
- Create: `test/external/jts/relate_skiplist.jl`
- Modify: `test/methods/relateng/xml_harness.jl` (activate the runner over all vendored files)
- Test: the XML suite itself

**Step 1:** Wire `run_relate_cases` (Task 14) to the real engine: `relate_fn = (a,b) -> GO.relate(GO.RelateNG(), a, b)`, `pattern_fn`, and predicate closures for every `BOOLEAN_OPS` name. Run over `test/data/jts/general/*.xml` first.

**Step 2:** Triage failures one file at a time (PP → PL → PA → LL → LA → AA → Boundary → Empty/GC), fixing engine bugs. Only after a documented analysis may a case go on the skiplist (legitimate reasons: FIXED precision-model cases; semantics divergence documented in the design doc). Every skiplist entry carries a comment with the case description and reason.

**Step 3:** Extend to `test/data/jts/validate/*.xml` (the big suites, ~11k lines of cases). Same triage discipline.

**Step 4:** Register the XML suite in `test/methods/relateng/runtests.jl` (guard the `validate/` set behind `if get(ENV, "GO_TEST_FULL_JTS", "true") == "true"` only if runtime proves prohibitive — measure first; target keeping it always-on).

**Step 5: Commit** — `Pass JTS relate XML test suites` (intermediate commits per fixed file are encouraged: `Fix <bug> found by TestRelateLL.xml`).

## Task 24: Cross-validation against existing GO predicates and LibGEOS

**Files:**
- Modify: `test/methods/geom_relations.jl`
- Test: same

**Step 1:** The existing suite compares `GO_f` vs `LG_f` per pair (`test/methods/geom_relations.jl:188-203`). Extend the comparison tuple: for every existing `@test_implementations` predicate check, add the RelateNG form, e.g.:

```julia
@test_implementations GO_f($g1, $g2) == GO.relateng_predicate(GO_f, $g1, $g2)
```

where a small mapping from the GO function to its `pred_*` factory lives in the test file. (Mechanically: add one line per existing comparison; keep the diff reviewable.)

**Step 2:** Run; divergences are triaged — a divergence where RelateNG agrees with LibGEOS and old-GO disagrees is a *known-gap win*: record it in a `KNOWN_OLD_GO_GAPS` list asserting the new behavior, don't regress the new engine to match the old.

**Step 3: Commit** — `Cross-validate RelateNG against existing predicates and LibGEOS`

## Task 25: Differential fuzzing vs LibGEOS

**Files:**
- Create: `test/methods/relateng/fuzz.jl`
- Test registration: relateng runner (seeded, bounded case count for CI)

**Step 1:** Generator: random polygon/line/point pairs reusing `test/data/polygon_generation.jl` and adversarial constructors — shared vertices, collinear edges, ulp-perturbed near-crossings (`nextfloat`/`prevfloat` on crossing configurations), zero-length lines, rings touching at points.

**Step 2:** Property: `string(GO.relate(GO.RelateNG(), a, b)) == LibGEOS.relate(lg(a), lg(b))` for N seeded cases (N≈500 in CI; overridable via ENV for deep runs). On divergence, print the WKT pair; triage: GEOS-version check (`LibGEOS.GEOS_VERSION >= v"3.13"` asserted up front so the oracle is RelateNG-native), exactness wins (we are right where FP GEOS is wrong — verify by hand with rational arithmetic, then add to a documented `EXACTNESS_WINS` fixture list asserting *our* answer), real bugs (fix).

**Step 3: Commit** — `Add differential fuzzing of RelateNG against LibGEOS`

---

# Stage 6 — Surface, docs, performance

## Task 26: Public API surface and named-predicate methods

**Files:**
- Modify: `relate_ng.jl`, `src/GeometryOps.jl` (exports)
- Test: `test/methods/relateng/relate_ng.jl` (append)

**Step 1 (RED):** API tests: `GO.relate(a, b) isa GO.DE9IM`; `GO.relate(a, b, "T*F**FFF*") isa Bool`; `GO.intersects(GO.RelateNG(), a, b) == GO.intersects(a, b)` for each of the 10 named predicates (opt-in form, design D4 — existing defaults untouched); `prepare`d equivalents.

**Step 2:** Implement the named-predicate methods:

```julia
intersects(alg::RelateNG, g1, g2) = relate_predicate(alg, pred_intersects(), g1, g2)
# ... one line each for disjoint, contains, within, covers, coveredby,
#     crosses, overlaps, touches, equals (equals → pred_equalstopo)
```

Export `relate`, `DE9IM`, `RelateNG`, `prepare`, `Mod2Boundary` (+ other rules). Check for export collisions (`prepare` vs `prepare_naturally` is fine; verify `relate` is genuinely new — Task 0's survey said yes).

**Step 3:** Run full relateng test battery + the whole GO test suite (`julia --project=docs -e 'using Pkg; Pkg.test("GeometryOps")'` or the established command) to catch collisions.

**Step 4: Commit** — `Export \`relate\` and RelateNG predicate methods`

## Task 27: Literate docs

**Files:**
- Modify: every `src/methods/geom_relations/relateng/*.jl` (top-of-file literate headers where still missing), `relate_ng.jl` (main docstring + `@example` blocks)

Follow the house literate convention (`src/methods/geom_relations/within.jl:1-53` as the template): `# # Relate` header, `export` line, `#= ... =#` block with a "What is relate?/DE-9IM" explainer and a small Makie `@example` (two overlapping polygons, show the matrix string), then `## Implementation`. Docstrings for `relate`, `RelateNG`, `DE9IM`, `prepare`, boundary rules. Check whether docs build includes new files automatically (look at `docs/make.jl` page generation) and add pages if manual.

Run docs build: `julia --project=docs docs/make.jl` — expect no errors (or at minimum no new errors; note preexisting failures).

**Commit** — `Add literate documentation for RelateNG`

## Task 28: Benchmarks and allocation discipline

**Files:**
- Create: `benchmarks/relateng.jl`
- Create: `test/methods/relateng/allocations.jl`

**Step 1: Benchmarks** (Chairmarks `@b`, matching `benchmarks/` house style): RelateNG vs old GO predicates vs LibGEOS (plain + prepared GEOS) across polygon sizes from the existing providers/`polygon_generation.jl`; two workload shapes: `intersects` (early-exit-friendly) and full `relate` (no exit). Print a comparison table; no CI gating.

**Step 2: Allocation tests:** after warmup, `@allocated GO.relate_predicate(alg, GO.pred_intersects(), a, b)` for mid-size polygon pairs — assert below a budget that excludes the unavoidable per-call setup (RelateGeometry, dicts) but catches per-segment-pair allocation regressions: establish the measured baseline, assert `<= 2x baseline` with the measured number recorded in a comment. Type-stability spot checks with `@inferred` on `rk_classify_intersection`, `update_dim!`, `relate_predicate`.

**Step 3:** Profile one large A/A `relate` case; if >20% of time is in the kernel boundary or dict ops, file observations in the follow-up register (F1 territory) — do not optimize speculatively in this task.

**Step 4: Commit** — `Add RelateNG benchmarks and allocation tests`

---

## Done criteria (whole plan)

1. All vendored JTS XML relate suites pass (skiplist documented and short).
2. Ported `RelateNGTest` + component JUnit suites pass.
3. Existing 63-pair predicate suite agrees (modulo documented old-GO gaps).
4. Seeded fuzz vs LibGEOS ≥ 3.13 clean (modulo documented exactness wins).
5. Kernel conformance suite green (the Spherical contract).
6. No FP-constructed intersection coordinate anywhere in `src/.../relateng/` (grep for `getIntersection`-style construction; `_exact_crossing_point` rationals are the only intersection-point computation and only inside `rk_nodes_coincide`).
7. Benchmarks recorded; hot path allocation-bounded.

Out of scope (follow-up register in the design doc): F1 fast coincidence filter, F2 Spherical kernel, F3 default flip + old-processor deprecation, F4 generalized prepared geometry, F5 extra GC fuzzing.
