# Prepared Geometry (F4) Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement the `Prepared{Geom, Preparations}` wrapper from
`docs/plans/2026-07-01-prepared-geometry-design.md`: capability-trait-indexed, eagerly
built, immutable preparations attached to geometries via recursive wrapping, with RelateNG
consuming them.

**Architecture:** Consumption interface (types, traits, `Prepared`, `get`/`getprep`,
GeoInterface forwarding) lives in **GeometryOpsCore** (`src/types/preparations.jl`). The
`prepare` builder and concrete v1 specs (`RingEdgeIndex`, `ChildTree`, `PointInArea`) live
in **GeometryOps** (`src/prepared/prepared.jl`) because extents are manifold-aware (they
need `rk_interaction_bounds` from the relateng kernel). RelateNG gains two seams: skip its
private extent-cache rebuild for matched-manifold `Prepared` inputs, and reuse a
`PointInAreaLike` preparation in its point locator.

**Tech Stack:** Julia, GeoInterface (`GI`), Extents, the in-repo `NaturalIndexing` and
`SpatialTreeInterface` modules, relateng kernel helpers (`rk_interaction_bounds`,
`_segment_extent`, `_to_kernel_point`), `IndexedPointInAreaLocator`.

---

## Deviations from the design doc (discovered during API fact-finding; Task 12 records them in the design doc)

1. **`SegmentIndex` (whole-geometry flattened edge tree) is dropped from v1.** Its
   `owners::Vector{NTuple{2,Int}}` table only makes sense relative to relateng's
   `extract_segment_strings` traversal order; building it independently risks silent
   owner-index mismatch. RelateNG keeps building its own `PreparedEdgeIndex` inside
   `prepare(::RelateNG, a)`. Revisit when a second consumer needs a flat segment tree.
2. **The clipping/`AutoAccelerator` seam is deferred.** Clipping accelerators index
   flattened edge tables (`NaturalIndex(edges)`), not child extents — same owners problem.
3. **The `tuples` keyword is deferred.** Callers compose: `prepare(GO.tuples(geom); ...)`.
4. **`build` is named `buildprep`** (a bare exported `build` invites collisions), and
   `appliesto`/`buildprep` are deliberately **not exported** — providers extend them via
   `import GeometryOpsCore: appliesto, buildprep`.
5. **Topmost-wins spec attachment:** a spec is consumed at the highest tree level whose
   trait matches `appliesto(spec)` and is not offered further down. E.g. `ChildTree()` on a
   MultiPolygon indexes polygon extents at the MP node; it does not additionally build
   ring-trees inside each polygon.

## Conventions for the executor

- **Worktree:** use superpowers:using-git-worktrees to create a worktree on a new branch
  `prepared-geometry` based on `relateng-spherical`. All paths below are relative to the
  worktree root.
- **Run single test files from the repo root** (never `Pkg.test()` in the loop; the full
  suite takes ~25 min): `julia --project=test -e 'include("test/methods/relateng/relate_ng.jl")'`.
  First run in a fresh depot precompiles (~1–2 min); later runs are seconds.
- `timeout` does not exist on this macOS — use `gtimeout` if you need one.
- **Commit style** (AGENTS.md): imperative, capitalized, no conventional-commit prefixes,
  backticks around identifiers, no trailing period.
- TDD: write the failing test, SEE it fail, implement, SEE it pass, commit.
- `GI` is `GeoInterface`; inside GeometryOpsCore it is already imported; test files are
  self-contained (they do their own `using`).

---

### Task 0: Worktree, baseline, commit the plan docs

**Step 1: Create worktree + branch**

Follow superpowers:using-git-worktrees. Base: `relateng-spherical`. Branch name:
`prepared-geometry`.

**Step 2: Verify baseline is green**

Run: `julia --project=test -e 'include("test/methods/relateng/relate_ng.jl")'`
Expected: all tests pass (this file includes the prepared-vs-unprepared fixture sweep).

**Step 3: Commit the design + plan docs**

The two docs were written in the main checkout; copy them into the worktree if the
worktree was created before they existed.

```bash
git add docs/plans/2026-07-01-prepared-geometry-design.md docs/plans/2026-07-01-prepared-geometry.md
git commit -m "Add prepared-geometry design and implementation plan"
```

---

### Task 1: Core preparation types + trait-indexed retrieval

**Files:**
- Create: `GeometryOpsCore/src/types/preparations.jl`
- Modify: `GeometryOpsCore/src/GeometryOpsCore.jl` (include line)
- Create: `test/core/preparations.jl`
- Modify: `test/runtests.jl` (wire the new file into the `"Core"` testset, mirroring the
  existing `test/core/*.jl` entries at lines ~8–13)

**Step 1: Write the failing test**

Create `test/core/preparations.jl`:

```julia
using Test
import GeometryOpsCore
import GeometryOpsCore: Prepared, AbstractPreparation, AbstractPreparationTrait,
    SpatialIndexLike, SpatialEdgeIndexLike, PointInAreaLike,
    preptrait, getprep, Planar, Spherical, manifold
import GeoInterface as GI
import GeoInterface: Extents

# Dummy preparations for interface tests
struct _DummyEdgeTree <: AbstractPreparation
    payload::Int
end
GeometryOpsCore.preptrait(::_DummyEdgeTree) = SpatialEdgeIndexLike()

struct _DummyChildTree <: AbstractPreparation end
GeometryOpsCore.preptrait(::_DummyChildTree) = SpatialIndexLike()

@testset "Preparation retrieval" begin
    ring = GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)])
    ext = GI.extent(ring; fallback = true)
    p = Prepared(ring, Planar(), (_DummyEdgeTree(42), _DummyChildTree()), ext)

    @test manifold(p) === Planar()
    @test parent(p) === ring

    # hit: first matching preparation, concretely typed
    et = get(p, SpatialEdgeIndexLike())
    @test et isa _DummyEdgeTree && et.payload == 42
    @test get(p, SpatialIndexLike()) isa _DummyChildTree

    # miss: nothing
    @test get(p, PointInAreaLike()) === nothing

    # plain geometries always miss — uniform call sites
    @test get(ring, SpatialEdgeIndexLike()) === nothing

    # manifold-checked retrieval: mismatch is a miss, not an error
    @test getprep(Planar(), p, SpatialEdgeIndexLike()) === et
    @test getprep(Spherical(), p, SpatialEdgeIndexLike()) === nothing
    @test getprep(Planar(), ring, SpatialEdgeIndexLike()) === nothing
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=test -e 'include("test/core/preparations.jl")'`
Expected: FAIL — `UndefVarError` / `ArgumentError` importing `Prepared` etc.

**Step 3: Write the implementation**

Create `GeometryOpsCore/src/types/preparations.jl`:

```julia
#=
# Prepared geometry

Interface types for prepared geometry: a GeoInterface-compatible wrapper that carries a
tuple of immutable, eagerly built acceleration structures ("preparations"), retrieved by
*capability* trait — algorithms only care that a geometry *has* an edge index, not how it
was built. Design: `docs/plans/2026-07-01-prepared-geometry-design.md` (repo root).

The consumption interface lives here in GeometryOpsCore so third-party packages can
provide preparations without depending on GeometryOps. The `prepare` builder and the
concrete preparations live in GeometryOps (extents are manifold-aware).
=#

export AbstractPreparation, AbstractPreparationTrait,
    SpatialIndexLike, SpatialEdgeIndexLike, PointInAreaLike,
    preptrait, PreparationSpec, Prepared, getprep, prepare

"""
    AbstractPreparation

Supertype of built, immutable acceleration structures attached to exactly one geometry via
[`Prepared`](@ref). Concrete subtypes implement [`preptrait`](@ref).
"""
abstract type AbstractPreparation end

"""
    AbstractPreparationTrait

Capability traits for preparations: what a preparation can answer, not what it is.
Retrieval via `get(prepared, trait)` returns the first preparation with that capability.
"""
abstract type AbstractPreparationTrait end

"Capability: an extent tree over child elements (rings of a polygon, polygons of a multi)."
struct SpatialIndexLike <: AbstractPreparationTrait end
"Capability: an extent tree over the segments/edges of a curve or geometry."
struct SpatialEdgeIndexLike <: AbstractPreparationTrait end
"Capability: a point-in-polygonal-area locator."
struct PointInAreaLike <: AbstractPreparationTrait end

"""
    preptrait(prep::AbstractPreparation)::AbstractPreparationTrait

The capability of `prep`. Every concrete preparation implements this.
"""
function preptrait end

"""
    PreparationSpec

Supertype of the *configs* passed to [`prepare`](@ref)'s `preps` tuple. A spec declares
which GeoInterface trait it attaches to via `appliesto(spec)` and how to build via
`buildprep(spec, manifold, geom) -> Union{AbstractPreparation, Nothing}` (`nothing` means
"nothing to index here", e.g. an empty geometry). Both are extended (not exported) via
`import GeometryOpsCore: appliesto, buildprep`.
"""
abstract type PreparationSpec end

function appliesto end
function buildprep end

"""
    Prepared(geom, manifold::Manifold, preps::Tuple, extent::Extents.Extent)

A geometry `geom` carrying preparations `preps`. Forwards GeoInterface, so it works
anywhere a geometry does. `extent` is always cached (manifold-aware: 3D on `Spherical`).
The wrapper and its preparation set are immutable; a miss on `get` means the caller takes
its unindexed path. Children of a `Prepared` built by `prepare` are themselves `Prepared`.
"""
struct Prepared{T <: GI.AbstractTrait, M <: Manifold, G,
                P <: Tuple{Vararg{AbstractPreparation}}, E <: Extents.Extent}
    geom::G
    manifold::M
    preps::P
    extent::E
    function Prepared(geom::G, m::M, preps::P, extent::E) where
            {M <: Manifold, G, P <: Tuple{Vararg{AbstractPreparation}}, E <: Extents.Extent}
        t = GI.trait(geom)
        t === nothing && throw(ArgumentError("`Prepared` requires a GeoInterface geometry; got $(typeof(geom))"))
        new{typeof(t), M, G, P, E}(geom, m, preps, extent)
    end
end

Base.parent(p::Prepared) = p.geom
manifold(p::Prepared) = p.manifold

"""
    get(p::Prepared, t::AbstractPreparationTrait) -> Union{AbstractPreparation, Nothing}
    get(geom, t::AbstractPreparationTrait) -> nothing

The first preparation on `p` with capability `t`, or `nothing`. Plain geometries always
return `nothing`, so call sites need no special-casing.
"""
Base.get(p::Prepared, t::AbstractPreparationTrait) = _findprep(t, p.preps)
Base.get(@nospecialize(_), ::AbstractPreparationTrait) = nothing

_findprep(::AbstractPreparationTrait, ::Tuple{}) = nothing
_findprep(t::AbstractPreparationTrait, preps::Tuple) =
    preptrait(first(preps)) === t ? first(preps) : _findprep(t, Base.tail(preps))

"""
    getprep(m::Manifold, geom, t::AbstractPreparationTrait) -> Union{AbstractPreparation, Nothing}

Manifold-checked retrieval: like `get(geom, t)` but a manifold mismatch is a miss —
preparations bake in manifold-specific state (kernel point types, 3D extents), so a
`Planar` edge tree must never be served to a `Spherical` algorithm.
"""
getprep(m::Manifold, p::Prepared, t::AbstractPreparationTrait) =
    manifold(p) === m ? get(p, t) : nothing
getprep(::Manifold, @nospecialize(_), ::AbstractPreparationTrait) = nothing

"""
    prepare(...)

Generic entry point for prepared-geometry optimizations. GeometryOps implements
`prepare(geom; preps, manifold)` (the wrapper builder) and `prepare(alg, geom)` methods
for algorithms (e.g. `RelateNG`).
"""
function prepare end
```

In `GeometryOpsCore/src/GeometryOpsCore.jl`, add the include directly after
`include("types/manifold.jl")`:

```julia
include("types/manifold.jl")
include("types/preparations.jl")
```

Wire the test file into `test/runtests.jl` next to the other `test/core/*.jl` entries,
matching their exact style (look at how `core/manifold.jl` is wired and copy it), e.g.:

```julia
@safetestset "Core preparations" begin include("core/preparations.jl") end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=test -e 'include("test/core/preparations.jl")'`
Expected: PASS.

**Step 5: Commit**

```bash
git add GeometryOpsCore/src/types/preparations.jl GeometryOpsCore/src/GeometryOpsCore.jl test/core/preparations.jl test/runtests.jl
git commit -m "Add preparation types and \`Prepared\` wrapper to GeometryOpsCore"
```

---

### Task 2: GeoInterface forwarding for `Prepared`

**Files:**
- Modify: `GeometryOpsCore/src/types/preparations.jl` (append)
- Modify: `test/core/preparations.jl` (append)

**Step 1: Write the failing test**

Append to `test/core/preparations.jl`:

```julia
@testset "GeoInterface forwarding" begin
    ring = GI.LinearRing([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 0.0)])
    rext = GI.extent(ring; fallback = true)
    pring = Prepared(ring, Planar(), (_DummyEdgeTree(1),), rext)

    @test GI.isgeometry(pring)
    @test GI.trait(pring) isa GI.LinearRingTrait
    @test GI.npoint(pring) == GI.npoint(ring)
    @test collect(GI.getpoint(pring)) == collect(GI.getpoint(ring))
    @test GI.coordinates(pring) == GI.coordinates(ring)
    @test GI.is3d(pring) == false
    @test GI.extent(pring) == rext          # served from the field
    @test Extents.extent(pring) == rext

    # a Prepared polygon whose ring is itself Prepared: GI hands back the prepared child
    poly = GI.Polygon([pring])
    pext = GI.extent(ring; fallback = true)
    ppoly = Prepared(poly, Planar(), (), pext)
    @test GI.trait(ppoly) isa GI.PolygonTrait
    @test GI.nring(ppoly) == 1
    @test GI.getring(ppoly, 1) === pring
    @test GI.getexterior(ppoly) === pring
    @test get(GI.getring(ppoly, 1), SpatialEdgeIndexLike()) isa _DummyEdgeTree
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=test -e 'include("test/core/preparations.jl")'`
Expected: FAIL — `MethodError` on `GI.npoint`/`GI.trait` etc. for `Prepared`.

**Step 3: Write the implementation**

Append to `GeometryOpsCore/src/types/preparations.jl` (adapted from the fixed forwarding
on the `as/prepare` branch, `git show origin/as/prepare:src/preparations/prepared_geometry.jl`):

```julia
#-- GeoInterface forwarding: a Prepared IS-A geometry wherever GI is spoken.
#-- The extent is served from the cached field, never recomputed.
GI.isgeometry(::Type{<:Prepared{T, M, G}}) where {T, M, G} = GI.isgeometry(G)
GI.trait(::Prepared{T}) where {T} = T()
GI.geomtrait(::Prepared{T}) where {T} = T()

Extents.extent(p::Prepared) = p.extent
GI.extent(::GI.AbstractTrait, p::Prepared) = p.extent
GI.crs(t::GI.AbstractTrait, p::Prepared) = GI.crs(t, parent(p))

for f in (:coordnames, :is3d, :ismeasured, :isempty, :coordinates, :ngeom, :getgeom)
    @eval GI.$f(t::GI.AbstractGeometryTrait, p::Prepared, args...) = GI.$f(t, parent(p), args...)
end
for f in (:npoint, :getpoint, :startpoint, :endpoint, :issimple, :isclosed, :isring)
    @eval GI.$f(t::GI.AbstractCurveTrait, p::Prepared, args...) = GI.$f(t, parent(p), args...)
end
for f in (:nring, :getring, :getexterior, :nhole, :gethole, :npoint, :getpoint)
    @eval GI.$f(t::GI.AbstractPolygonTrait, p::Prepared, args...) = GI.$f(t, parent(p), args...)
end
for f in (:nring, :getring, :npoint, :getpoint)
    @eval GI.$f(t::GI.AbstractMultiPolygonTrait, p::Prepared, args...) = GI.$f(t, parent(p), args...)
end
for f in (:npoint, :getpoint, :issimple)
    @eval GI.$f(t::GI.AbstractMultiPointTrait, p::Prepared, args...) = GI.$f(t, parent(p), args...)
    @eval GI.$f(t::GI.AbstractMultiCurveTrait, p::Prepared, args...) = GI.$f(t, parent(p), args...)
end

#-- disambiguation against GI's trait-specific defaults (e.g. npoint(::LineTrait) = 2)
for T in (:LineTrait, :TriangleTrait, :PentagonTrait, :HexagonTrait, :RectangleTrait, :QuadTrait)
    @eval GI.npoint(t::GI.$T, p::Prepared) = GI.npoint(t, parent(p))
end
for T in (:TriangleTrait, :RectangleTrait, :QuadTrait, :PentagonTrait, :HexagonTrait)
    @eval GI.nring(t::GI.$T, p::Prepared) = GI.nring(t, parent(p))
end
```

If the test surfaces a genuine method ambiguity the lists above don't cover, resolve it
the same way (a specific `(concrete trait, Prepared)` method that forwards) — do not
broaden signatures.

**Step 4: Run test to verify it passes**

Run: `julia --project=test -e 'include("test/core/preparations.jl")'`
Expected: PASS (both testsets).

**Step 5: Commit**

```bash
git add GeometryOpsCore/src/types/preparations.jl test/core/preparations.jl
git commit -m "Forward GeoInterface through \`Prepared\`"
```

---

### Task 3: Surface the Core names in GeometryOps; `RelateNG.prepare` extends the Core generic

**Files:**
- Modify: `src/GeometryOps.jl` (the `import GeometryOpsCore:` block at lines ~5–19 and the
  `export` line at line ~20)
- Modify: `src/methods/geom_relations/relateng/relate_ng.jl:3` (export line comment only —
  see below)

**Step 1: Write the failing test**

Create `test/prepared/prepared_geometry.jl`:

```julia
using Test
import GeometryOps as GO
import GeometryOps: Prepared, prepare, getprep,
    SpatialIndexLike, SpatialEdgeIndexLike, PointInAreaLike
import GeometryOpsCore
import GeoInterface as GI
import GeoInterface: Extents

@testset "prepare is one generic function" begin
    # the algorithm method and the (future) geometry method share the Core binding
    @test GO.prepare === GeometryOpsCore.prepare
    poly = GI.Polygon([GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)])])
    prep = GO.prepare(GO.RelateNG(), poly)
    @test prep isa GO.PreparedRelate
end
```

Wire into `test/runtests.jl` after the RelateNG entry (line ~34):

```julia
@safetestset "Prepared geometry" begin include("prepared/prepared_geometry.jl") end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: FAIL — `GO.prepare === GeometryOpsCore.prepare` is `false` (relateng currently
*defines* its own `prepare` function rather than extending Core's), or an import error for
`Prepared` from GeometryOps.

**Step 3: Write the implementation**

In `src/GeometryOps.jl`:

1. Add to the `import GeometryOpsCore:` list (after the `manifold, best_manifold,` line):

```julia
                AbstractPreparation, AbstractPreparationTrait,
                SpatialIndexLike, SpatialEdgeIndexLike, PointInAreaLike,
                preptrait, PreparationSpec, Prepared, getprep, prepare,
                appliesto, buildprep,
```

2. Extend the `export` line (line ~20) with the user-facing names:

```julia
export TraitTarget, Manifold, Planar, Spherical, Geodesic, apply, applyreduce, flatten, reconstruct, rebuild, unwrap, get_geometries
export Prepared, prepare, getprep, preptrait, SpatialIndexLike, SpatialEdgeIndexLike, PointInAreaLike
```

Because `prepare` is now *imported*, the existing `function prepare(alg::RelateNG, a)` at
`relate_ng.jl:665` automatically becomes a method of `GeometryOpsCore.prepare` — no code
change needed there. The `export relate, RelateNG, prepare` at `relate_ng.jl:3` stays (a
duplicate export of the same binding is harmless).

**Step 4: Run tests to verify they pass**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: PASS.

Then confirm nothing in relateng broke (this exercises the whole prepared fixture sweep):

Run: `julia --project=test -e 'include("test/methods/relateng/relate_ng.jl")'`
Expected: PASS.

**Step 5: Commit**

```bash
git add src/GeometryOps.jl test/prepared/prepared_geometry.jl test/runtests.jl
git commit -m "Re-export prepared-geometry interface from GeometryOps and unify \`prepare\`"
```

---

### Task 4: The `prepare` builder (recursive wrapping, no specs yet)

**Files:**
- Create: `src/prepared/prepared.jl`
- Modify: `src/GeometryOps.jl` (add `include("prepared/prepared.jl")` on the line directly
  after `include("methods/geom_relations/relateng/relate_ng.jl")` — it uses relateng
  kernel helpers)
- Modify: `test/prepared/prepared_geometry.jl` (append)

**Step 0: Verify helper names**

The builder uses relateng-internal helpers. Confirm they exist with these exact names
before writing code; if a name differs, use the actual one:

```bash
grep -n "_to_kernel_point\|_kernel_point_type" src/methods/geom_relations/relateng/kernel*.jl | head
grep -n "function rk_interaction_bounds\|rk_interaction_bounds(::\|rk_interaction_bounds(m" src/methods/geom_relations/relateng/kernel*.jl | head
grep -n "_segment_extent_type\|_segment_extent(" src/methods/geom_relations/relateng/edge_intersector.jl
```

Expected: `_to_kernel_point(::Planar, p)` / `(::Spherical, p)` (or equivalent),
`rk_interaction_bounds(m, geom)`, `_segment_extent(m, p, q)`, `_segment_extent_type(m)`.

**Step 1: Write the failing test**

Append to `test/prepared/prepared_geometry.jl`:

```julia
# shared fixtures for the remaining testsets
const _PG_RING_OUTER = GI.LinearRing([(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)])
const _PG_RING_HOLE  = GI.LinearRing([(4.0, 4.0), (6.0, 4.0), (6.0, 6.0), (4.0, 6.0), (4.0, 4.0)])
const _PG_POLY  = GI.Polygon([_PG_RING_OUTER, _PG_RING_HOLE])
const _PG_POLY2 = GI.Polygon([GI.LinearRing([(20.0, 0.0), (25.0, 0.0), (25.0, 5.0), (20.0, 0.0)])])
const _PG_MP    = GI.MultiPolygon([_PG_POLY, _PG_POLY2])
# mid-latitude spherical quad (lon/lat degrees)
const _PG_SPH_POLY = GI.Polygon([GI.LinearRing([(10.0, 40.0), (20.0, 40.0), (20.0, 50.0), (10.0, 50.0), (10.0, 40.0)])])

@testset "prepare: recursive wrapping, no specs" begin
    p = prepare(_PG_POLY)
    @test p isa Prepared
    @test GI.trait(p) isa GI.PolygonTrait
    @test GeometryOpsCore.manifold(p) === GO.Planar()
    @test p.preps === ()

    # children are themselves Prepared, with per-node extents
    r1 = GI.getring(p, 1)
    @test r1 isa Prepared && GI.trait(r1) isa GI.LinearRingTrait
    @test GI.extent(r1) == GI.extent(_PG_RING_OUTER; fallback = true)
    @test GI.extent(p) == GI.extent(_PG_POLY; fallback = true)

    # geometry content is unchanged
    @test GI.coordinates(p) == GI.coordinates(_PG_POLY)

    # points pass through unwrapped
    pt = GI.Point(1.0, 2.0)
    @test prepare(pt) === pt

    # multipolygon: every level wrapped
    mp = prepare(_PG_MP)
    @test mp isa Prepared && GI.getgeom(mp, 1) isa Prepared
    @test GI.getring(GI.getgeom(mp, 1), 1) isa Prepared

    # spherical: extents are 3D unit-sphere boxes
    ps = prepare(_PG_SPH_POLY; manifold = GO.Spherical())
    @test GI.extent(ps) isa Extents.Extent{(:X, :Y, :Z)}
    @test GI.extent(GI.getring(ps, 1)) isa Extents.Extent{(:X, :Y, :Z)}
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: FAIL — `MethodError: no method matching prepare(::GI.Polygon...)`.

**Step 3: Write the implementation**

Create `src/prepared/prepared.jl`:

```julia
#=
# Prepared geometry: builder and v1 preparations

`prepare(geom; preps, manifold)` rebuilds the GeoInterface tree once, bottom-up, wrapping
each node in `GeometryOpsCore.Prepared` with a manifold-aware cached extent
(`rk_interaction_bounds`), and attaching each spec's built preparation at the *highest*
tree level whose trait matches `appliesto(spec)` (topmost-wins). Modeled on relateng's
`_relate_cache_extents` (`relate_geometry.jl`), which it subsumes for `Prepared` inputs.

Design: `docs/plans/2026-07-01-prepared-geometry-design.md`.
=#

export RingEdgeIndex, ChildTree, PointInArea

"""
    prepare(geom; preps::Tuple = (), manifold::Manifold = Planar())

Wrap `geom` (and, recursively, its parts) in [`Prepared`](@ref), eagerly building the
preparations requested in `preps` (a tuple of `PreparationSpec`s such as
[`RingEdgeIndex`](@ref), [`ChildTree`](@ref), [`PointInArea`](@ref)).

Every wrapped node caches its manifold-aware extent. A spec is consumed at the highest
matching level of the tree. Preparations do not survive coordinate transformations —
re-`prepare` after `apply`-style rebuilds.
"""
prepare(geom; preps::Tuple = (), manifold::Manifold = Planar()) =
    _prepare(manifold, preps, GI.trait(geom), geom)

#-- split the spec tuple: consumed at this node vs offered to children (topmost-wins)
_specs_here(specs::Tuple, trait) = filter(s -> trait isa appliesto(s), specs)
_specs_below(specs::Tuple, trait) = filter(s -> !(trait isa appliesto(s)), specs)

#-- points and multipoints pass through: their extent is themselves
_prepare(::Manifold, specs, ::Union{GI.AbstractPointTrait, GI.AbstractMultiPointTrait}, geom) = geom

#-- curve leaves: the only level where coordinates are read
function _prepare(m::Manifold, specs, trait::GI.AbstractCurveTrait, line)
    GI.isempty(line) && return line
    return _wrap(m, _specs_here(specs, trait), line, rk_interaction_bounds(m, line))
end

function _prepare(m::Manifold, specs, trait::GI.AbstractPolygonTrait, poly)
    GI.isempty(poly) && return poly
    below = _specs_below(specs, trait)
    rings = [_prepare(m, below, GI.trait(r), r) for r in GI.getring(poly)]
    rebuilt = GI.geointerface_geomtype(trait)(rings; crs = GI.crs(poly))
    ext = _union_prepared_extents(m, rings)
    ext === nothing && return rebuilt
    return _wrap(m, _specs_here(specs, trait), rebuilt, ext)
end

#-- collections (covers Multi* types too)
function _prepare(m::Manifold, specs, trait::GI.AbstractGeometryCollectionTrait, geom)
    GI.isempty(geom) && return geom
    below = _specs_below(specs, trait)
    children = [_prepare(m, below, GI.trait(g), g) for g in GI.getgeom(geom)]
    rebuilt = GI.geointerface_geomtype(trait)(children; crs = GI.crs(geom))
    ext = _union_prepared_extents(m, children)
    ext === nothing && return rebuilt
    return _wrap(m, _specs_here(specs, trait), rebuilt, ext)
end

#-- unknown traits pass through untouched
_prepare(::Manifold, specs, ::GI.AbstractTrait, geom) = geom

function _wrap(m::Manifold, specs::Tuple, geom, ext)
    return Prepared(geom, m, _build_preps(specs, m, geom), ext)
end

#-- build each consumed spec; a spec may decline (empty geometry) by returning `nothing`
_build_preps(::Tuple{}, m, geom) = ()
function _build_preps(specs::Tuple, m, geom)
    p = buildprep(first(specs), m, geom)
    rest = _build_preps(Base.tail(specs), m, geom)
    return p === nothing ? rest : (p, rest...)
end

function _union_prepared_extents(m::Manifold, children)
    ext = nothing
    for c in children
        ce = if c isa Prepared
            c.extent
        elseif GI.isempty(c)
            nothing
        else
            rk_interaction_bounds(m, c)
        end
        ce === nothing && continue
        ext = ext === nothing ? ce : Extents.union(ext, ce)
    end
    return ext
end
```

Note: if `rk_interaction_bounds` has no method for bare point children inside collections,
mirror whatever `_relate_extent` (`relate_geometry.jl:104-117`) does for atomic point
elements instead of inventing a new path.

Add the include to `src/GeometryOps.jl` directly after the relateng block:

```julia
include("methods/geom_relations/relateng/relate_ng.jl")
include("prepared/prepared.jl")
```

**Step 4: Run test to verify it passes**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: PASS.

**Step 5: Commit**

```bash
git add src/prepared/prepared.jl src/GeometryOps.jl test/prepared/prepared_geometry.jl
git commit -m "Add recursive \`prepare\` builder for prepared geometry"
```

---

### Task 5: `RingEdgeIndex` — per-ring segment extent tree

**Files:**
- Modify: `src/prepared/prepared.jl` (append)
- Modify: `test/prepared/prepared_geometry.jl` (append)

**Step 1: Write the failing test**

```julia
@testset "RingEdgeIndex" begin
    p = prepare(_PG_POLY; preps = (GO.RingEdgeIndex(),))
    @test p.preps === ()                       # consumed at ring level, not polygon level
    ring = GI.getring(p, 1)
    prep = get(ring, SpatialEdgeIndexLike())
    @test prep isa GO.SpatialEdgeIndex
    tree = prep.tree
    @test tree isa GO.NaturalIndexing.NaturalIndex

    # the tree indexes segments in ring order: query a box covering only segment 1
    # (outer ring segment 1 runs (0,0)->(10,0))
    hits = GO.SpatialTreeInterface.query(tree, Extents.Extent(X = (4.0, 5.0), Y = (-0.1, 0.1)))
    @test hits == [1]

    # spherical: builds 3D arc-extent trees without error
    ps = prepare(_PG_SPH_POLY; manifold = GO.Spherical(), preps = (GO.RingEdgeIndex(),))
    sprep = get(GI.getring(ps, 1), SpatialEdgeIndexLike())
    @test sprep isa GO.SpatialEdgeIndex
    @test GO.SpatialTreeInterface.node_extent(sprep.tree) isa Extents.Extent{(:X, :Y, :Z)}

    # manifold-checked retrieval
    @test getprep(GO.Planar(), GI.getring(p, 1), SpatialEdgeIndexLike()) === prep
    @test getprep(GO.Spherical(), GI.getring(p, 1), SpatialEdgeIndexLike()) === nothing
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: FAIL — `UndefVarError: RingEdgeIndex`.

**Step 3: Write the implementation**

Append to `src/prepared/prepared.jl`:

```julia
"""
    SpatialEdgeIndex(tree)

Built preparation: an extent tree (anything implementing SpatialTreeInterface) over the
segments of the wrapped curve, in traversal order. Capability: `SpatialEdgeIndexLike`.
"""
struct SpatialEdgeIndex{T} <: AbstractPreparation
    tree::T
end
preptrait(::SpatialEdgeIndex) = SpatialEdgeIndexLike()

"""
    RingEdgeIndex(; nodecapacity = 16)

Spec: build a `NaturalIndex` over each linear ring's segment extents (manifold-aware:
3D arc extents on `Spherical`). Replaces the `NaturallyIndexedRing` experiment.
"""
struct RingEdgeIndex <: PreparationSpec
    nodecapacity::Int
end
RingEdgeIndex(; nodecapacity::Integer = 16) = RingEdgeIndex(Int(nodecapacity))
appliesto(::RingEdgeIndex) = GI.LinearRingTrait

function buildprep(spec::RingEdgeIndex, m::Manifold, ring)
    exts = _curve_segment_extents(m, ring)
    isempty(exts) && return nothing
    return SpatialEdgeIndex(NaturalIndex(exts; nodecapacity = spec.nodecapacity))
end

function _curve_segment_extents(m::Manifold, curve)
    exts = _segment_extent_type(m)[]
    n = GI.npoint(curve)
    n < 2 && return exts
    sizehint!(exts, n - 1)
    prev = _to_kernel_point(m, GI.getpoint(curve, 1))
    for i in 2:n
        cur = _to_kernel_point(m, GI.getpoint(curve, i))
        push!(exts, _segment_extent(m, prev, cur))
        prev = cur
    end
    return exts
end
```

(Use the actual helper names confirmed in Task 4 Step 0.)

Also export the built-prep types: change the export line at the top of
`src/prepared/prepared.jl` to:

```julia
export RingEdgeIndex, ChildTree, PointInArea
export SpatialEdgeIndex, SpatialIndex, PointInAreaIndex
```

(`ChildTree`/`SpatialIndex`/`PointInArea`/`PointInAreaIndex` come in Tasks 6–7; exporting
ahead is fine since the names will exist by the time the module loads — if you prefer,
add each name in its own task.)

**Step 4: Run test to verify it passes**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: PASS.

**Step 5: Commit**

```bash
git add src/prepared/prepared.jl test/prepared/prepared_geometry.jl
git commit -m "Add \`RingEdgeIndex\` preparation"
```

---

### Task 6: `ChildTree` — extent tree over child elements

**Files:**
- Modify: `src/prepared/prepared.jl` (append)
- Modify: `test/prepared/prepared_geometry.jl` (append)

**Step 1: Write the failing test**

```julia
@testset "ChildTree" begin
    # on a multipolygon: consumed at the MP node (topmost-wins), indexing polygon extents
    mp = prepare(_PG_MP; preps = (GO.ChildTree(),))
    prep = get(mp, SpatialIndexLike())
    @test prep isa GO.SpatialIndex
    # polygon 2 lives at x in (20, 25): a query box there hits only child 2
    @test GO.SpatialTreeInterface.query(prep.tree, Extents.Extent(X = (21.0, 22.0), Y = (0.0, 1.0))) == [2]
    # children did not also get trees
    @test get(GI.getgeom(mp, 1), SpatialIndexLike()) === nothing

    # on a bare polygon: indexes ring extents
    p = prepare(_PG_POLY; preps = (GO.ChildTree(),))
    rprep = get(p, SpatialIndexLike())
    @test rprep isa GO.SpatialIndex
    # the hole (ring 2) occupies (4..6, 4..6)
    @test GO.SpatialTreeInterface.query(rprep.tree, Extents.Extent(X = (4.5, 5.5), Y = (4.5, 5.5))) == [1, 2]
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: FAIL — `UndefVarError: ChildTree`.

**Step 3: Write the implementation**

Append to `src/prepared/prepared.jl`:

```julia
"""
    SpatialIndex(tree)

Built preparation: an extent tree over the wrapped geometry's *child elements* (rings of a
polygon, polygons of a multipolygon, members of a collection), in child order.
Capability: `SpatialIndexLike`.
"""
struct SpatialIndex{T} <: AbstractPreparation
    tree::T
end
preptrait(::SpatialIndex) = SpatialIndexLike()

"""
    ChildTree(; nodecapacity = 32)

Spec: build a `NaturalIndex` over child-element extents. Attaches (topmost-wins) to
polygons, multi-geometries, and geometry collections.
"""
struct ChildTree <: PreparationSpec
    nodecapacity::Int
end
ChildTree(; nodecapacity::Integer = 32) = ChildTree(Int(nodecapacity))
appliesto(::ChildTree) = Union{GI.PolygonTrait, GI.AbstractMultiPolygonTrait,
    GI.AbstractMultiCurveTrait, GI.GeometryCollectionTrait}

function buildprep(spec::ChildTree, m::Manifold, geom)
    exts = _segment_extent_type(m)[]
    for c in GI.getgeom(geom)
        GI.isempty(c) && continue
        push!(exts, c isa Prepared ? c.extent : rk_interaction_bounds(m, c))
    end
    isempty(exts) && return nothing
    return SpatialIndex(NaturalIndex(exts; nodecapacity = spec.nodecapacity))
end
```

(Note: `_segment_extent_type(m)` doubles as "the manifold's extent type" — segment and
interaction bounds share it. If empty children make child indices ambiguous for a caller,
that caller filters; v1 keeps build simple and skips empties.)

**Step 4: Run test to verify it passes**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: PASS.

**Step 5: Commit**

```bash
git add src/prepared/prepared.jl test/prepared/prepared_geometry.jl
git commit -m "Add \`ChildTree\` preparation"
```

---

### Task 7: `PointInArea` — planar point-in-polygonal-area locator

**Files:**
- Modify: `src/prepared/prepared.jl` (append)
- Modify: `test/prepared/prepared_geometry.jl` (append)

**Step 1: Write the failing test**

```julia
@testset "PointInArea" begin
    p = prepare(_PG_POLY; preps = (GO.PointInArea(),))
    prep = get(p, PointInAreaLike())
    @test prep isa GO.PointInAreaIndex
    # ground truth against the unindexed locator: interior, hole, boundary, exterior
    @test GO.locate(prep.locator, (1.0, 1.0)) == GO.LOC_INTERIOR
    @test GO.locate(prep.locator, (5.0, 5.0)) == GO.LOC_EXTERIOR   # inside the hole
    @test GO.locate(prep.locator, (5.0, 0.0)) == GO.LOC_BOUNDARY
    @test GO.locate(prep.locator, (11.0, 5.0)) == GO.LOC_EXTERIOR

    # multipolygon input: consumed at the MP node (relateng's element granularity)
    mp = prepare(_PG_MP; preps = (GO.PointInArea(),))
    @test get(mp, PointInAreaLike()) isa GO.PointInAreaIndex
    @test get(GI.getgeom(mp, 1), PointInAreaLike()) === nothing

    # spherical: a clear, early error (Task 16 of the spherical plan was skipped)
    @test_throws ArgumentError prepare(_PG_SPH_POLY;
        manifold = GO.Spherical(), preps = (GO.PointInArea(),))
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: FAIL — `UndefVarError: PointInArea`.

**Step 3: Write the implementation**

Append to `src/prepared/prepared.jl`:

```julia
"""
    PointInAreaIndex(locator)

Built preparation: an `IndexedPointInAreaLocator` (sorted leaves) over the wrapped
polygonal geometry. Capability: `PointInAreaLike`. `Planar` only.
"""
struct PointInAreaIndex{L} <: AbstractPreparation
    locator::L
end
preptrait(::PointInAreaIndex) = PointInAreaLike()

"""
    PointInArea(; exact = True())

Spec: build an indexed point-in-area locator for a polygon or multipolygon. `Planar` only —
the spherical point-in-area index is deliberately unimplemented (relateng falls back to an
O(n) `rk_point_in_ring` walk on `Spherical`).
"""
struct PointInArea{E} <: PreparationSpec
    exact::E
end
PointInArea(; exact = True()) = PointInArea(exact)
appliesto(::PointInArea) = Union{GI.PolygonTrait, GI.MultiPolygonTrait}

buildprep(spec::PointInArea, m::Planar, geom) =
    PointInAreaIndex(IndexedPointInAreaLocator(m, geom; exact = spec.exact, sort_leaves = true))
buildprep(::PointInArea, m::Manifold, _) =
    throw(ArgumentError("`PointInArea` requires the `Planar` manifold; got `$(typeof(m))`. \
        The spherical point-in-area index is not implemented — spherical algorithms use an \
        O(n) ring walk instead, so simply omit this spec."))
```

If `GO.locate`/`GO.LOC_INTERIOR` etc. are named differently, check
`src/methods/geom_relations/relateng/indexed_point_in_area.jl:307-316` and
`point_locator.jl` for the actual names (`locate`, `LOC_INTERIOR`, `LOC_BOUNDARY`,
`LOC_EXTERIOR` are the expected ones) and use those in the test.

**Step 4: Run test to verify it passes**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: PASS.

**Step 5: Commit**

```bash
git add src/prepared/prepared.jl test/prepared/prepared_geometry.jl
git commit -m "Add \`PointInArea\` preparation"
```

---

### Task 8: RelateNG seam 1 — trust `Prepared` extents

**Files:**
- Modify: `src/methods/geom_relations/relateng/relate_geometry.jl:136` (the 1-arg
  `_relate_cache_extents` dispatcher)
- Modify: `test/prepared/prepared_geometry.jl` (append)

**Step 1: Write the failing test**

```julia
@testset "relateng seam: Prepared extents are trusted" begin
    m = GO.Planar()
    pp = prepare(_PG_POLY; preps = (GO.RingEdgeIndex(),))
    rg = GO.RelateGeometry(m, pp; exact = GO.True())
    # matched manifold: the Prepared tree is used as-is (no rewrap)
    @test rg.geom === pp

    # mismatched manifold: falls back to the rewrap path (correctness over reuse)
    ps = prepare(_PG_SPH_POLY)   # Planar-prepared (2D extents)
    rgs = GO.RelateGeometry(GO.Spherical(), ps; exact = GO.True())
    @test rgs.geom !== ps
    @test GI.extent(rgs.geom) isa Extents.Extent{(:X, :Y, :Z)}
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: FAIL — `rg.geom === pp` is false (relateng rewraps the Prepared tree today).

**Step 3: Write the implementation**

In `relate_geometry.jl`, replace the 1-arg dispatcher

```julia
_relate_cache_extents(m::Manifold, geom) = _relate_cache_extents(m, GI.trait(geom), geom)
```

with

```julia
#-- a Prepared input built on the same manifold already carries exactly the extents this
#-- rebuild would compute (rk_interaction_bounds at every level) — use it as-is. A
#-- mismatched manifold falls through to the rewrap: wrong-manifold extents must not leak.
function _relate_cache_extents(m::Manifold, geom)
    geom isa Prepared && GeometryOpsCore.manifold(geom) === m && return geom
    return _relate_cache_extents(m, GI.trait(geom), geom)
end
```

Do NOT extend `_has_stored_extent` — for mismatched manifolds the trait-dispatched rewrap
must treat `Prepared` nodes as *not* having usable stored extents, which the existing
`WrapperGeometry`-only check already guarantees.

**Step 4: Run tests to verify they pass**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: PASS.

Run: `julia --project=test -e 'include("test/methods/relateng/relate_geometry.jl")'`
Expected: PASS (no regression in the extent-cache tests).

**Step 5: Commit**

```bash
git add src/methods/geom_relations/relateng/relate_geometry.jl test/prepared/prepared_geometry.jl
git commit -m "Trust matched-manifold \`Prepared\` extents in relateng"
```

---

### Task 9: RelateNG seam 2 — reuse `PointInAreaLike` in the point locator

**Files:**
- Modify: `src/methods/geom_relations/relateng/point_locator.jl` (the
  `locate_on_polygonal` / `_get_poly_locator` region, lines ~522–556)
- Modify: `test/prepared/prepared_geometry.jl` (append)

**Step 1: Write the failing test**

```julia
@testset "relateng seam: PointInArea reuse" begin
    m = GO.Planar()
    pp = prepare(_PG_POLY; preps = (GO.PointInArea(),))
    prep = get(pp, PointInAreaLike())

    # unprepared locator over a Prepared element uses the preparation immediately
    loc = GO.RelatePointLocator(m, pp; exact = GO.True())
    @test GO.locate(loc, (1.0, 1.0)) == GO.LOC_INTERIOR   # forces the polygonal path
    @test loc.poly_locator[1] === prep.locator             # identity: reused, not rebuilt

    # and in prepared mode too
    ploc = GO.RelatePointLocator(m, pp; exact = GO.True(), is_prepared = true)
    GO.locate(ploc, (1.0, 1.0))
    @test ploc.poly_locator[1] === prep.locator
end
```

If `locate(loc::RelatePointLocator, p)` is not the public entry (check
`point_locator.jl` for the actual query function relateng calls — e.g.
`locate_with_dim` or similar), call whatever function routes through
`locate_on_polygonal`, or call `GO.locate_on_polygonal(loc, (1.0,1.0), false, nothing, 1)`
directly with its real signature.

**Step 2: Run test to verify it fails**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: FAIL — identity test fails (a fresh `IndexedPointInAreaLocator` is built), or in
the unprepared case the direct ring walk is used for the first 8 queries so
`poly_locator[1] === nothing`.

**Step 3: Write the implementation**

In `point_locator.jl`, modify the planar branch of `locate_on_polygonal` (lines ~528–536)
so a preparation counts as "use the index now":

```julia
    if loc.m isa Planar
        use_index = loc.is_prepared ||
            getprep(loc.m, polygonal, PointInAreaLike()) !== nothing
        if !use_index
            count = (loc.poly_query_count[index] += Int32(1))
            use_index = count > _LAZY_INDEX_QUERY_THRESHOLD
        end
        if use_index
            return locate(_get_poly_locator(loc, index), p)
        end
    end
```

and modify `_get_poly_locator` (lines ~548–556) to reuse the prepared locator when its
manifold and `exact` types match the locator slot:

```julia
function _get_poly_locator(loc::RelatePointLocator, index::Int)
    locator = loc.poly_locator[index]
    if locator === nothing
        polygonal = loc.polygons[index]
        prep = getprep(loc.m, polygonal, PointInAreaLike())
        locator = if prep isa PointInAreaIndex &&
                prep.locator isa IndexedPointInAreaLocator{typeof(loc.m), typeof(loc.exact)}
            prep.locator
        else
            IndexedPointInAreaLocator(loc.m, polygonal;
                exact = loc.exact, sort_leaves = loc.is_prepared)
        end
        loc.poly_locator[index] = locator
    end
    return locator
end
```

**Step 4: Run tests to verify they pass**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: PASS.

Run: `julia --project=test -e 'include("test/methods/relateng/point_locator.jl")'`
Expected: PASS (no regression).

**Step 5: Commit**

```bash
git add src/methods/geom_relations/relateng/point_locator.jl test/prepared/prepared_geometry.jl
git commit -m "Reuse \`PointInAreaLike\` preparations in the relate point locator"
```

---

### Task 10: End-to-end equality — `relate` over `Prepared` inputs

**Files:**
- Modify: `test/prepared/prepared_geometry.jl` (append)

**Step 1: Write the test (expected to pass — this is the safety net, watch it actually run)**

```julia
@testset "relate: plain == Prepared, planar & spherical" begin
    line = GI.LineString([(-1.0, -1.0), (5.0, 5.0), (12.0, 5.0)])
    pt   = GI.Point(5.0, 1.0)
    pairs = [
        (_PG_POLY, _PG_POLY2),   # disjoint
        (_PG_POLY, GI.Polygon([GI.LinearRing([(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 5.0)])])), # overlap
        (_PG_POLY, line),
        (_PG_POLY, pt),
        (_PG_MP, _PG_POLY),
    ]
    alg = GO.RelateNG()
    for (a, b) in pairs
        expected = string(GO.relate(alg, a, b))
        pa = prepare(a; preps = (GO.RingEdgeIndex(), GO.ChildTree(), GO.PointInArea()))
        # plain path with a Prepared input
        @test string(GO.relate(alg, pa, b)) == expected
        # algorithm-prepared path on top of a Prepared input
        pr = GO.prepare(alg, pa)
        @test string(GO.relate(pr, b)) == expected
        # B side Prepared as well
        @test string(GO.relate(alg, pa, prepare(b; preps = (GO.RingEdgeIndex(),)))) == expected
    end

    # spherical (no PointInArea on the sphere)
    salg = GO.RelateNG(; manifold = GO.Spherical())
    sb = GI.Polygon([GI.LinearRing([(15.0, 45.0), (25.0, 45.0), (25.0, 55.0), (15.0, 55.0), (15.0, 45.0)])])
    sexpected = string(GO.relate(salg, _PG_SPH_POLY, sb))
    spa = prepare(_PG_SPH_POLY; manifold = GO.Spherical(), preps = (GO.RingEdgeIndex(),))
    @test string(GO.relate(salg, spa, sb)) == sexpected
    @test string(GO.relate(GO.prepare(salg, spa), sb)) == sexpected
end
```

**Step 2: Run the test**

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: PASS. If any pair disagrees, STOP and debug with
superpowers:systematic-debugging — an equality failure here means a seam changed
semantics, which is never acceptable; do not "fix" the test.

**Step 3: Run the relateng sweeps to confirm no regression**

Run: `julia --project=test -e 'include("test/methods/relateng/relate_ng.jl")'`
Expected: PASS.
Run: `julia --project=test -e 'include("test/methods/relateng/spherical_end_to_end.jl")'`
Expected: PASS.

**Step 4: Commit**

```bash
git add test/prepared/prepared_geometry.jl
git commit -m "Test \`relate\` equality over \`Prepared\` inputs"
```

---

### Task 11: Delete the `NaturallyIndexedRing` experiment

It has zero call sites outside its own file (verified by grep during planning) and its
docstring says it exists only to prototype what `Prepared` now provides.

**Files:**
- Modify: `src/utils/NaturalIndexing.jl`

**Step 1: Delete**

In `src/utils/NaturalIndexing.jl`:
1. Delete the `NaturallyIndexedRing` struct, its constructors, its GI methods, and
   `prepare_naturally` (the block at lines ~199–243, from the section comment above the
   struct through end of file — read the file first to get the exact block).
2. Change the export line (line ~10) from
   `export NaturalIndex, NaturallyIndexedRing, prepare_naturally` to
   `export NaturalIndex`.
3. Delete the module-header line
   `import ..GeometryOps as GO # TODO: only needed for NaturallyIndexedRing, remove when that is removed.`
   (line ~8) — that TODO is now done.

**Step 2: Verify**

```bash
grep -rn "NaturallyIndexedRing\|prepare_naturally" src/ test/ ext/ 2>/dev/null
```
Expected: no output.

Run: `julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'`
Expected: PASS (package still loads; `RingEdgeIndex` is the replacement).

**Step 3: Commit**

```bash
git add src/utils/NaturalIndexing.jl
git commit -m "Remove \`NaturallyIndexedRing\` in favor of \`RingEdgeIndex\` preparation"
```

---

### Task 12: Versions, design-doc amendments

**Files:**
- Modify: `GeometryOpsCore/Project.toml` (`version = "0.1.10"` → `"0.1.11"`)
- Modify: `Project.toml` (GO `version = "0.1.40"` → `"0.1.41"`; compat
  `GeometryOpsCore = "=0.1.10"` → `"=0.1.11"`)
- Modify: `docs/plans/2026-07-01-prepared-geometry-design.md`

**Step 1: Bump versions** as above.

**Step 2: Amend the design doc.** Update the status line to
`**Status: ACCEPTED — implemented on branch \`prepared-geometry\`; deviations below.**`
and append this section verbatim:

```markdown
## Amendments (discovered during implementation, 2026-07-01)

- **`SegmentIndex` dropped from v1.** A flattened whole-geometry edge tree carries an
  `owners` table whose indices are only meaningful relative to relateng's
  `extract_segment_strings` traversal order; building it independently risks silent
  owner mismatches. RelateNG keeps building its own `PreparedEdgeIndex` inside
  `prepare(::RelateNG, a)`. Revisit when a second consumer needs flat segment trees.
- **Clipping `AutoAccelerator` seam deferred** for the same reason (clipping indexes
  flattened edge tables, not child extents).
- **`tuples` keyword deferred** — compose as `prepare(GO.tuples(geom); ...)`.
- **`build` is spelled `buildprep`**; `appliesto`/`buildprep` are unexported (providers
  extend via `import GeometryOpsCore: appliesto, buildprep`).
- **Topmost-wins attachment:** a spec is consumed at the highest matching level and not
  offered further down (so `ChildTree` on a MultiPolygon indexes polygons at the top
  node only). Multi-level trees = call `prepare` with the spec on the sub-geometry.
- **Mismatched-manifold `Prepared` inputs to relateng** take the ordinary rewrap path
  (correct extents rebuilt; preparations shadowed) rather than erroring.
```

**Step 3: Verify the workspace still resolves**

Run: `julia --project=test -e 'import Pkg; Pkg.resolve(); using GeometryOps'`
Expected: resolves and loads without error.

**Step 4: Commit**

```bash
git add GeometryOpsCore/Project.toml Project.toml docs/plans/2026-07-01-prepared-geometry-design.md
git commit -m "Bump GeometryOpsCore to 0.1.11 and record prepared-geometry design amendments"
```

---

### Task 13: Final verification sweep

**Step 1: Run the targeted suites** (each from repo root; a few minutes total):

```bash
julia --project=test -e 'include("test/core/preparations.jl")'
julia --project=test -e 'include("test/prepared/prepared_geometry.jl")'
julia --project=test -e 'include("test/methods/relateng/relate_ng.jl")'
julia --project=test -e 'include("test/methods/relateng/relate_geometry.jl")'
julia --project=test -e 'include("test/methods/relateng/point_locator.jl")'
julia --project=test -e 'include("test/methods/relateng/kernel_conformance.jl")'
julia --project=test -e 'include("test/methods/relateng/spherical_end_to_end.jl")'
```

Expected: all PASS. (`kernel_conformance.jl` is ~77k assertions; it runs in seconds once
precompiled.)

**Step 2: Kick off the slow differential suite in the background** (final gate, ~25 min —
do not block on it interactively):

```bash
julia --project=test -e 'include("test/methods/relateng/xml_suite.jl")'
```

Run it as a background task and check the result when it finishes. Expected: PASS.

**Step 3:** Use superpowers:verification-before-completion, then report done with the
evidence (test counts), and hand off per superpowers:finishing-a-development-branch.
