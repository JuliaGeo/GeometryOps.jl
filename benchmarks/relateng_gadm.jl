# # RelateNG GADM real-data benchmarks
#
#=
Benchmarks for the RelateNG DE-9IM engine on **GADM** full-resolution country
boundaries — the high-resolution complement to `benchmarks/relateng_realdata.jl`
(Natural Earth). GADM level-0 polygons are one to two orders of magnitude denser
than Natural Earth (Canada is ~3.9M vertices / 24k rings here, vs ~68k / 412 at
NE 10m), so this file exercises `prepare`, its indexed point-in-area locators,
and the spherical validation self-join at real production scale.

This is a **separate, independently runnable file** because GADM downloads are
heavy: GADM.jl fetches a per-country GeoPackage on first use (tens of MB each),
so the first ever run needs network access and writes to the DataDeps cache;
cached runs are offline. **Pre-download** the countries used below before a timed
run — the file loads:

    CAN                          (point-in-area target, largest sampled country)
    EGY SDN FRA ITA GRC TUR      (neighbour-pair `touches`)

e.g. `julia --project=docs -e 'import GADM; foreach(GADM.get, ["CAN","EGY","SDN","FRA","ITA","GRC","TUR"])'`.

Three workload groups, each a table of per-operation medians:

1. Point-in-area at GADM scale on the largest sampled country (Canada): the
   planar `intersects` mix — unprepared, extent-stamped (`GO.tuples(x;
   calc_extent = true)`), prepared (`GO.prepare`), LibGEOS plain and prepared —
   per query, plus the `GO.prepare` build split into *index build*
   (`validate = false`) vs the default validating build, so the validation
   overhead at ~4M vertices is a recorded number.
2. Spherical at GADM scale on the same country: prepared point-in-area through
   the longitude-interval indexed locator, and the prepared build with and
   without the ring self-crossing validation — the locator and validator's
   first serious-scale exercise.
3. `touches` on real land-border neighbour pairs, planar RelateNG vs LibGEOS.

Run with `julia --project=docs benchmarks/relateng_gadm.jl`. Prints comparison
tables; no CI gating (representative output is in the comment block at the
bottom). Everything is seeded (`MersenneTwister(7)`), so reruns measure
identical workloads. Total runtime is ~5–10 min on a warm data cache (excluding
package precompilation), dominated by the ~6 s spherical `prepare` builds.
=#

import GeometryOps as GO,
    GeoInterface as GI,
    LibGEOS as LG
import Extents
import GADM
using Chairmarks
using Printf
using Random

# As in `benchmarks/relateng_realdata.jl`: give each package its native geometry.
lg_and_go(geometry) = (GI.convert(LG, geometry), GO.tuples(geometry))

prettytime(s) =
    s < 1e-6 ? @sprintf("%8.1f ns", s * 1e9) :
    s < 1e-3 ? @sprintf("%8.1f μs", s * 1e6) :
    s < 1.0  ? @sprintf("%8.1f ms", s * 1e3) :
               @sprintf("%8.2f s ", s)

median_time(trial) = Chairmarks.median(trial).time
# Each evaluation sweeps a workload of `n` operations; report per-op time.
per_op(trial, n) = median_time(trial) / n

function print_table(title, firstcol, colnames, rows)
    printstyled(title; color = :green, bold = true)
    println()
    @printf("%-30s", firstcol)
    foreach(c -> @printf(" │ %18s", c), colnames)
    println()
    println("─"^(30 + 21 * length(colnames)))
    for (label, times) in rows
        @printf("%-30s", label)
        foreach(t -> @printf(" │ %18s", prettytime(t)), times)
        println()
    end
    println()
end

# Load a GADM level-0 country by ISO-3 code and return its (Multi)Polygon.
# GADM.jl returns a Tables.jl feature collection; level-0 is a single feature.
function gadm_geom(code)
    tbl = try
        GADM.get(code)
    catch err
        error("""
            Could not load GADM country `$code`.
            GADM.jl downloads a per-country GeoPackage on first use — this
            machine either needs network access or a pre-warmed DataDeps cache.
            Underlying error: $(sprint(showerror, err))""")
    end
    return GI.geometry(GI.getfeature(tbl, 1))
end

const ALG  = GO.RelateNG()
const SALG = GO.RelateNG(; manifold = GO.Spherical())

# ## Workload 1: point-in-area on the largest sampled country (Canada)

big_raw = gadm_geom("CAN")
lg_big, go_big = lg_and_go(big_raw)
go_big_ext = GO.tuples(big_raw; calc_extent = true)  # extent-stamped variant
nrings = GO.applyreduce(x -> 1, +, GI.LinearRingTrait, go_big; init = 0)

# Seed-7 point mix: uniform points in Canada's extent (about half hit).
ext = GI.extent(go_big)
rng = MersenneTwister(7)
pts = [GI.Point((ext.X[1] + rand(rng) * (ext.X[2] - ext.X[1]),
                 ext.Y[1] + rand(rng) * (ext.Y[2] - ext.Y[1]))) for _ in 1:10_000]
lg_pts = [LG.readgeom("POINT($(GI.x(p)) $(GI.y(p)))") for p in pts]

prep = GO.prepare(ALG, go_big)
lg_prep = LG.prepareGeom(lg_big)

nhit = count(p -> GO.relate_predicate(prep, GO.pred_intersects(), p), pts)
@printf("Point-in-area target: Canada — %d vertices, %d rings; point mix: %d/%d hits\n\n",
    GI.npoint(go_big), nrings, nhit, length(pts))

# Unprepared point queries reconvert and re-extent the whole geometry per call,
# so at ~4M vertices they are milliseconds — a small point count suffices.
t_unprep  = @be count(p -> GO.intersects($ALG, $go_big, p), $(pts[1:20])) seconds=1
t_stamped = @be count(p -> GO.intersects($ALG, $go_big_ext, p), $(pts[1:50])) seconds=1
t_prep    = @be count(p -> GO.relate_predicate($prep, GO.pred_intersects(), p), $pts) seconds=1
t_lg      = @be count(p -> LG.intersects($lg_big, p), $(lg_pts[1:1000])) seconds=1
t_lg_prep = @be count(p -> LG.intersects($lg_prep, p), $lg_pts) seconds=1
t_build   = @be GO.prepare($ALG, $go_big) seconds=2

pq_unprep, pq_prep = per_op(t_unprep, 20), per_op(t_prep, length(pts))
print_table("point-in-area intersects (Canada, GADM full res, seed-7 point mix, per query)",
    "workload",
    ["unprepared", "extent-stamped", "prepared", "LibGEOS", "LibGEOS prepared"],
    ["Canada, per point" =>
        [pq_unprep, per_op(t_stamped, 50), pq_prep,
         per_op(t_lg, 1000), per_op(t_lg_prep, length(pts))]])

build = median_time(t_build)
@printf("GO.prepare build: %s → amortized against unprepared after ~%.1f queries\n\n",
    strip(prettytime(build)), build / (pq_unprep - pq_prep))

# Build cost split: index build (`validate = false`) vs the default validating
# build (the ring self-crossing self-join). Planar `prepare` defaults to
# `validate = false`; timing both isolates the validation overhead at scale.
t_build_noval = @be GO.prepare($ALG, $go_big; validate = false) seconds=2
t_build_val   = @be GO.prepare($ALG, $go_big; validate = true)  seconds=2
print_table("planar prepare build cost split (Canada, per build)",
    "build",
    ["index only (validate=false)", "validating (validate=true)"],
    ["GO.prepare" => [median_time(t_build_noval), median_time(t_build_val)]])
@printf("planar validation overhead: %s (%.0f%% of the validating build)\n\n",
    strip(prettytime(median_time(t_build_val) - median_time(t_build_noval))),
    100 * (median_time(t_build_val) - median_time(t_build_noval)) / median_time(t_build_val))

# ## Workload 2: spherical point-in-area + build on the same country
#
# Prepared spherical point-in-area runs the indexed locator (longitude-interval
# edge index + meridian-arc crossing parity). Spherical `prepare` defaults to
# `validate = true` (the ring self-crossing check that motivates default
# validation); the build split records what that validation costs at ~4M edges.

sprep = GO.prepare(SALG, go_big; validate = false)
nhit_s = count(p -> GO.relate_predicate(sprep, GO.pred_intersects(), p), pts[1:200])
t_s_prep = @be count(p -> GO.relate_predicate($sprep, GO.pred_intersects(), p), $(pts[1:200])) seconds=1
t_s_build_noval = @be GO.prepare($SALG, $go_big; validate = false) seconds=1
t_s_build_val   = @be GO.prepare($SALG, $go_big; validate = true)  seconds=1

print_table("spherical RelateNG at GADM scale (Canada, per op)",
    "workload",
    ["Spherical", "Planar"],
    ["point-in-area prepared, per point" => [per_op(t_s_prep, 200), pq_prep],
     "prepare build (index only)"        => [median_time(t_s_build_noval), median_time(t_build_noval)],
     "prepare build (validating)"        => [median_time(t_s_build_val),   median_time(t_build_val)]])
@printf("spherical prepared point mix: %d/200 hits; spherical validation overhead: %s (%.0f%% of the validating build)\n\n",
    nhit_s,
    strip(prettytime(median_time(t_s_build_val) - median_time(t_s_build_noval))),
    100 * (median_time(t_s_build_val) - median_time(t_s_build_noval)) / median_time(t_s_build_val))

# ## Workload 3: `touches` on real land-border neighbour pairs

const PAIRS = [("EGY", "SDN"), ("FRA", "ITA"), ("GRC", "TUR")]
pair_geoms = Dict(c => lg_and_go(gadm_geom(c)) for c in unique(Iterators.flatten(PAIRS)))

rows = Pair{String, Vector{Float64}}[]
for (ca, cb) in PAIRS
    lga, goa = pair_geoms[ca]
    lgb, gob = pair_geoms[cb]
    pt_ng = @be GO.touches($ALG, $goa, $gob) seconds=1
    pt_lg = @be LG.touches($lga, $lgb) seconds=1
    push!(rows, "$(ca)–$(cb) ($(GI.npoint(goa))/$(GI.npoint(gob)) verts)" =>
        [median_time(pt_ng), median_time(pt_lg)])
end
print_table("touches on real land-border pairs (GADM full res, per predicate)",
    "pair", ["RelateNG", "LibGEOS"], rows)

#=
Representative output (2026-07-15, Apple M4 Pro, macOS — Darwin 25.5.0;
Julia 1.12.6, GEOS 3.14.1, GADM.jl 1.2.0 (GADM 4.1), GeometryOps @ 2e6dd8aef;
`julia --project=docs benchmarks/relateng_gadm.jl`, ~4 min wall on a warm data
cache excluding package precompilation):

Point-in-area target: Canada — 3892522 vertices, 24480 rings; point mix: 4676/10000 hits

point-in-area intersects (Canada, GADM full res, seed-7 point mix, per query)
workload                       │         unprepared │     extent-stamped │           prepared │            LibGEOS │   LibGEOS prepared
───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Canada, per point              │             2.9 ms │           550.1 μs │             1.8 μs │           710.2 μs │             1.2 μs

GO.prepare build: 40.1 ms → amortized against unprepared after ~13.7 queries

planar prepare build cost split (Canada, per build)
build                          │ index only (validate=false) │ validating (validate=true)
────────────────────────────────────────────────────────────────────────
GO.prepare                     │            35.7 ms │           973.3 ms

planar validation overhead: 937.6 ms (96% of the validating build)

spherical RelateNG at GADM scale (Canada, per op)
workload                       │          Spherical │             Planar
────────────────────────────────────────────────────────────────────────
point-in-area prepared, per point │             2.3 μs │             1.8 μs
prepare build (index only)     │            5.47 s  │            35.7 ms
prepare build (validating)     │            7.06 s  │           973.3 ms

spherical prepared point mix: 104/200 hits; spherical validation overhead: 1.58 s (22% of the validating build)

touches on real land-border pairs (GADM full res, per predicate)
pair                           │           RelateNG │            LibGEOS
────────────────────────────────────────────────────────────────────────
EGY–SDN (117090/49006 verts)   │             1.4 ms │             1.9 ms
FRA–ITA (216353/303163 verts)  │             6.1 ms │            10.8 ms
GRC–TUR (436532/230727 verts)  │            31.1 ms │            55.9 ms

Reading notes:

- Point-in-area at ~4M vertices: the prepared indexed locator is the whole
  story. Prepared RelateNG (1.8 μs) tracks prepared LibGEOS (1.2 μs) to within
  ~1.5x and is ~1600x faster than unprepared RelateNG (2.9 ms) and ~390x
  faster than plain LibGEOS (710 μs) — at this scale the per-call envelope /
  reconversion cost that unprepared and plain-GEOS pay every query dwarfs the
  point test. `GO.prepare` (35.7 ms index build) pays for itself after ~14
  unprepared-equivalent queries. Extent-stamping the input
  (`GO.tuples(x; calc_extent = true)`) removes the per-call extent pass and
  cuts unprepared from 2.9 ms to 550 μs, but the indexed prepared path is three
  orders of magnitude below either.
- Planar validation at scale: the ring self-crossing self-join over ~3.9M edges
  is 938 ms — 96% of the validating build. Planar `prepare` therefore leaves it
  off by default (index build alone is 36 ms); opting in (`validate = true`) is
  a ~1 s one-time cost that a planar-invalid input would otherwise surface only
  as a wrong answer.
- Spherical: the longitude-interval indexed locator makes prepared spherical
  point-in-area 2.3 μs/query — essentially at parity with planar (1.8 μs) even
  at 3.9M vertices, confirming the locator's per-query cost scales with
  edges-crossing-the-meridian, not total edges. The spherical build is the real
  cost: 5.47 s index build (dominated by per-vertex kernel-point conversion of
  3.9M vertices) plus 1.58 s validation (22% of the 7.06 s validating build).
  Spherical `prepare` validates by default because the silent-class defect it
  catches inverts containment; at GADM scale that safety costs ~1.6 s on top of
  a build already dominated by index construction.
- `touches` on real land borders: RelateNG is faster than LibGEOS at GADM full
  resolution (1.4 vs 1.9 ms, 6.1 vs 10.8 ms, 31.1 vs 55.9 ms) — the reverse of
  the Natural Earth 10m result (~1.2x behind), because full-res borders resolve
  the predicate from the boundary node topology before escalation, where
  RelateNG's edge indexing amortizes better than GEOS's per-call relate.
=#
