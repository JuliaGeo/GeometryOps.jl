# # RelateNG real-data benchmarks
#
#=
Benchmarks for the RelateNG DE-9IM engine on real Natural Earth data
(complementing `benchmarks/relateng.jl`, which uses synthetic random
polygons). Data comes from NaturalEarth.jl — the artifacts are downloaded
and cached on first use, so the first ever run needs network access; cached
runs are offline.

Five workload groups, each printed as a table of per-operation medians:

1. Point-in-area `intersects` on the largest 10m country (Canada) over a
   seeded uniform point mix in its extent: unprepared, unprepared with
   extent-stamped input (`GO.tuples(x; calc_extent = true)` — the per-call
   extent pass dominates unprepared point queries, and pre-stamped extents
   short-circuit it; this column tracks a planned optimization), prepared
   (`GO.prepare`), and LibGEOS plain/prepared. Plus `GO.prepare` build cost
   and the break-even query count.
2. Pairwise country `intersects` at 110m on a seeded sample of ordered
   pairs: RelateNG vs the old GO per-pair processors vs LibGEOS.
3. `touches` and full `relate` on real border-sharing neighbor pairs at 10m
   (the extent-densest 24 countries; pairs that actually touch per LibGEOS).
4. Rivers x countries `intersects`/`crosses` at 10m on a seeded sample of
   extent-intersecting pairs (`crosses` exercises the self-noding path).
5. Spherical variants (`GO.RelateNG(; manifold = GO.Spherical())`) of the
   point-in-area and pairwise workloads, against their planar equivalents.

Run with `julia --project=docs benchmarks/relateng_realdata.jl`. Prints
comparison tables; no CI gating (representative output is recorded in the
comment block at the bottom of this file). Everything is seeded — reruns
measure identical workloads.
=#

import GeometryOps as GO,
    GeoInterface as GI,
    LibGEOS as LG
import Extents
using NaturalEarth, GeoJSON
using Chairmarks
using Printf
using Random

# As in `benchmarks.jl`: give each package its native geometry.
lg_and_go(geometry) = (GI.convert(LG, geometry), GO.tuples(geometry))

const LG_CTX = LG.get_global_context()
# LibGEOS has no high-level `relate`; the generated wrapper returns the
# DE-9IM string itself (cf. benchmarks/relateng.jl, test/methods/relateng/fuzz.jl).
lg_relate(la, lb) = LG.GEOSRelate_r(LG_CTX, la, lb)

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

# Load a NaturalEarth layer and return (names, geometries), skipping features
# with no geometry and empty geometries (an empty MultiLineString in the 10m
# rivers layer breaks `GO.tuples`).
function ne_geoms(name, scale)
    fc = try
        naturalearth(name, scale)
    catch err
        error("""
            Could not load NaturalEarth layer `$name` at $(scale)m.
            NaturalEarth.jl downloads each layer on first use — this machine
            either needs network access or a pre-warmed artifact cache.
            Underlying error: $(sprint(showerror, err))""")
    end
    names = String[]
    geoms = []
    for f in fc
        g = GeoJSON.geometry(f)
        g === nothing && continue
        GI.npoint(g) == 0 && continue
        props = GeoJSON.properties(f)  # Dict{Symbol, Any}
        push!(names, string(get(props, :NAME, get(props, :name, "?"))))
        push!(geoms, g)
    end
    return names, geoms
end

const ALG = GO.RelateNG()

# ## Workload 1: point-in-area on the largest 10m country

names10, geoms10_raw = ne_geoms("admin_0_countries", 10)
go10 = [GO.tuples(g) for g in geoms10_raw]
exts10 = [GI.extent(g) for g in go10]

i_big = argmax(GI.npoint.(go10))
big_raw = geoms10_raw[i_big]
lg_big, go_big = lg_and_go(big_raw)
go_big_ext = GO.tuples(big_raw; calc_extent = true)  # extent-stamped variant
nrings = GO.applyreduce(x -> 1, +, GI.LinearRingTrait, go_big; init = 0)

# The profiling campaign's seed-7 point mix: uniform points in the target's
# extent (about half hit). Sub-ranges below are prefixes of the same mix.
ext = exts10[i_big]
rng = MersenneTwister(7)
pts = [GI.Point((ext.X[1] + rand(rng) * (ext.X[2] - ext.X[1]),
                 ext.Y[1] + rand(rng) * (ext.Y[2] - ext.Y[1]))) for _ in 1:10_000]
lg_pts = [LG.readgeom("POINT($(GI.x(p)) $(GI.y(p)))") for p in pts]

prep = GO.prepare(ALG, go_big)
lg_prep = LG.prepareGeom(lg_big)

nhit = count(p -> GO.relate_predicate(prep, GO.pred_intersects(), p), pts)
@printf("Point-in-area target: %s — %d vertices, %d rings; point mix: %d/%d hits\n\n",
    names10[i_big], GI.npoint(go_big), nrings, nhit, length(pts))

t_unprep  = @be count(p -> GO.intersects($ALG, $go_big, p), $(pts[1:100])) seconds=1
t_stamped = @be count(p -> GO.intersects($ALG, $go_big_ext, p), $(pts[1:200])) seconds=1
t_prep    = @be count(p -> GO.relate_predicate($prep, GO.pred_intersects(), p), $pts) seconds=1
t_lg      = @be count(p -> LG.intersects($lg_big, p), $(lg_pts[1:1000])) seconds=1
t_lg_prep = @be count(p -> LG.intersects($lg_prep, p), $lg_pts) seconds=1
t_build   = @be GO.prepare($ALG, $go_big) seconds=1

pq_unprep, pq_prep = per_op(t_unprep, 100), per_op(t_prep, length(pts))
print_table("point-in-area intersects (largest 10m country, seed-7 point mix, per query)",
    "workload",
    ["unprepared", "extent-stamped", "prepared", "LibGEOS", "LibGEOS prepared"],
    ["$(names10[i_big]), per point" =>
        [pq_unprep, per_op(t_stamped, 200), pq_prep,
         per_op(t_lg, 1000), per_op(t_lg_prep, length(lg_pts))]])

build = median_time(t_build)
@printf("GO.prepare build: %s → amortized against unprepared after ~%.1f queries\n\n",
    strip(prettytime(build)), build / (pq_unprep - pq_prep))

# ## Workload 2: pairwise country `intersects` at 110m

names110, geoms110_raw = ne_geoms("admin_0_countries", 110)
go110 = [GO.tuples(g) for g in geoms110_raw]
lg110 = [GI.convert(LG, g) for g in geoms110_raw]
n110 = length(go110)

pairs_all = [(i, j) for i in 1:n110 for j in 1:n110 if i != j]
pair_sample = shuffle(Xoshiro(42), pairs_all)[1:3000]
npairhit = count(((i, j),) -> LG.intersects(lg110[i], lg110[j]), pair_sample)
@printf("110m: %d countries (%d vertices); %d sampled ordered pairs, %d intersecting\n\n",
    n110, sum(GI.npoint, go110), length(pair_sample), npairhit)

t_ng  = @be count(((i, j),) -> GO.intersects($ALG, $go110[i], $go110[j]), $pair_sample) seconds=1
t_old = @be count(((i, j),) -> GO.intersects($go110[i], $go110[j]), $pair_sample) seconds=1
t_lgp = @be count(((i, j),) -> LG.intersects($lg110[i], $lg110[j]), $pair_sample) seconds=1

print_table("pairwise intersects (110m countries, seeded pair sample, per pair)",
    "workload",
    ["RelateNG", "GO old", "LibGEOS"],
    ["$(length(pair_sample)) ordered pairs" =>
        per_op.([t_ng, t_old, t_lgp], length(pair_sample))])

# ## Workload 3: `touches` + full `relate` on real 10m border pairs
#
# The campaign's 10m neighbor set: among the 24 countries with the most
# extent-overlaps, the ordered pairs that actually share a border (per
# LibGEOS `touches`). Real borders are near-collinear point soups — this is
# the exact-predicate stress test.

deg = zeros(Int, length(go10))
for i in eachindex(go10), j in eachindex(go10)
    if i != j && Extents.intersects(exts10[i], exts10[j])
        deg[i] += 1
    end
end
top = sortperm(deg; rev = true)[1:24]
cand10 = [(i, j) for i in top for j in top if i != j && Extents.intersects(exts10[i], exts10[j])]
lg10 = Dict(i => GI.convert(LG, geoms10_raw[i]) for i in unique(Iterators.flatten(cand10)))
nbrs = [(i, j) for (i, j) in cand10 if LG.touches(lg10[i], lg10[j])]
@printf("10m neighbor set: %d extent-hit ordered pairs among top-24 countries, %d touching\n\n",
    length(cand10), length(nbrs))

t_t_ng = @be count(((i, j),) -> GO.touches($ALG, $go10[i], $go10[j]), $nbrs) seconds=1
t_t_lg = @be count(((i, j),) -> LG.touches($lg10[i], $lg10[j]), $nbrs) seconds=1
t_r_ng = @be sum(((i, j),) -> length(string(GO.relate($ALG, $go10[i], $go10[j]))), $nbrs) seconds=1
t_r_lg = @be sum(((i, j),) -> length(lg_relate($lg10[i], $lg10[j])), $nbrs) seconds=1

print_table("real border-sharing neighbor pairs (10m, $(length(nbrs)) pairs, per pair)",
    "predicate",
    ["RelateNG", "LibGEOS"],
    ["touches"     => per_op.([t_t_ng, t_t_lg], length(nbrs)),
     "full relate" => per_op.([t_r_ng, t_r_lg], length(nbrs))])

# ## Workload 4: rivers x countries at 10m

rnames, rgeoms_raw = ne_geoms("rivers_lake_centerlines", 10)
gor = [GO.tuples(g) for g in rgeoms_raw]
rexts = [GI.extent(g) for g in gor]
rc = [(ri, ci) for ri in eachindex(gor) for ci in eachindex(go10)
      if Extents.intersects(rexts[ri], exts10[ci])]
rc_sample = shuffle(Xoshiro(42), rc)[1:150]
lgr = Dict(ri => GI.convert(LG, rgeoms_raw[ri]) for ri in unique(first.(rc_sample)))
lgc = Dict(ci => GI.convert(LG, geoms10_raw[ci]) for ci in unique(last.(rc_sample)))
@printf("rivers x countries: %d rivers (%d vertices), %d extent-hit pairs, %d sampled\n\n",
    length(gor), sum(GI.npoint, gor), length(rc), length(rc_sample))

t_i_ng = @be count(((ri, ci),) -> GO.intersects($ALG, $gor[ri], $go10[ci]), $rc_sample) seconds=1
t_i_lg = @be count(((ri, ci),) -> LG.intersects($lgr[ri], $lgc[ci]), $rc_sample) seconds=1
t_c_ng = @be count(((ri, ci),) -> GO.crosses($ALG, $gor[ri], $go10[ci]), $rc_sample) seconds=1
t_c_lg = @be count(((ri, ci),) -> LG.crosses($lgr[ri], $lgc[ci]), $rc_sample) seconds=1

print_table("rivers x countries (10m, seeded sample of extent-hit pairs, per pair)",
    "predicate",
    ["RelateNG", "LibGEOS"],
    ["intersects" => per_op.([t_i_ng, t_i_lg], length(rc_sample)),
     "crosses"    => per_op.([t_c_ng, t_c_lg], length(rc_sample))])

# ## Workload 5: spherical variants
#
# NOTE: as of this file's creation, spherical `intersects` has a known
# containment bug (disjoint continent-scale pairs report true), so the
# spherical columns below measure the always-true early-exit path — see the
# caveat in the representative-output block.

const SALG = GO.RelateNG(; manifold = GO.Spherical())
sprep = GO.prepare(SALG, go_big)

t_s_unprep = @be count(p -> GO.intersects($SALG, $go_big, p), $(pts[1:5])) seconds=1
t_s_prep   = @be count(p -> GO.relate_predicate($sprep, GO.pred_intersects(), p), $(pts[1:20])) seconds=1
t_s_build  = @be GO.prepare($SALG, $go_big) seconds=1
t_s_pairs  = @be count(((i, j),) -> GO.intersects($SALG, $go110[i], $go110[j]), $(pair_sample[1:50])) seconds=1
t_p_pairs  = @be count(((i, j),) -> GO.intersects($ALG, $go110[i], $go110[j]), $(pair_sample[1:50])) seconds=1

print_table("spherical RelateNG vs planar (same real-data workloads, per op)",
    "workload",
    ["Spherical", "Planar"],
    ["point-in-area unprepared"  => [per_op(t_s_unprep, 5), pq_unprep],
     "point-in-area prepared"    => [per_op(t_s_prep, 20), pq_prep],
     "prepare build"             => [median_time(t_s_build), build],
     "pairwise intersects 110m"  => [per_op(t_s_pairs, 50), per_op(t_p_pairs, 50)]])

#=
Representative output (2026-07-14, Apple M4 Pro, macOS — Darwin 25.5.0;
Julia 1.12.6, GEOS 3.14.1, GeometryOps @ ba903d1ac — after the spherical
winding-independence fix and the engine type-erasure/precompile work;
`julia --project=docs benchmarks/relateng_realdata.jl`, ~95 s wall
excluding package precompilation):

Point-in-area target: Canada — 68193 vertices, 412 rings; point mix: 4672/10000 hits

point-in-area intersects (largest 10m country, seed-7 point mix, per query)
workload                       │         unprepared │     extent-stamped │           prepared │            LibGEOS │   LibGEOS prepared
───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Canada, per point              │           861.9 μs │            55.4 μs │             1.0 μs │            26.8 μs │           675.7 ns

GO.prepare build: 1.4 ms → amortized against unprepared after ~1.6 queries

110m: 177 countries (10654 vertices); 3000 sampled ordered pairs, 70 intersecting

pairwise intersects (110m countries, seeded pair sample, per pair)
workload                       │           RelateNG │             GO old │            LibGEOS
─────────────────────────────────────────────────────────────────────────────────────────────
3000 ordered pairs             │             2.8 μs │             3.3 μs │           116.5 ns

10m neighbor set: 180 extent-hit ordered pairs among top-24 countries, 18 touching

real border-sharing neighbor pairs (10m, 18 pairs, per pair)
predicate                      │           RelateNG │            LibGEOS
────────────────────────────────────────────────────────────────────────
touches                        │             1.7 ms │           617.9 μs
full relate                    │             1.6 ms │           620.0 μs

rivers x countries: 1454 rivers (256386 vertices), 5503 extent-hit pairs, 150 sampled

rivers x countries (10m, seeded sample of extent-hit pairs, per pair)
predicate                      │           RelateNG │            LibGEOS
────────────────────────────────────────────────────────────────────────
intersects                     │           305.4 μs │            35.0 μs
crosses                        │           995.4 μs │            79.8 μs

spherical RelateNG vs planar (same real-data workloads, per op)
workload                       │          Spherical │             Planar
────────────────────────────────────────────────────────────────────────
point-in-area unprepared       │           109.8 ms │           861.9 μs
point-in-area prepared         │            17.6 ms │             1.0 μs
prepare build                  │           101.6 ms │             1.4 ms
pairwise intersects 110m       │           169.8 μs │             2.9 μs

Reading notes:

- Point-in-area: the unprepared 33x gap vs LibGEOS is almost entirely the
  per-call extent-caching pass over all 68k vertices — pre-stamping extents
  (`GO.tuples(x; calc_extent = true)`) recovers 16x of it, and `GO.prepare`
  (indexed point-in-area locator) closes to ~2x of prepared GEOS. Prepare
  pays for itself after ~2 unprepared-equivalent queries.
- Pairwise 110m: RelateNG edges out the old GO processors and sits ~25x
  behind LibGEOS on these mostly-disjoint pairs (LibGEOS resolves most of
  them from the envelope alone).
- Real 10m borders (near-collinear point soups) force constant escalation
  to exact predicates: ~3x behind GEOS on `touches`/full `relate`.
- Rivers x countries: `crosses` is ~3x `intersects` — it requires
  self-noding, which rebuilds and re-traverses edge indexes per evaluation.
- Spherical: measured after the winding-independence fix (real answers, no
  always-true early exit). The ~100–17,000x gaps vs planar are dominated by
  kernel-point conversion — every query reconverts all ring vertices to
  `UnitSphericalPoint` — and by the absence of a spherical indexed
  point-in-area locator (prepared point queries re-scan converted rings).
  Both are addressed by docs/plans/2026-07-14-spherical-indexed-locator.md;
  re-record this table as those layers land.
=#
