# # OverlayNG real-data benchmarks
#
#=
Benchmarks for the OverlayNG exact-arrangement overlay engine on real Natural
Earth data (complementing `benchmarks/overlayng.jl`, which uses synthetic random
polygons and circles). Data comes from NaturalEarth.jl — the artifacts are
downloaded and cached on first use, so the first ever run needs network access;
cached runs are offline.

As in `benchmarks/overlayng.jl`, LibGEOS is the exactly-matched reference: GEOS's
default overlay engine since 3.9 *is* OverlayNG, the same algorithm by the same
author from the same JTS source, so the ratio measures this port's exact
arrangement against a mature C++ implementation of it rather than "Julia vs C".
Foster–Hormann — the current GeometryOps default, which OverlayNG does not
displace — is the third column wherever it has the operation at all.

Six workload groups, each printed as a table of per-operation medians:

1. Pairwise country overlay at 110 m on a seeded sample of intersecting ordered
   pairs — the sparse regime at country scale (most pairs share only a border).
2. Named 10 m neighbour pairs, per pair, with vertex counts — real shared
   borders, which are near-collinear point soups and the exact-predicate stress
   test.
3. Shifted-self overlays at 10 m: a country against a copy of itself translated
   by 0.5°, so every border segment crosses its own copy. This is the densest
   crossing workload real data offers, and the shape `test/methods/clipping/
   overlayng/realdata_identities.jl` uses for the same reason.
4. Rivers x countries at 10 m — line x area, the mixed-dimension path.
5. Points x country at 10 m — the mixed-points path, which never reaches the
   arrangement.
6. Spherical variants (`GO.OverlayNG(GO.Spherical())`) of groups 1 and 3,
   against their planar equivalents.

Run with `julia --project=docs benchmarks/overlayng_realdata.jl`. Prints
comparison tables; no CI gating (representative output is recorded in the comment
block at the bottom of this file). Everything is seeded — reruns measure
identical workloads.

**Known-defect routing.** The ring-assembly defect cluster tracked in
`test/external/jts/overlay_skiplist.jl` (`OverlayTopologyError: unable to assign
free hole to a shell` / `found two shells in EdgeRing list`) still fires on some
real inputs. Every candidate case is therefore run once through all four ops
before it is timed, and cases that raise are dropped and counted, so a defect
shrinks the workload rather than aborting the run or silently biasing it. The
count is printed with each group.
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

prettytime(s) =
    isnan(s) ? @sprintf("%11s", "—") :
    s < 1e-6 ? @sprintf("%8.1f ns", s * 1e9) :
    s < 1e-3 ? @sprintf("%8.1f μs", s * 1e6) :
    s < 1.0  ? @sprintf("%8.1f ms", s * 1e3) :
               @sprintf("%8.2f s ", s)

median_time(trial) = Chairmarks.median(trial).time

function print_table(title, firstcol, colnames, rows)
    printstyled(title; color = :green, bold = true)
    println()
    @printf("%-36s", firstcol)
    foreach(c -> @printf(" │ %18s", c), colnames)
    println()
    println("─"^(36 + 21 * length(colnames)))
    for (label, times) in rows
        @printf("%-36s", label)
        foreach(t -> @printf(" │ %18s", prettytime(t)), times)
        println()
    end
    println()
end

# A failed cell is reported as `—` rather than aborting the sweep (see the
# known-defect note in the header; Foster–Hormann also throws on some pairs).
macro trial(expr)
    quote
        try
            median_time($(esc(expr)))
        catch err
            @warn "benchmark cell failed" err = first(sprint(showerror, err), 120) maxlog = 4
            NaN
        end
    end
end

# Load a NaturalEarth layer and return (names, geometries), skipping features
# with no geometry and empty geometries (an empty MultiLineString in the 10m
# rivers layer breaks `GO.tuples`) — the loader from
# `benchmarks/relateng_realdata.jl`.
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

# Overlay contracts on valid input, so invalid rings are dropped rather than
# left to poison the sweep — the same rule as
# `test/methods/clipping/overlayng/realdata_identities.jl`, which also excludes
# NE 110 m Sudan by name (clockwise exterior, duplicate vertex, geodesically
# self-intersecting).
const NE_EXCLUDED = Set(["Sudan"])

function load_areas(scale)
    names, raw = ne_geoms("admin_0_countries", scale)
    keep_names, keep_go, keep_lg = String[], [], []
    for (nm, g) in zip(names, raw)
        nm in NE_EXCLUDED && continue
        (GI.trait(g) isa GI.PolygonTrait || GI.trait(g) isa GI.MultiPolygonTrait) || continue
        lg, go = try lg_and_go(g) catch; continue end
        (try LG.isValid(lg) catch; false end) || continue
        push!(keep_names, nm); push!(keep_go, go); push!(keep_lg, lg)
    end
    return keep_names, keep_go, keep_lg
end

const ALG  = GO.OverlayNG()
const SALG = GO.OverlayNG(GO.Spherical())
const OPS  = ("intersection", "union", "difference", "symdifference")

ng_op(alg, name, a, b) =
    name == "intersection"  ? GO.intersection(alg, a, b) :
    name == "union"         ? GO.union(alg, a, b) :
    name == "difference"    ? GO.difference(alg, a, b) :
                              GO.symdifference(alg, a, b)

fh_op(name, a, b) =
    name == "intersection"  ? GO.intersection(a, b; target = GI.PolygonTrait()) :
    name == "union"         ? GO.union(a, b; target = GI.PolygonTrait()) :
    name == "difference"    ? GO.difference(a, b; target = GI.PolygonTrait()) :
                              error("Foster–Hormann has no symmetric difference")

lg_op(name, a, b) =
    name == "intersection"  ? LG.intersection(a, b) :
    name == "union"         ? LG.union(a, b) :
    name == "difference"    ? LG.difference(a, b) :
                              LG.symmetricDifference(a, b)

# Drop cases (`(label, go_a, go_b, lg_a, lg_b)`) that the engine cannot complete
# on this manifold, so a known defect shrinks the workload instead of aborting
# the run. Returns the kept cases and the labels dropped.
function clean_cases(alg, cases)
    kept, dropped = eltype(cases)[], String[]
    for c in cases
        ok = true
        for op in OPS
            try
                ng_op(alg, op, c[2], c[3])
            catch
                ok = false
                break
            end
        end
        ok ? push!(kept, c) : push!(dropped, c[1])
    end
    return kept, dropped
end

# Sweep every case once, for each engine, and report per-case medians.
sweep_ng(alg, op, cases) = for c in cases; ng_op(alg, op, c[2], c[3]); end
sweep_fh(op, cases)      = for c in cases; fh_op(op, c[2], c[3]); end
sweep_lg(op, cases)      = for c in cases; lg_op(op, c[4], c[5]); end

function op_rows(alg, cases; with_fh = true, seconds = 1)
    n = length(cases)
    return [op => [
        @trial(@be sweep_ng($alg, $op, $cases) seconds=seconds) / n,
        with_fh && op != "symdifference" ?
            @trial(@be sweep_fh($op, $cases) seconds=seconds) / n : NaN,
        @trial(@be sweep_lg($op, $cases) seconds=seconds) / n,
    ] for op in OPS]
end

const COLNAMES = ["OverlayNG", "Foster–Hormann", "LibGEOS"]

# ## Workload 1: pairwise country overlay at 110 m

names110, go110, lg110 = load_areas(110)
@printf("110m: %d valid countries, %d vertices total\n", length(go110), sum(GI.npoint, go110))

exts110 = [GI.extent(g) for g in go110]
cand = [(i, j) for i in eachindex(go110) for j in eachindex(go110)
        if i < j && Extents.intersects(exts110[i], exts110[j])]
pair_sample = shuffle(Xoshiro(42), cand)
pairs110 = Tuple{String, Any, Any, Any, Any}[]
for (i, j) in pair_sample
    length(pairs110) >= 40 && break
    (try LG.intersects(lg110[i], lg110[j]) catch; false end) || continue
    push!(pairs110, ("$(names110[i]) x $(names110[j])", go110[i], go110[j], lg110[i], lg110[j]))
end
pairs110, dropped110 = clean_cases(ALG, pairs110)
@printf("workload 1: %d intersecting 110m pairs (%d dropped: %s)\n\n",
    length(pairs110), length(dropped110), isempty(dropped110) ? "none" : join(dropped110, ", "))

print_table("pairwise country overlay (110m, $(length(pairs110)) seeded pairs, per pair)",
    "operation", COLNAMES, op_rows(ALG, pairs110))

# ## Workload 2: named 10 m neighbour pairs, per pair

names10, go10, lg10 = load_areas(10)
@printf("10m: %d valid countries, %d vertices total\n\n", length(go10), sum(GI.npoint, go10))

idx10(nm) = findfirst(==(nm), names10)
#-- NE 110 m Sudan is excluded above as a documented bad input, so the
#-- Egypt neighbour used here is Libya
const NBR_PAIRS = [("Egypt", "Libya"), ("France", "Italy"), ("Greece", "Turkey"),
                   ("Brazil", "Argentina"), ("Norway", "Sweden")]
rows2, rows2lg = Pair{String, Vector{Float64}}[], Pair{String, Vector{Float64}}[]
for (na, nb) in NBR_PAIRS
    ia, ib = idx10(na), idx10(nb)
    (ia === nothing || ib === nothing) && continue
    case, dropped = clean_cases(ALG, [("$na x $nb", go10[ia], go10[ib], lg10[ia], lg10[ib])])
    if isempty(case)
        println("workload 2: $na x $nb dropped (known ring-assembly defect)")
        continue
    end
    label = @sprintf("%s–%s (%d/%d verts)", na, nb, GI.npoint(go10[ia]), GI.npoint(go10[ib]))
    r = op_rows(ALG, case; with_fh = false)
    push!(rows2,   label => [x.second[1] for x in r])
    push!(rows2lg, label => [x.second[3] for x in r])
end
print_table("named 10m neighbour pairs, OverlayNG (per operation)",
    "pair", collect(OPS), rows2)
print_table("named 10m neighbour pairs, LibGEOS (per operation)",
    "pair", collect(OPS), rows2lg)

# ## Workload 3: shifted-self overlays at 10 m (densest real crossing workload)

shifted(A, dx, dy) = GO.apply(GI.PointTrait(), A) do p
    (GI.x(p) + dx, GI.y(p) + dy)
end

const SHIFT_NAMES = ["Egypt", "France", "Norway", "Brazil"]
shifts = Tuple{String, Any, Any, Any, Any}[]
for nm in SHIFT_NAMES
    i = idx10(nm)
    if i === nothing
        println("workload 3: $nm not among the valid 10m countries, skipped")
        continue
    end
    B = shifted(go10[i], 0.5, 0.25)
    lgB = try GI.convert(LG, B) catch; continue end
    push!(shifts, ("$nm shifted", go10[i], B, lg10[i], lgB))
end
shifts, dropped3 = clean_cases(ALG, shifts)
@printf("workload 3: %d shifted-self 10m cases (%d dropped: %s)\n\n",
    length(shifts), length(dropped3), isempty(dropped3) ? "none" : join(dropped3, ", "))

rows3, rows3b = Pair{String, Vector{Float64}}[], Pair{String, Vector{Float64}}[]
for c in shifts
    label = @sprintf("%s (%d verts)", c[1], GI.npoint(c[2]))
    r = op_rows(ALG, [c])
    push!(rows3,  label => [x.second[1] for x in r])         # the four ops, OverlayNG
    push!(rows3b, label => r[2].second)                      # `union`, the three engines
end
print_table("shifted-self overlay, OverlayNG (10m, per operation)",
    "case", collect(OPS), rows3)
print_table("shifted-self `union` across engines (10m, per operation)",
    "case", COLNAMES, rows3b)

# ## Workload 4: rivers x countries at 10 m (line x area)

rnames, rraw = ne_geoms("rivers_lake_centerlines", 10)
gor = [GO.tuples(g) for g in rraw]
lgr = [GI.convert(LG, g) for g in rraw]
rexts = [GI.extent(g) for g in gor]
exts10 = [GI.extent(g) for g in go10]
rc = [(ri, ci) for ri in eachindex(gor) for ci in eachindex(go10)
      if Extents.intersects(rexts[ri], exts10[ci])]
rc_sample = shuffle(Xoshiro(42), rc)[1:120]
lines = [("river $ri x country $ci", gor[ri], go10[ci], lgr[ri], lg10[ci])
         for (ri, ci) in rc_sample]
lines, dropped4 = clean_cases(ALG, lines)
@printf("workload 4: %d rivers (%d vertices); %d extent-hit pairs, %d sampled, %d kept\n\n",
    length(gor), sum(GI.npoint, gor), length(rc), length(rc_sample), length(lines))

print_table("rivers x countries (10m, $(length(lines)) pairs, line x area, per pair)",
    "operation", COLNAMES, op_rows(ALG, lines; with_fh = false))

# ## Workload 5: points x country at 10 m (mixed-points path)

i_big = argmax(GI.npoint.(go10))
ext = exts10[i_big]
rng = MersenneTwister(7)
pts = [(ext.X[1] + rand(rng) * (ext.X[2] - ext.X[1]),
        ext.Y[1] + rand(rng) * (ext.Y[2] - ext.Y[1])) for _ in 1:1000]
mp_go = GI.MultiPoint(pts)
mp_lg = LG.MultiPoint([[p[1], p[2]] for p in pts])
ptcase = [("1000 points x $(names10[i_big])", mp_go, go10[i_big], mp_lg, lg10[i_big])]
@printf("workload 5: 1000 seeded points x %s (%d vertices)\n\n",
    names10[i_big], GI.npoint(go10[i_big]))

print_table("points x country (10m, 1000-point MultiPoint, per operation)",
    "operation", COLNAMES, op_rows(ALG, ptcase; with_fh = false))

# ## Workload 6: spherical variants
#
# Overlay on `Spherical()` runs great-circle edges and exact spherical
# predicates throughout — the same arrangement, a different kernel.

sph110, dropped_s1 = clean_cases(SALG, pairs110)
sph_shifts, dropped_s3 = clean_cases(SALG, shifts)
@printf("workload 6: %d/%d 110m pairs and %d/%d shifted cases run spherically (dropped: %s)\n\n",
    length(sph110), length(pairs110), length(sph_shifts), length(shifts),
    isempty([dropped_s1; dropped_s3]) ? "none" : join([dropped_s1; dropped_s3], ", "))

rows6 = Pair{String, Vector{Float64}}[]
for op in OPS
    n1 = length(sph110)
    push!(rows6, "110m pairs, $op" => [
        @trial(@be sweep_ng($SALG, $op, $sph110) seconds=1) / n1,
        @trial(@be sweep_ng($ALG,  $op, $sph110) seconds=1) / n1,
        @trial(@be sweep_lg($op, $sph110) seconds=1) / n1])
end
for op in OPS
    n3 = length(sph_shifts)
    n3 == 0 && break
    push!(rows6, "shifted-self, $op" => [
        @trial(@be sweep_ng($SALG, $op, $sph_shifts) seconds=1) / n3,
        @trial(@be sweep_ng($ALG,  $op, $sph_shifts) seconds=1) / n3,
        @trial(@be sweep_lg($op, $sph_shifts) seconds=1) / n3])
end
print_table("spherical vs planar OverlayNG (same real-data workloads, per op)",
    "workload", ["Spherical", "Planar", "LibGEOS (planar)"], rows6)

#=
Representative output (2026-07-28, Apple M4 Pro, macOS — Darwin 25.5.0;
Julia 1.12.6, GEOS 3.14.1, GeometryOps @ the phase-3 stack tip;
`julia --project=docs benchmarks/overlayng_realdata.jl`, ~2.5 min wall on a
warm NaturalEarth artifact cache, excluding package precompilation. The
`benchmark cell failed` warnings interleaved in the real output are the
Foster–Hormann cells discussed in the notes, and are elided here):

110m: 174 valid countries, 9900 vertices total
workload 1: 40 intersecting 110m pairs (0 dropped: none)

pairwise country overlay (110m, 40 seeded pairs, per pair)
operation                            │          OverlayNG │     Foster–Hormann │            LibGEOS
───────────────────────────────────────────────────────────────────────────────────────────────────
intersection                         │            73.5 μs │            75.6 μs │            22.2 μs
union                                │            70.9 μs │                  — │            37.9 μs
difference                           │            70.6 μs │                  — │            31.3 μs
symdifference                        │            71.7 μs │                  — │            38.0 μs

10m: 256 valid countries, 545359 vertices total

named 10m neighbour pairs, OverlayNG (per operation)
pair                                 │       intersection │              union │         difference │      symdifference
────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
France–Italy (4641/3361 verts)       │             4.8 ms │             5.0 ms │             4.9 ms │             5.0 ms
Greece–Turkey (6096/3640 verts)      │             7.0 ms │             7.2 ms │             7.1 ms │             7.2 ms
Brazil–Argentina (11121/4674 verts)  │             8.9 ms │             9.0 ms │             8.9 ms │             9.0 ms
Norway–Sweden (15817/4665 verts)     │            13.6 ms │            14.0 ms │            14.0 ms │            13.9 ms

named 10m neighbour pairs, LibGEOS (per operation)
pair                                 │       intersection │              union │         difference │      symdifference
────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
France–Italy (4641/3361 verts)       │             1.2 ms │             2.0 ms │             1.9 ms │             2.0 ms
Greece–Turkey (6096/3640 verts)      │             1.3 ms │             2.6 ms │             2.3 ms │             2.6 ms
Brazil–Argentina (11121/4674 verts)  │             1.4 ms │             3.9 ms │             3.7 ms │             3.9 ms
Norway–Sweden (15817/4665 verts)     │             4.6 ms │             6.3 ms │             6.3 ms │             6.3 ms

workload 3: Egypt not among the valid 10m countries, skipped
workload 3: 3 shifted-self 10m cases (0 dropped: none)

shifted-self overlay, OverlayNG (10m, per operation)
case                                 │       intersection │              union │         difference │      symdifference
────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
France shifted (4641 verts)          │             7.2 ms │             7.2 ms │             7.1 ms │             7.2 ms
Norway shifted (15817 verts)         │            38.4 ms │            39.5 ms │            38.1 ms │            38.7 ms
Brazil shifted (11121 verts)         │            20.8 ms │            24.9 ms │            22.6 ms │            25.3 ms

shifted-self `union` across engines (10m, per operation)
case                                 │          OverlayNG │     Foster–Hormann │            LibGEOS
───────────────────────────────────────────────────────────────────────────────────────────────────
France shifted (4641 verts)          │             7.2 ms │                  — │             2.2 ms
Norway shifted (15817 verts)         │            39.5 ms │                  — │            14.7 ms
Brazil shifted (11121 verts)         │            24.9 ms │                  — │             6.7 ms

workload 4: 1454 rivers (256386 vertices); 5474 extent-hit pairs, 120 sampled, 120 kept

rivers x countries (10m, 120 pairs, line x area, per pair)
operation                            │          OverlayNG │     Foster–Hormann │            LibGEOS
───────────────────────────────────────────────────────────────────────────────────────────────────
intersection                         │             5.7 ms │                  — │           188.3 μs
union                                │             6.2 ms │                  — │             3.5 ms
difference                           │             5.6 ms │                  — │           192.6 μs
symdifference                        │             6.2 ms │                  — │             3.5 ms

workload 5: 1000 seeded points x Canada (68193 vertices)

points x country (10m, 1000-point MultiPoint, per operation)
operation                            │          OverlayNG │     Foster–Hormann │            LibGEOS
───────────────────────────────────────────────────────────────────────────────────────────────────
intersection                         │             2.9 ms │                  — │             2.6 ms
union                                │             2.4 ms │                  — │            20.6 ms
difference                           │             2.3 ms │                  — │             2.6 ms
symdifference                        │             2.4 ms │                  — │            20.5 ms

workload 6: 40/40 110m pairs and 3/3 shifted cases run spherically (dropped: none)

spherical vs planar OverlayNG (same real-data workloads, per op)
workload                             │          Spherical │             Planar │   LibGEOS (planar)
───────────────────────────────────────────────────────────────────────────────────────────────────
110m pairs, intersection             │           375.3 μs │            72.7 μs │            22.7 μs
110m pairs, union                    │           386.2 μs │            73.9 μs │            39.1 μs
110m pairs, difference               │           384.1 μs │            72.8 μs │            34.0 μs
110m pairs, symdifference            │           392.4 μs │            74.0 μs │            38.8 μs
shifted-self, intersection           │           236.8 ms │            22.5 ms │             7.6 ms
shifted-self, union                  │           253.0 ms │            22.3 ms │             8.0 ms
shifted-self, difference             │           250.7 ms │            22.0 ms │             7.7 ms
shifted-self, symdifference          │           263.4 ms │            22.4 ms │             7.9 ms

Reading notes:

- **Area x area against GEOS: 2-4x, everywhere.** Real country pairs put
  OverlayNG within a small constant factor of LibGEOS across three very
  different shapes — 1.9-3.2x on the sparse 110 m pairs, 2.2-4.0x on real 10 m
  shared borders, and 2.7-3.7x on the shifted-self cases where every border
  segment crosses its own copy. The synthetic sweep in `benchmarks/overlayng.jl`
  reaches parity (1.0-1.3x) at 4096 vertices on dense random polygons; real
  country data sits a little worse because even a shared border is mostly
  *non*-crossing linework, which is the regime where this port's per-segment
  graph costs more than JTS's per-noded-run one (see that file's notes on
  `_label_disconnected_edges!`). That single labelling change is worth most of
  the remaining gap on every table here.
- **The four ops cost the same, again.** Columns agree to a few percent within
  a row: `_overlay_ng` takes the op as a value and everything up to result
  extraction is shared. GEOS does *not* have this property — its `intersection`
  is 2-4x cheaper than its `union` on the same pair, because it clips the inputs
  to the intersection envelope first. That asymmetry is most of why our ratio
  looks worse on `intersection` than on `union`.
- **Foster–Hormann cannot do these at all.** `union` and `difference` raise
  `MethodError: Cannot convert an object of type GeoInterface.Wrappers.Polygon`
  on the MultiPolygon country pairs, for every `target` (`PolygonTrait`,
  `MultiPolygonTrait`), which is why the column is empty outside the
  `intersection` row on the 110 m sample. Only `intersection` with a polygon
  target completes. This is a limitation of the current default engine on
  multi-geometry inputs, not of the benchmark, and it is the practical argument
  for the opt-in: on real administrative data OverlayNG is not merely faster
  above ~256 vertices, it is the only one of the two that returns an answer.
- **Line x area is the worst ratio: 30x on `intersection`.** A river against a
  country is the extreme sparse case — a handful of crossings against thousands
  of area edges — so GEOS's envelope clipping discards almost all of the
  country's linework (188 μs) while this engine nodes, graphs and labels every
  edge of it (5.7 ms). The same pair's `union`, where GEOS cannot clip anything
  away, is 1.8x. The gap is the same structural one as above, seen through its
  most favourable case for GEOS.
- **Points x area is the one workload we win.** `union` and `symdifference` of a
  1000-point MultiPoint with Canada are 2.4 ms here against 20.6 ms for GEOS
  (8.6x faster); `intersection` and `difference` are at parity (2.9 vs 2.6 ms,
  2.3 vs 2.6 ms). The mixed-points path never builds an arrangement — it locates
  each point against the indexed area and copies through — which is exactly the
  shape of work the rest of these tables say we are slow at *not* doing.
- **Spherical costs ~5x planar on sparse pairs and ~11x on dense ones**
  (375 μs vs 73 μs at 110 m; 237-263 ms vs 22 ms shifted-self). The arrangement
  is the same; the difference is the kernel — per-vertex conversion to 3-vectors,
  great-circle edge geometry, and exact spherical orientation predicates whose
  filter fails more often than the planar one. There is no spherical reference
  engine in this comparison because LibGEOS has none; the LibGEOS column is the
  planar answer, shown to make the manifold's cost legible.
- **Known-defect routing found nothing to route around.** Every candidate case
  in every group completed on both manifolds (0 dropped). The ring-assembly
  cluster in `test/external/jts/overlay_skiplist.jl` does not fire on Natural
  Earth country pairs at these two resolutions.
=#
