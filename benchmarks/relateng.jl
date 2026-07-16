# # RelateNG benchmarks
#
#=
Benchmarks for the RelateNG DE-9IM engine (Task 28 of the relateng plan):
RelateNG vs the old GO per-pair predicate processors vs LibGEOS (plain and
prepared — GEOS >= 3.13 runs RelateNG natively), across polygon sizes, for
two workload shapes:

- `intersects`: early-exit-friendly — the engine can stop at the first
  interior/interior interaction it proves.
- full `relate`: no early exit — the complete DE-9IM matrix is computed.

Run with `julia --project=docs benchmarks/relateng.jl`. Prints a comparison
table; no CI gating (representative output is recorded in the comment block
at the bottom of this file).
=#

import GeometryOps as GO,
    GeoInterface as GI,
    LibGEOS as LG
using Chairmarks
using Printf
using Random

include(joinpath(@__DIR__, "..", "test", "data", "polygon_generation.jl"))

# As in `benchmarks.jl`: give each package its native geometry.
lg_and_go(geometry) = (GI.convert(LG, geometry), GO.tuples(geometry))

# A valid (LibGEOS-checked) random polygon with `nverts` vertices centered at
# `(x, y)`. The generator does not guarantee non-self-intersecting rings, so
# draw until valid (low spikiness makes rejection rare).
function valid_random_poly(x, y, nverts, rng)
    while true
        poly = GI.Polygon(generate_random_poly(x, y, nverts, 2.0, 0.3, 0.1, rng))
        LG.isValid(GI.convert(LG, poly)) && return poly
    end
end

const LG_CTX = LG.get_global_context()
# LibGEOS has no high-level `relate`; the generated wrapper returns the
# DE-9IM string itself (cf. test/methods/relateng/fuzz.jl).
lg_relate(la, lb) = LG.GEOSRelate_r(LG_CTX, la, lb)

prettytime(s) =
    s < 1e-6 ? @sprintf("%8.1f ns", s * 1e9) :
    s < 1e-3 ? @sprintf("%8.1f μs", s * 1e6) :
    s < 1.0  ? @sprintf("%8.1f ms", s * 1e3) :
               @sprintf("%8.2f s ", s)

median_time(trial) = Chairmarks.median(trial).time

function print_table(title, colnames, rows)
    printstyled(title; color = :green, bold = true)
    println()
    @printf("%8s", "nverts")
    foreach(c -> @printf(" │ %18s", c), colnames)
    println()
    println("─"^(8 + 21 * length(colnames)))
    for (n, times) in rows
        @printf("%8d", n)
        foreach(t -> @printf(" │ %18s", prettytime(t)), times)
        println()
    end
    println()
end

const ALG = GO.RelateNG()
const NVERTS = [2^4, 2^6, 2^8, 2^10, 2^12]

intersects_rows = Vector{Pair{Int, Vector{Float64}}}()
relate_rows = Vector{Pair{Int, Vector{Float64}}}()

for nverts in NVERTS
    rng = Xoshiro(42)
    # Two overlapping random polygons (centers 2.0 apart, radius 2.0), so
    # `intersects` is true and the full matrix is the overlaps pattern.
    a = valid_random_poly(0.0, 0.0, nverts, rng)
    b = valid_random_poly(2.0, 0.0, nverts, rng)
    lg_a, go_a = lg_and_go(a)
    lg_b, go_b = lg_and_go(b)
    prep_go = GO.prepare(ALG, go_a)
    prep_lg = LG.prepareGeom(lg_a)

    #-- Workload 1: `intersects` (early-exit-friendly)
    t_old      = @be GO.intersects($go_a, $go_b) seconds=1
    t_ng       = @be GO.intersects($ALG, $go_a, $go_b) seconds=1
    t_ng_prep  = @be GO.relate_predicate($prep_go, GO.pred_intersects(), $go_b) seconds=1
    t_lg       = @be LG.intersects($lg_a, $lg_b) seconds=1
    t_lg_prep  = @be LG.intersects($prep_lg, $lg_b) seconds=1
    push!(intersects_rows, nverts => median_time.([t_old, t_ng, t_ng_prep, t_lg, t_lg_prep]))

    #-- Workload 2: full `relate` (no early exit; old GO has no `relate`)
    t_ng       = @be GO.relate($ALG, $go_a, $go_b) seconds=1
    t_ng_prep  = @be GO.relate($prep_go, $go_b) seconds=1
    t_lg       = @be lg_relate($lg_a, $lg_b) seconds=1
    push!(relate_rows, nverts => median_time.([t_ng, t_ng_prep, t_lg]))
end

print_table("intersects (early exit, overlapping random polygons)",
    ["GO old", "RelateNG", "RelateNG prepared", "LibGEOS", "LibGEOS prepared"],
    intersects_rows)

print_table("relate (full DE-9IM matrix, overlapping random polygons)",
    ["RelateNG", "RelateNG prepared", "LibGEOS"],
    relate_rows)

#=
Representative output (2026-06-11, Apple M4 Pro, macOS — Darwin 25.5.0;
Julia 1.12.6, GEOS 3.14.1; `julia --project=docs benchmarks/relateng.jl`):

intersects (early exit, overlapping random polygons)
  nverts │             GO old │           RelateNG │  RelateNG prepared │            LibGEOS │   LibGEOS prepared
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
      16 │             1.7 μs │            11.6 μs │            11.0 μs │             2.1 μs │           827.4 ns
      64 │            23.0 μs │           137.8 μs │            70.4 μs │             4.7 μs │             2.1 μs
     256 │           328.4 μs │           562.7 μs │           282.6 μs │            13.8 μs │             7.5 μs
    1024 │             5.2 ms │            38.7 μs │            26.2 μs │             2.0 μs │           106.4 ns
    4096 │            79.8 ms │           155.5 μs │           105.1 μs │             6.7 μs │           469.1 ns

relate (full DE-9IM matrix, overlapping random polygons)
  nverts │           RelateNG │  RelateNG prepared │            LibGEOS
───────────────────────────────────────────────────────────────────────
      16 │            38.4 μs │            37.5 μs │             3.6 μs
      64 │           158.7 μs │            90.2 μs │             5.7 μs
     256 │           681.1 μs │           399.0 μs │            23.7 μs
    1024 │             3.7 ms │             2.6 ms │           201.1 μs
    4096 │            36.3 ms │            37.1 ms │             3.5 ms

Reading notes:

- The `intersects` columns are non-monotonic by design of the workload, not
  noise: above the accelerator threshold (1024+ vertices here) the STRtree
  edge index kicks in AND the early exit fires on the first proven
  interior/interior interaction, so the 1024/4096 rows are *cheaper* than
  256 (which takes the below-threshold nested loop over ~65k segment pairs
  to the first hit). GEOS short-circuits the same way throughout, and its
  prepared `intersects` resolves these cases from the point-in-area index
  alone (sub-μs).
- Old GO `intersects` scales quadratically here (5.2 ms / 80 ms at
  1024/4096) — RelateNG overtakes it between 256 and 1024 vertices and is
  ~500x faster at 4096.
- Full `relate` has no early exit, so every engine pays for the complete
  node analysis; RelateNG is within ~10x of native GEOS RelateNG at the
  large end. A flat profile of the 4096-vertex A/A `relate` case shows the
  time split between exact segment-pair classification (~25%, mostly
  AdaptivePredicates orientation — productive work, negligible call-boundary
  overhead), per-call segment-extent-table construction for the STRtree
  (~28%), and node analysis (Dict-backed `evaluate_nodes!` is ~8%). No
  rational-fallback (`rk_nodes_coincide`) frames appear at all — no
  F1-level hotspot.
=#
