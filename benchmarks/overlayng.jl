# # OverlayNG benchmarks
#
#=
Benchmarks for the OverlayNG exact-arrangement overlay engine (phase 3 of the
overlay plan): OverlayNG vs the default Foster–Hormann clipping engine vs
LibGEOS, across polygon sizes, for the four overlay operations.

**LibGEOS is the exactly-matched reference here.** GEOS's default overlay engine
since 3.9 *is* OverlayNG — the same algorithm by the same author, ported from
the same JTS source. The GO/GEOS ratio in the last table is therefore not a
"Julia vs C" number in the usual loose sense: it is the price of this port's
design decisions (exact predicates on every uncertain filter, symbolic crossing
nodes, coordinates rounded once at emission) against a mature C++ implementation
of the same algorithm that instead uses a floating-point precision model and a
snapping retry ladder.

Two workload shapes, each over five polygon sizes:

- **random overlapping polygons** — two irregular random polygons whose centres
  are one radius apart, so their boundaries cross O(n) times. This is the dense
  regime: work is dominated by noding and by the arrangement.
- **circle vs shifted circle** — a regular circle against a copy of itself
  translated by 0.6 radii, which crosses exactly twice at any vertex count. This
  is the sparse regime: work is dominated by ingest and by copying the
  non-crossing arcs through to the result. It is the same shape `benchmarks.jl`
  uses for the Foster–Hormann suite, so the numbers are comparable to it.

Run with `julia --project=docs benchmarks/overlayng.jl`. Prints comparison
tables; no CI gating (representative output is recorded in the comment block at
the bottom of this file). Everything is seeded — reruns measure identical
workloads.
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
# draw until valid (low spikiness makes rejection rare). Overlay contracts on
# valid input, so this filter is a correctness requirement here, not hygiene.
function valid_random_poly(x, y, nverts, rng)
    while true
        poly = GI.Polygon(generate_random_poly(x, y, nverts, 2.0, 0.3, 0.1, rng))
        LG.isValid(GI.convert(LG, poly)) && return poly
    end
end

# A closed regular circle of `nverts` vertices, centred at `(x, y)`.
circle_poly(x, y, nverts) = GO.ClosedRing()(GI.Polygon([[
    (x + cos(θ), y + sin(θ)) for θ in LinRange(0, 2π, nverts)]]))

prettytime(s) =
    isnan(s)   ? @sprintf("%11s", "—") :
    s < 1e-6   ? @sprintf("%8.1f ns", s * 1e9) :
    s < 1e-3   ? @sprintf("%8.1f μs", s * 1e6) :
    s < 1.0    ? @sprintf("%8.1f ms", s * 1e3) :
                 @sprintf("%8.2f s ", s)

prettyratio(r) = isnan(r) ? @sprintf("%11s", "—") : @sprintf("%9.1fx", r)

median_time(trial) = Chairmarks.median(trial).time

# Every engine here can fail on some input: Foster–Hormann has no
# `symdifference` and throws on some pairs, and OverlayNG has a known
# ring-assembly defect cluster (`OverlayTopologyError`) still being fixed. A
# failed cell is reported as `—` rather than aborting the sweep.
macro trial(expr)
    quote
        try
            median_time($(esc(expr)))
        catch err
            @warn "benchmark cell failed" err = first(sprint(showerror, err), 120)
            NaN
        end
    end
end

function print_table(title, colnames, rows; fmt = prettytime)
    printstyled(title; color = :green, bold = true)
    println()
    @printf("%8s", "nverts")
    foreach(c -> @printf(" │ %18s", c), colnames)
    println()
    println("─"^(8 + 21 * length(colnames)))
    for (n, times) in rows
        @printf("%8d", n)
        foreach(t -> @printf(" │ %18s", fmt(t)), times)
        println()
    end
    println()
end

const ALG = GO.OverlayNG()
const NVERTS = [2^4, 2^6, 2^8, 2^10, 2^12]

# One size of one workload: OverlayNG, Foster–Hormann and LibGEOS on each of the
# four ops. Returns `(op => [ng, fh, lg])` rows. Foster–Hormann is the current
# default engine and takes a `target`; it has no symmetric difference at all,
# which is why that cell is structurally empty.
function measure(go_a, go_b, lg_a, lg_b)
    tgt = GI.PolygonTrait()
    return [
        "intersection" => [
            @trial(@be GO.intersection($ALG, $go_a, $go_b) seconds=1),
            @trial(@be GO.intersection($go_a, $go_b; target = $tgt) seconds=1),
            @trial(@be LG.intersection($lg_a, $lg_b) seconds=1)],
        "union" => [
            @trial(@be GO.union($ALG, $go_a, $go_b) seconds=1),
            @trial(@be GO.union($go_a, $go_b; target = $tgt) seconds=1),
            @trial(@be LG.union($lg_a, $lg_b) seconds=1)],
        "difference" => [
            @trial(@be GO.difference($ALG, $go_a, $go_b) seconds=1),
            @trial(@be GO.difference($go_a, $go_b; target = $tgt) seconds=1),
            @trial(@be LG.difference($lg_a, $lg_b) seconds=1)],
        "symdifference" => [
            @trial(@be GO.symdifference($ALG, $go_a, $go_b) seconds=1),
            NaN,
            @trial(@be LG.symmetricDifference($lg_a, $lg_b) seconds=1)],
    ]
end

const COLNAMES = ["OverlayNG", "Foster–Hormann", "LibGEOS"]

# `results[workload][op]` is a vector of `nverts => [ng, fh, lg]` rows.
results = Dict{String, Dict{String, Vector{Pair{Int, Vector{Float64}}}}}()
for w in ("random", "circle")
    results[w] = Dict(op => Pair{Int, Vector{Float64}}[]
                      for op in ("intersection", "union", "difference", "symdifference"))
end

for nverts in NVERTS
    #-- workload 1: two overlapping irregular random polygons (dense crossings)
    rng = Xoshiro(42)
    a = valid_random_poly(0.0, 0.0, nverts, rng)
    b = valid_random_poly(2.0, 0.0, nverts, rng)
    lg_a, go_a = lg_and_go(a)
    lg_b, go_b = lg_and_go(b)
    for (op, ts) in measure(go_a, go_b, lg_a, lg_b)
        push!(results["random"][op], nverts => ts)
    end

    #-- workload 2: a circle against a copy shifted by 0.6 radii (2 crossings)
    lg_c1, go_c1 = lg_and_go(circle_poly(0.0, 0.0, nverts))
    lg_c2, go_c2 = lg_and_go(circle_poly(0.6, 0.0, nverts))
    for (op, ts) in measure(go_c1, go_c2, lg_c1, lg_c2)
        push!(results["circle"][op], nverts => ts)
    end
end

for (w, wlabel) in (("random", "overlapping random polygons"),
                    ("circle", "circle vs circle shifted 0.6 radii"))
    for op in ("intersection", "union", "difference", "symdifference")
        print_table("$op ($wlabel)", COLNAMES, results[w][op])
    end
end

#-- the headline: what the exact arrangement costs against the same algorithm
#-- in C++ (GEOS >= 3.9 *is* OverlayNG)
for (w, wlabel) in (("random", "overlapping random polygons"),
                    ("circle", "circle vs circle shifted 0.6 radii"))
    rows = Pair{Int, Vector{Float64}}[]
    for (i, nverts) in enumerate(NVERTS)
        push!(rows, nverts => [results[w][op][i].second[1] / results[w][op][i].second[3]
                               for op in ("intersection", "union", "difference", "symdifference")])
    end
    print_table("OverlayNG ÷ LibGEOS ($wlabel)",
        ["intersection", "union", "difference", "symdifference"], rows; fmt = prettyratio)
end

#=
Representative output (2026-07-28, Apple M4 Pro, macOS — Darwin 25.5.0;
Julia 1.12.6, GEOS 3.14.1, GeometryOps @ the phase-3 stack tip;
`julia --project=docs benchmarks/overlayng.jl`, ~2 min wall):

intersection (overlapping random polygons)
  nverts │          OverlayNG │     Foster–Hormann │            LibGEOS
───────────────────────────────────────────────────────────────────────
      16 │            17.4 μs │             3.8 μs │             8.3 μs
      64 │            75.7 μs │            58.6 μs │            16.6 μs
     256 │           335.3 μs │           916.5 μs │            76.3 μs
    1024 │             1.9 ms │            14.8 ms │           875.9 μs
    4096 │            16.6 ms │           232.8 ms │            13.1 ms

union (overlapping random polygons)
  nverts │          OverlayNG │     Foster–Hormann │            LibGEOS
───────────────────────────────────────────────────────────────────────
      16 │            17.6 μs │             3.9 μs │             6.6 μs
      64 │            75.7 μs │            59.0 μs │            14.2 μs
     256 │           335.0 μs │           920.2 μs │            82.7 μs
    1024 │             2.0 ms │            14.9 ms │             1.3 ms
    4096 │            17.9 ms │           234.1 ms │            18.1 ms

difference (overlapping random polygons)
  nverts │          OverlayNG │     Foster–Hormann │            LibGEOS
───────────────────────────────────────────────────────────────────────
      16 │            17.1 μs │             3.9 μs │             7.4 μs
      64 │            74.9 μs │            58.7 μs │            15.8 μs
     256 │           331.5 μs │           932.0 μs │            82.9 μs
    1024 │             1.9 ms │            14.6 ms │             1.1 ms
    4096 │            16.3 ms │           232.0 ms │            14.8 ms

symdifference (overlapping random polygons)
  nverts │          OverlayNG │     Foster–Hormann │            LibGEOS
───────────────────────────────────────────────────────────────────────
      16 │            17.5 μs │                  — │             7.2 μs
      64 │            76.4 μs │                  — │            15.3 μs
     256 │           345.5 μs │                  — │            89.1 μs
    1024 │             1.9 ms │                  — │             1.2 ms
    4096 │            16.6 ms │                  — │            16.3 ms

intersection (circle vs circle shifted 0.6 radii)
  nverts │          OverlayNG │     Foster–Hormann │            LibGEOS
───────────────────────────────────────────────────────────────────────
      16 │            18.9 μs │             4.2 μs │             8.0 μs
      64 │            80.5 μs │            60.8 μs │            13.4 μs
     256 │           322.8 μs │           917.9 μs │            29.3 μs
    1024 │             1.4 ms │            14.5 ms │            87.2 μs
    4096 │             5.9 ms │           230.5 ms │           306.5 μs

union (circle vs circle shifted 0.6 radii)
  nverts │          OverlayNG │     Foster–Hormann │            LibGEOS
───────────────────────────────────────────────────────────────────────
      16 │            18.6 μs │             4.2 μs │             7.2 μs
      64 │            80.0 μs │            60.9 μs │             9.1 μs
     256 │           318.1 μs │           920.7 μs │            15.8 μs
    1024 │             1.4 ms │            14.5 ms │            38.1 μs
    4096 │             5.8 ms │           232.1 ms │           118.2 μs

difference (circle vs circle shifted 0.6 radii)
  nverts │          OverlayNG │     Foster–Hormann │            LibGEOS
───────────────────────────────────────────────────────────────────────
      16 │            19.2 μs │             4.3 μs │             8.1 μs
      64 │            79.2 μs │            61.1 μs │            11.2 μs
     256 │           317.6 μs │           914.5 μs │            22.9 μs
    1024 │             1.4 ms │            14.7 ms │            62.3 μs
    4096 │             5.8 ms │           231.1 ms │           211.7 μs

symdifference (circle vs circle shifted 0.6 radii)
  nverts │          OverlayNG │     Foster–Hormann │            LibGEOS
───────────────────────────────────────────────────────────────────────
      16 │            19.7 μs │                  — │             7.7 μs
      64 │            80.8 μs │                  — │            10.0 μs
     256 │           321.4 μs │                  — │            17.4 μs
    1024 │             1.4 ms │                  — │            43.2 μs
    4096 │             5.9 ms │                  — │           138.4 μs

OverlayNG ÷ LibGEOS (overlapping random polygons)
  nverts │       intersection │              union │         difference │      symdifference
────────────────────────────────────────────────────────────────────────────────────────────
      16 │               2.1x │               2.7x │               2.3x │               2.4x
      64 │               4.6x │               5.3x │               4.8x │               5.0x
     256 │               4.4x │               4.1x │               4.0x │               3.9x
    1024 │               2.2x │               1.6x │               1.8x │               1.6x
    4096 │               1.3x │               1.0x │               1.1x │               1.0x

OverlayNG ÷ LibGEOS (circle vs circle shifted 0.6 radii)
  nverts │       intersection │              union │         difference │      symdifference
────────────────────────────────────────────────────────────────────────────────────────────
      16 │               2.4x │               2.6x │               2.4x │               2.6x
      64 │               6.0x │               8.8x │               7.1x │               8.1x
     256 │              11.0x │              20.1x │              13.9x │              18.5x
    1024 │              15.8x │              36.0x │              21.7x │              32.4x
    4096 │              19.3x │              49.1x │              27.3x │              42.9x

Reading notes:

- **Dense regime: parity with GEOS.** On overlapping random polygons OverlayNG
  converges to 1.0-1.3x LibGEOS at 4096 vertices (16.6 vs 13.1 ms for
  intersection, 17.9 vs 18.1 ms for union, i.e. a hair *faster* on union). This
  is the headline result of the port: an exact arrangement in Julia costs
  essentially nothing against a mature C++ implementation of the same algorithm
  once there is real topology to compute. The small-size rows (2-5x) are
  per-call overhead, not scaling — see below.
- **Against the default engine.** OverlayNG overtakes Foster–Hormann between 64
  and 256 vertices and is 14x faster at 4096 (16.6 vs 232.8 ms), because
  Foster–Hormann's per-pair processors are quadratic in the vertex count while
  the arrangement is not. It is nonetheless ~4x *slower* below ~64 vertices, and
  that is the honest trade of the opt-in: the arrangement has a fixed setup cost
  (segment-string extraction, node table, edge index, half-edge graph) that a
  16-vertex clip does not amortize.
- **Sparse regime: 19-49x behind GEOS, and this is the one real gap.** Two
  circles crossing exactly twice give OverlayNG 5.8 ms at 4096 vertices against
  GEOS's 118 μs. It is not ingest: extracting both inputs to segment strings is
  23.8 μs of that 5.8 ms (0.4%). The split at n = 4096 (union) is

      segment strings   23.8 μs   |  arrangement (noding)  1.6 ms
      graph build      564.6 μs   |  labelling             3.0 ms
      extraction       188.4 μs   |  (LibGEOS total       121.2 μs)

  and inside labelling, 2.8 of the 3.0 ms is pass 5, `_label_disconnected_edges!`
  (JTS `labelDisconnectedEdges`) — the per-edge point-in-area location. The
  reason is structural: this port makes *every original vertex a node* (design
  §3 amendments 1–2 — half-edge directions are original-coordinate differences,
  which is what makes the angular sort exact), so two 4096-vertex circles
  produce 8194 nodes and 16392 half-edges, where JTS/GEOS produce one `Edge` per
  noded run between crossings, i.e. about four. Label propagation therefore
  reaches only the half-edges incident to the two real crossings (16376 of 16392
  are still unknown when pass 5 starts) and every remaining edge is located by
  an indexed PIP. The dense regime hides this because there the two edge counts
  converge. The fix is local to the labeller — propagate through degree-2 nodes,
  i.e. carry a label along a whole run rather than re-locating each segment of
  it — and is worth roughly the whole sparse-regime gap. It is *not* the
  `RingClipper` / `LineLimiter` input clipping that JTS uses and this design
  deliberately drops (that constructs coordinates, which the exactness invariant
  forbids, and it would not help here anyway: the overlap envelope of these two
  circles covers most of both).
- **The four ops cost the same.** Within a workload the four columns agree to a
  few percent, because `_overlay_ng` takes the op as a *value* and the whole
  pipeline up to result extraction is shared; only the marking predicate and the
  emitted ring set differ. Foster–Hormann has no symmetric difference at all,
  which is the structurally empty column.
- **Per-call floor.** At 16 vertices OverlayNG is ~17 μs against LibGEOS's ~7 μs
  and Foster–Hormann's ~4 μs. That is the setup cost of the arrangement, and it
  puts a floor under any workload that overlays many tiny geometries — the
  amortization story there is a prepared/reusable arrangement (S1 kept ingest
  separable for exactly this), not a faster pipeline.
=#
