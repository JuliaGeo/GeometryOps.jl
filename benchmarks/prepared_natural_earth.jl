# # Prepared-geometry point-in-polygon benchmark on Natural Earth data
#
# Measures whether GeometryOps' prepared geometry (per-ring edge trees + cached
# extent) speeds up `contains(polygon, point)` on real-world country borders,
# and by how much as a function of polygon size and Natural Earth zoom level.
#
# Usage (one scale per process):
#
#     julia --project=test benchmarks/prepared_natural_earth.jl 110
#     julia --project=test benchmarks/prepared_natural_earth.jl 50
#     julia --project=test benchmarks/prepared_natural_earth.jl 10
#
# Method (see the header comments in each section for the fairness rationale):
#   * the plain baseline uses `GO.tuples`-converted polygons (the strongest
#     plain layout); the prepared variants are built from the RAW GeoJSON
#     geometries, since `prepare` materializes into native layout itself —
#     that conversion cost is part of the prepare timing;
#   * each polygon gets a fixed 15x15 grid of points over its bbox, shared by
#     every backend;
#   * a correctness gate asserts the inside-count agrees across all backends
#     before any timing happens;
#   * timings are min-of-3 `@elapsed` over the whole workload, after one warmup.

import GeometryOps as GO
import GeoInterface as GI
using NaturalEarth
using Printf, Statistics

const Unsorted = GO.FlexibleRTrees.Unsorted
const HPR      = GO.FlexibleRTrees.HPR

# 15x15 = 225 points per polygon; a mix of inside/outside is realistic.
const GRID_N = 15
const REPS   = 3

# Backends under test.  The first is the plain (unprepared) baseline; the rest
# build a `Prepared` wrapper whose edge trees use the named tree backend.
# `NaturalIndex` is the default preparation, so `GO.prepare(p)` selects it.
prepare_naturalindex(p) = GO.prepare(p)
prepare_unsorted(p)     = GO.prepare(p; preps = (GO.EdgeTree(Unsorted()),))
prepare_hpr(p)          = GO.prepare(p; preps = (GO.EdgeTree(HPR()),))

const PREP_BACKENDS = [
    ("NaturalIndex", prepare_naturalindex),
    ("Unsorted",     prepare_unsorted),
    ("HPR",          prepare_hpr),
]

# Vertex-count buckets for size-stratified reporting.
const BUCKET_LABELS = ["<100", "100-1k", "1k-10k", ">10k"]
bucket_of(v) = v < 100 ? 1 : v < 1_000 ? 2 : v < 10_000 ? 3 : 4

# Total exterior + hole vertices of a polygon.
function polygon_vertices(poly)
    v = GI.npoint(GI.getexterior(poly))
    for hole in GI.gethole(poly)
        v += GI.npoint(hole)
    end
    return v
end

# A fixed n*n grid of Float64 points over the polygon's bounding box.
function bbox_grid(poly, n)
    ext = GI.extent(poly)
    xmin, xmax = Float64.(ext.X)
    ymin, ymax = Float64.(ext.Y)
    xs = range(xmin, xmax; length = n)
    ys = range(ymin, ymax; length = n)
    return [(x, y) for x in xs for y in ys]
end

# The timed workload: every point of every grid tested against its polygon.
# Returns the total inside-count so the compiler can't elide the work.
function run_queries(geoms, grids)
    total = 0
    @inbounds for i in eachindex(geoms)
        g = geoms[i]
        grid = grids[i]
        for pt in grid
            total += GO.contains(g, pt)::Bool
        end
    end
    return total
end

# min-of-REPS elapsed seconds, after one warmup call.
function timed(f)
    f()
    return minimum(@elapsed(f()) for _ in 1:REPS)
end

build_all(build, polys) = map(build, polys)

function main(scale)
    println("=" ^ 72)
    println("Prepared point-in-polygon benchmark - Natural Earth admin_0_countries @ $(scale)m")
    println("=" ^ 72)

    # --- Load and convert (once) -------------------------------------------
    fc = naturalearth("admin_0_countries", scale)
    raw = collect(GO.flatten(GI.PolygonTrait, fc))
    # Skip degenerate polygons (need >= 4 points to form a closed ring).
    keep = findall(p -> GI.npoint(GI.getexterior(p)) >= 4, raw)
    n_skipped = length(raw) - length(keep)
    raw = raw[keep]
    # Plain baseline: tuple-converted (the strongest plain layout).  Prepared
    # variants are built from `raw` — materialization is `prepare`'s job.
    polys = GO.tuples.(raw)

    vcounts = polygon_vertices.(polys)
    grids   = [bbox_grid(p, GRID_N) for p in polys]
    nqueries = sum(length, grids)

    println()
    @printf("polygons: %d (skipped %d degenerate)\n", length(polys), n_skipped)
    @printf("vertices: total=%d  median=%d  max=%d\n",
            sum(vcounts), round(Int, median(vcounts)), maximum(vcounts))
    @printf("grid: %dx%d = %d points/polygon  ->  %d total point-in-polygon queries\n",
            GRID_N, GRID_N, GRID_N^2, nqueries)
    bucket_counts = [count(==(b), bucket_of.(vcounts)) for b in 1:4]
    println("size buckets ", BUCKET_LABELS, " -> ", bucket_counts, " polygons")

    # --- Build prepared variants + time preparation ------------------------
    println("\nPreparing geometries (timing preparation itself)...")
    prepared = Dict{String,Vector}()
    prep_time = Dict{String,Float64}()
    for (name, build) in PREP_BACKENDS
        prep_time[name] = timed(() -> build_all(build, raw))
        prepared[name] = build_all(build, raw)
    end

    # --- Correctness gate --------------------------------------------------
    # Every backend must agree on the inside-count for every polygon before we
    # trust (or time) any of them.  A mismatch is a bug, reported and fatal.
    println("Correctness gate: comparing inside-counts across all backends...")
    mismatches = 0
    for i in eachindex(polys)
        grid = grids[i]
        base = count(pt -> GO.contains(polys[i], pt), grid)
        for (name, _) in PREP_BACKENDS
            c = count(pt -> GO.contains(prepared[name][i], pt), grid)
            if c != base
                mismatches += 1
                @printf("  MISMATCH polygon #%d (%d verts): plain=%d %s=%d\n",
                        i, vcounts[i], base, name, c)
                if mismatches <= 5
                    ext = GI.extent(polys[i])
                    println("    extent = ", ext)
                end
            end
        end
    end
    if mismatches > 0
        error("Correctness gate FAILED: $mismatches polygon/backend mismatches. " *
              "This is a bug report - see the MISMATCH lines above. Aborting benchmark.")
    end
    println("  OK - all backends agree on all $(length(polys)) polygons.")

    # --- Whole-workload timing --------------------------------------------
    println("\nTiming query workload (min of $REPS, whole dataset)...")
    plain_time = timed(() -> run_queries(polys, grids))
    query_time = Dict{String,Float64}("plain" => plain_time)
    for (name, _) in PREP_BACKENDS
        query_time[name] = timed(() -> run_queries(prepared[name], grids))
    end

    # --- Report: per-scale backend table -----------------------------------
    println("\n### Backend comparison @ $(scale)m")
    println("| backend      | total query (s) | per-query (ns) | speedup vs plain | prepare (s) |")
    println("|--------------|-----------------|----------------|------------------|-------------|")
    perq(t) = t / nqueries * 1e9
    @printf("| %-12s | %15.4f | %14.1f | %16s | %11s |\n",
            "plain", plain_time, perq(plain_time), "1.00x", "-")
    for (name, _) in PREP_BACKENDS
        t = query_time[name]
        @printf("| %-12s | %15.4f | %14.1f | %15.2fx | %11.4f |\n",
                name, t, perq(t), plain_time / t, prep_time[name])
    end

    # --- Report: speedup by polygon size (plain vs NaturalIndex) -----------
    # Per-bucket timing, so small polygons (where preparation can't pay off)
    # are separated from large ones (where it should).
    println("\n### Speedup by polygon size - NaturalIndex vs plain @ $(scale)m")
    println("| size bucket | polygons | queries | plain (s) | NaturalIndex (s) | speedup |")
    println("|-------------|----------|---------|-----------|------------------|---------|")
    for b in 1:4
        idx = findall(==(b), bucket_of.(vcounts))
        isempty(idx) && continue
        bpolys = polys[idx]
        bpreps = prepared["NaturalIndex"][idx]
        bgrids = grids[idx]
        bq = sum(length, bgrids)
        tp = timed(() -> run_queries(bpolys, bgrids))
        tn = timed(() -> run_queries(bpreps, bgrids))
        @printf("| %-11s | %8d | %7d | %9.4f | %16.4f | %6.2fx |\n",
                BUCKET_LABELS[b], length(idx), bq, tp, tn, tp / tn)
    end

    # --- Break-even analysis ----------------------------------------------
    # How many queries before preparation pays for itself (NaturalIndex):
    #   plain cost for N queries  = N * t_plain
    #   prepared cost for N       = T_prepare + N * t_prep
    #   break-even N*             = T_prepare / (t_plain - t_prep)
    t_plain = plain_time / nqueries
    t_prep  = query_time["NaturalIndex"] / nqueries
    Tprep   = prep_time["NaturalIndex"]
    println("\n### Break-even (NaturalIndex) @ $(scale)m")
    if t_plain > t_prep
        Nstar = Tprep / (t_plain - t_prep)
        @printf("  prepare cost = %.4f s total (%.2f us/polygon)\n",
                Tprep, Tprep / length(polys) * 1e6)
        @printf("  per-query saving = %.1f ns  ->  break-even at %.0f total queries",
                (t_plain - t_prep) * 1e9, Nstar)
        @printf(" (~%.1f queries/polygon)\n", Nstar / length(polys))
    else
        println("  prepared is not faster on average at this scale; no break-even.")
    end

    # --- Anomalies: polygons where prepared was slower ---------------------
    # Time each backend per polygon-size bucket already covers the trend; here
    # we flag the whole-dataset case explicitly.
    slower = [name for (name, _) in PREP_BACKENDS if query_time[name] > plain_time]
    if isempty(slower)
        println("\nNo backend was slower than plain on the whole $(scale)m workload.")
    else
        println("\nBackends slower than plain on the whole $(scale)m workload: ", join(slower, ", "))
    end
    println()
end

# --- Entry point -----------------------------------------------------------
scale = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 110
scale in (110, 50, 10) || error("scale must be one of 110, 50, 10 (got $scale)")
main(scale)
