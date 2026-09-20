# # Buffer TTFB (fresh-process first-call) probe
#
#=
First-call latency ("time to first buffer") for the native planar
[`buffer`](@ref), built exactly like `benchmarks/overlayng_ttfx.jl`: each probe
instance spawns a *fresh* Julia process, loads GeometryOps, builds one tiny
synthetic geometry, and times the first and second `buffer` call on it — so the
first call is almost pure compile time and the second is steady state.

This is the measurement that decides what belongs in the package's
PrecompileTools workload (`src/precompile.jl`). Precompilation is a trade: every
workload instance buys first-call latency with package build time and image
size, so an instance is worth adding only when this probe says it is expensive
*and* it is on a path real users hit. The recorded numbers at the bottom of this
file are the before/after pair for the buffer workload block.

Run with `julia --project=docs benchmarks/buffer_ttfb.jl`; `JULIA_EXE` or
trailing ARGS point the children at another binary, as in
`benchmarks/overlayng_ttfx.jl`. Children reuse the depot's precompile caches, so
rerun on a cold cache for clean numbers. No CI gating.
=#

using Printf

const JULIA_CMD =
    !isempty(ARGS)           ? Cmd(String.(ARGS)) :
    haskey(ENV, "JULIA_EXE") ? Cmd(String.(split(ENV["JULIA_EXE"]))) :
                               Cmd([joinpath(Sys.BINDIR, "julia")])
const PROJECT = Base.active_project()

# The child program: time package load, then the first and second `buffer` of
# `ARGS = (shape, distance)`. Geometries are tiny — compile time depends on
# types, not sizes — but each one reaches a different code path: `poly` and
# `holed` run the whole offset-curve-plus-winding-overlay pipeline, `line` and
# `point` change the generator's ingest, `mpoly` and `collection` add the
# flattening walk, `vector` adds `apply`'s container recursion, and `empty` (a
# negative distance on a line) returns before the engine is entered at all.
const CHILD_CODE = raw"""
const t0 = time_ns()
import GeometryOps as GO
import GeoInterface as GI
const t_load = (time_ns() - t0) / 1e9
p = GI.Polygon([[(0.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0), (0.0, 0.0)]])
h = GI.Polygon([[(0.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0), (0.0, 0.0)],
                [(2.5, 0.8), (3.5, 0.8), (3.5, 3.2), (2.5, 3.2), (2.5, 0.8)]])
q = GI.Polygon([[(8.0, 0.0), (10.0, 0.0), (10.0, 2.0), (8.0, 0.0)]])
l = GI.LineString([(0.0, 0.0), (3.0, 0.5), (4.0, 3.0), (1.5, 2.0)])
pt = GI.Point((1.0, 1.0))
geoms = Dict(
    "poly"       => p,
    "holed"      => h,
    "line"       => l,
    "point"      => pt,
    "mpoly"      => GI.MultiPolygon([p, q]),
    "collection" => GI.GeometryCollection([pt, l, p]),
    "vector"     => [p, l],
    "empty"      => l,
)
g = geoms[ARGS[1]]
d = parse(Float64, ARGS[2])
t1 = @elapsed GO.buffer(g, d)
t2 = @elapsed GO.buffer(g, d)
println("TTFB_RESULT ", VERSION, " ", t_load, " ", t1, " ", t2)
"""

function probe(shape, d)
    cmd = `$JULIA_CMD --startup-file=no --project=$PROJECT -e $CHILD_CODE $shape $d`
    buf = IOBuffer()
    ok = success(pipeline(cmd; stdout = buf, stderr = buf))
    out = String(take!(buf))
    m = match(r"TTFB_RESULT (\S+) (\S+) (\S+) (\S+)", out)
    (ok && m !== nothing) || error("child process failed for buffer($shape, $d):\n$out")
    return (; version = m[1],
              t_load = parse(Float64, m[2]),
              t_first = parse(Float64, m[3]),
              t_second = parse(Float64, m[4]))
end

prettytime(s) =
    s < 1e-6 ? @sprintf("%8.1f ns", s * 1e9) :
    s < 1e-3 ? @sprintf("%8.1f μs", s * 1e6) :
    s < 1.0  ? @sprintf("%8.1f ms", s * 1e3) :
               @sprintf("%8.2f s ", s)

const INSTANCES = (("poly", 0.5), ("holed", -0.5), ("line", 0.5), ("point", 1.0),
                   ("mpoly", 0.3), ("collection", 0.3), ("vector", 0.3), ("empty", -1.0))

results = [("buffer($shape, $d)" => probe(shape, d)) for (shape, d) in INSTANCES]

println("child: julia $(last(results).second.version) (`$(join(JULIA_CMD.exec, ' '))`)")
println("project: $PROJECT")
println()
printstyled("fresh-process first call (compile) vs second call (steady state)";
    color = :green, bold = true)
println()
@printf("%-28s", "instance")
foreach(c -> @printf(" │ %18s", c), ["package load", "first call", "second call"])
println()
println("─"^(28 + 21 * 3))
for (label, r) in results
    @printf("%-28s", label)
    foreach(t -> @printf(" │ %18s", prettytime(t)), [r.t_load, r.t_first, r.t_second])
    println()
end
println()
@printf("summed first call: %s over %d instances\n\n",
    strip(prettytime(sum(r.t_first for (_, r) in results))), length(results))

#=
Representative output (2026-09-17, 16-core x86-64 Linux; Julia 1.13.0,
GeometryOps @ the native-buffer branch tip; warm precompile caches;
`julia --project=test benchmarks/buffer_ttfb.jl`, ~15 s per table).

## Before the buffer workload block

    instance                     │       package load │         first call │        second call
    buffer(poly, 0.5)            │           447.8 ms │           977.2 ms │            50.6 μs
    buffer(holed, -0.5)          │           457.3 ms │           976.3 ms │            71.1 μs
    buffer(line, 0.5)            │           467.0 ms │           941.9 ms │            88.4 μs
    buffer(point, 1.0)           │           445.9 ms │           893.1 ms │            54.9 μs
    buffer(mpoly, 0.3)           │           479.6 ms │           974.2 ms │            80.2 μs
    buffer(collection, 0.3)      │           445.4 ms │            1.00 s  │           135.9 μs
    buffer(vector, 0.3)          │           438.5 ms │            1.11 s  │           112.8 μs
    buffer(empty, -1.0)          │           460.9 ms │           942.0 ms │            11.1 μs

    summed first call: 7.82 s over 8 instances

## After

    buffer(poly, 0.5)            │           449.7 ms │           139.7 μs │            43.3 μs
    buffer(holed, -0.5)          │           447.4 ms │           164.9 μs │            61.7 μs
    buffer(line, 0.5)            │           446.6 ms │             4.0 ms │            62.5 μs
    buffer(point, 1.0)           │           457.7 ms │             3.7 ms │            36.6 μs
    buffer(mpoly, 0.3)           │           459.4 ms │             3.9 ms │            67.9 μs
    buffer(collection, 0.3)      │           445.2 ms │           101.0 ms │           165.3 μs
    buffer(vector, 0.3)          │           449.3 ms │            95.3 ms │           178.7 μs
    buffer(empty, -1.0)          │           439.1 ms │             3.7 ms │             9.1 μs

    summed first call: 211.8 ms over 8 instances

Within a single session the eight shapes share nearly everything, so the same
calls in one process cost 1.43 s before the block and 0.25 s after — read the
fresh-process tables for what one instance is worth and the single-session
figure for what a script pays.

## Ledger

Package precompile 35.7 s before, 36.9 s after, against a 33.3-36.9 s
run-to-run spread — the block's build cost sits below the noise floor.
Pkgimage 48.61 MB -> 49.77 MB (+1.16 MB), which is ~5 ms of extra package load.

| Workload instance | Verdict | Buys |
|:--|:--|:--|
| `buffer(poly, 0.5)`, through the keyword entry point | IN | 977 ms -> 140 μs, and with it the whole engine core for every other input type. ~1.1 MB, the bulk of the block. |
| `buffer(alg, line, 0.5)` | IN | 942 ms -> 4.0 ms, for the generator's line-and-cap ingest. Nearly free after the areal row. |
| `buffer(alg, point, 0.5)` | IN | 893 ms -> 3.7 ms, for `_offset_circle`. Nearly free. |
| `buffer(alg, mpoly, 0.5)` | IN | 974 ms -> 3.9 ms for +0.06 MB, covering the multi-part flattening walk. |
| A vector of geometries, and a GeometryCollection | OUT | ~1.05 s -> ~98 ms each comes free with the rows above; the residual is `apply`'s container recursion, which every operation shares and none of them caches. |
| A second `quadsegs` | OUT | `quadsegs` is a field, not a type parameter, so varying it costs no instances. |
| `exact = False()` | OUT | A type parameter, so it doubles the engine instances — for a configuration the docstring tells users not to select unless they are measuring the cost of exactness. |
| `Spherical()` | N/A | No native spherical buffer exists. |

## Reading notes

- Steady state is unchanged by precompilation, as it must be: 9 μs for an empty
  result and 37-180 μs on these tiny inputs.
- The residual 3.7-4.0 ms on the line, point, multipolygon and empty rows is the
  generator's per-input-type ingest, which scales with the number of type
  combinations rather than with engine size. Chasing it means precompiling a
  combinatorial matrix; the returns stop here.
- `buffer(empty, -1.0)` — a negative distance on a line — returns an empty
  `MultiPolygon` without building an arrangement, so its 942 ms first call is
  the entry point and the generator's input cleaning alone. That it drops to
  3.7 ms without being in the workload is the clearest single reading of how
  much the areal row caches.
=#
