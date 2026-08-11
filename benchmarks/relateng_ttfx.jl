# # RelateNG TTFX (fresh-process first-call) probe
#
#=
First-call latency ("time to first X") for RelateNG predicates. Each probe
instance spawns a *fresh* Julia process, loads GeometryOps, builds tiny
synthetic geometries, and times the first and second call of one predicate
on one geometry-type pair — so the first call is almost pure compile time
and the second call is steady state.

This is the tool that exposed the Julia 1.12 compile-time blowup: steady-
state runtime does *not* regress on 1.12, but first-call compilation of
polygon-target instances regresses 10-36x (e.g. `crosses(poly, poly)`
0.20 s on 1.11.9 vs 7.40 s on 1.12.6). A test suite touches many such
predicate x type-pair instances, which is what turned the package test
suite from ~8 min into ~48 min on 1.12 CI. The instance set here is the
worst tail of the campaign's full 3 x 7 x 7 sweep: the three boolean
predicates against polygonal/line targets.

Run with `julia --project=docs benchmarks/relateng_ttfx.jl`. The child
processes use the same project as the parent and, by default, the same
julia binary. To compare julia versions, point the children at another
binary via the `JULIA_EXE` environment variable or trailing ARGS (both
forms work with juliaup channel selectors):

    JULIA_EXE="julia +1.11" julia --project=docs benchmarks/relateng_ttfx.jl
    julia --project=docs benchmarks/relateng_ttfx.jl julia +1.12

The project must be instantiated for the child's julia version (for a
non-default version: copy `Manifest.toml` to `Manifest-v1.X.toml`, then
`julia +1.X --project=docs -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'`).
Children reuse the depot's precompile caches — on a cold cache the first
instance additionally pays package precompilation, so rerun for clean
numbers. Representative output for both julia channels is recorded in the
comment block at the bottom of this file; no CI gating.
=#

using Printf

const JULIA_CMD =
    !isempty(ARGS)           ? Cmd(String.(ARGS)) :
    haskey(ENV, "JULIA_EXE") ? Cmd(String.(split(ENV["JULIA_EXE"]))) :
                               Cmd([joinpath(Sys.BINDIR, "julia")])
const PROJECT = Base.active_project()

# The child program: time package load, then the first and second call of
# `ARGS = (predicate, type-of-a, type-of-b)`. Geometries are the tiny
# overlapping shapes from the profiling campaign's compile probe — compile
# time depends on types, not sizes. The call is wrapped in try/catch exactly
# like the campaign probe (some predicate/type combos may throw; the compile
# cost is what is being measured).
const CHILD_CODE = raw"""
const t0 = time_ns()
import GeometryOps as GO
import GeoInterface as GI
const t_load = (time_ns() - t0) / 1e9
p1 = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
p2 = GI.Polygon([[(2.0, 2.0), (5.0, 2.0), (5.0, 5.0), (2.0, 5.0), (2.0, 2.0)]])
geoms = Dict(
    "poly"  => p1,
    "mpoly" => GI.MultiPolygon([p1, p2]),
    "line"  => GI.LineString([(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)]),
)
preds = Dict("intersects" => GO.intersects, "touches" => GO.touches, "crosses" => GO.crosses)
f, a, b = preds[ARGS[1]], geoms[ARGS[2]], geoms[ARGS[3]]
alg = GO.RelateNG()
t1 = @elapsed try f(alg, a, b) catch end
t2 = @elapsed try f(alg, a, b) catch end
println("TTFX_RESULT ", VERSION, " ", t_load, " ", t1, " ", t2)
"""

function probe(pred, an, bn)
    cmd = `$JULIA_CMD --startup-file=no --project=$PROJECT -e $CHILD_CODE $pred $an $bn`
    buf = IOBuffer()
    ok = success(pipeline(cmd; stdout = buf, stderr = buf))
    out = String(take!(buf))
    m = match(r"TTFX_RESULT (\S+) (\S+) (\S+) (\S+)", out)
    (ok && m !== nothing) || error("child process failed for $pred($an, $bn):\n$out")
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

const PREDS = ("intersects", "touches", "crosses")
const PAIRS = (("poly", "poly"), ("mpoly", "poly"), ("line", "poly"))

results = Pair{String, NamedTuple}[]
for pred in PREDS, (an, bn) in PAIRS
    push!(results, "$pred($an, $bn)" => probe(pred, an, bn))
end

println("child: julia $(last(results).second.version) (`$(join(JULIA_CMD.exec, ' '))`)")
println("project: $PROJECT")
println()
printstyled("fresh-process first call (compile) vs second call (steady state)";
    color = :green, bold = true)
println()
@printf("%-24s", "instance")
foreach(c -> @printf(" │ %18s", c), ["package load", "first call", "second call"])
println()
println("─"^(24 + 21 * 3))
for (label, r) in results
    @printf("%-24s", label)
    foreach(t -> @printf(" │ %18s", prettytime(t)), [r.t_load, r.t_first, r.t_second])
    println()
end
println()

#=
Representative output (2026-07-14, Apple M4 Pro, macOS — Darwin 25.5.0;
GeometryOps @ ba903d1ac — after the engine type-erasure + PrecompileTools
workload; both julia channels installed via juliaup; warm precompile caches).

`julia --project=docs benchmarks/relateng_ttfx.jl`:

child: julia 1.12.6 (`.../julia-1.12.6+0.aarch64.apple.darwin14/.../bin/julia`)

fresh-process first call (compile) vs second call (steady state)
instance                 │       package load │         first call │        second call
───────────────────────────────────────────────────────────────────────────────────────
intersects(poly, poly)   │           631.6 ms │           322.5 μs │            17.1 μs
intersects(mpoly, poly)  │           374.4 ms │           243.3 ms │            25.7 μs
intersects(line, poly)   │           363.1 ms │           195.0 ms │            10.9 μs
touches(poly, poly)      │           374.9 ms │           436.6 μs │            22.0 μs
touches(mpoly, poly)     │           369.1 ms │           256.9 ms │            28.0 μs
touches(line, poly)      │           368.5 ms │           201.8 ms │            19.0 μs
crosses(poly, poly)      │           375.7 ms │            29.0 μs │             3.8 μs
crosses(mpoly, poly)     │           369.7 ms │           249.5 ms │            26.2 μs
crosses(line, poly)      │           382.0 ms │           204.9 ms │            26.2 μs

`JULIA_EXE="julia +1.11" julia --project=docs benchmarks/relateng_ttfx.jl`:

child: julia 1.11.9 (`julia +1.11`)

fresh-process first call (compile) vs second call (steady state)
instance                 │       package load │         first call │        second call
───────────────────────────────────────────────────────────────────────────────────────
intersects(poly, poly)   │           396.4 ms │           311.4 μs │             8.6 μs
intersects(mpoly, poly)  │           404.1 ms │            95.6 ms │            32.5 μs
intersects(line, poly)   │           408.3 ms │           106.2 ms │            21.1 μs
touches(poly, poly)      │           413.4 ms │             9.1 ms │            27.8 μs
touches(mpoly, poly)     │           407.7 ms │           140.8 ms │            53.6 μs
touches(line, poly)      │           413.8 ms │           127.6 ms │            40.3 μs
crosses(poly, poly)      │           420.4 ms │           105.7 μs │             6.4 μs
crosses(mpoly, poly)     │           413.6 ms │            72.8 ms │            28.5 μs
crosses(line, poly)      │           448.0 ms │           125.3 ms │            43.0 μs

Reading notes:

- Summed first-call time: ~1.35 s on 1.12.6 vs ~0.68 s on 1.11.9. Instances
  covered by the PrecompileTools workload (the `(poly, poly)` rows and
  `crosses(poly, poly)` in particular) resolve in the μs range; the rest pay
  ~100–260 ms of residual per-type-pair inference.
- History: at 099acacd1, before the engine core was type-erased (it was
  re-inferred per input geometry-type pair) and before the package had a
  PrecompileTools workload, this probe read 66.5 s summed first-call on
  1.12.6 vs 21.5 s on 1.11.9 — up to 8.6x per instance
  (`crosses(poly, poly)`: 9.61 s vs 1.12 s), with package load and
  steady-state equivalent across versions. That first-call pathology is what
  made CI's RelateNG testsets take 48 min on Julia 1.12 vs 7.65 min on 1.11.
  This probe is the regression gate for it: if the 1.12 column grows back
  into the seconds, the engine has re-grown geometry-type-parameterized
  internals.
=#
