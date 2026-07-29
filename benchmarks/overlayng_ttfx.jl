# # OverlayNG TTFX (fresh-process first-call) probe
#
#=
First-call latency ("time to first X") for the four OverlayNG overlay
operations, built exactly like `benchmarks/relateng_ttfx.jl`: each probe
instance spawns a *fresh* Julia process, loads GeometryOps, builds tiny
synthetic geometries, and times the first and second call of one operation on
one (manifold, geometry-type-pair) instance — so the first call is almost pure
compile time and the second call is steady state.

This is the measurement that decides what belongs in the package's
PrecompileTools workload (`src/precompile.jl`). Precompilation is a trade:
every workload instance buys first-call latency with package build time and
image size, so an instance is only worth adding if this probe says it is
expensive *and* it is on a path real users hit. The recorded numbers at the
bottom of this file are the before/after pair for the OverlayNG workload block.

Run with `julia --project=docs benchmarks/overlayng_ttfx.jl`. The child
processes use the same project as the parent and, by default, the same julia
binary; `JULIA_EXE` or trailing ARGS point them elsewhere, exactly as in
`relateng_ttfx.jl`:

    JULIA_EXE="julia +1.11" julia --project=docs benchmarks/overlayng_ttfx.jl
    julia --project=docs benchmarks/overlayng_ttfx.jl julia +1.12

Children reuse the depot's precompile caches — on a cold cache the first
instance additionally pays package precompilation, so rerun for clean numbers.
No CI gating.
=#

using Printf

const JULIA_CMD =
    !isempty(ARGS)           ? Cmd(String.(ARGS)) :
    haskey(ENV, "JULIA_EXE") ? Cmd(String.(split(ENV["JULIA_EXE"]))) :
                               Cmd([joinpath(Sys.BINDIR, "julia")])
const PROJECT = Base.active_project()

# The child program: time package load, then the first and second call of
# `ARGS = (op, type-of-a, type-of-b, manifold)`. The geometries are tiny —
# compile time depends on types, not sizes — but they are *valid* (multi
# components disjoint) and they genuinely interact, so the call runs the whole
# pipeline rather than exiting through an early empty-result path.
const CHILD_CODE = raw"""
const t0 = time_ns()
import GeometryOps as GO
import GeoInterface as GI
const t_load = (time_ns() - t0) / 1e9
p1 = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
p2 = GI.Polygon([[(2.0, 2.0), (5.0, 2.0), (5.0, 5.0), (2.0, 5.0), (2.0, 2.0)]])
p3 = GI.Polygon([[(6.0, 0.0), (8.0, 0.0), (8.0, 2.0), (6.0, 2.0), (6.0, 0.0)]])
geoms = Dict(
    "poly"   => p1,
    "mpoly"  => GI.MultiPolygon([p1, p3]),
    "line"   => GI.LineString([(-1.0, 1.0), (1.5, 1.5), (4.0, 1.0)]),
    "mline"  => GI.MultiLineString([[(-1.0, 1.0), (4.0, 1.0)], [(1.0, -1.0), (1.0, 4.0)]]),
    "point"  => GI.Point((1.0, 1.0)),
    "mpoint" => GI.MultiPoint([(1.0, 1.0), (9.0, 9.0)]),
)
ops = Dict("intersection" => GO.intersection, "union" => GO.union,
           "difference" => GO.difference, "symdifference" => GO.symdifference)
manifolds = Dict("planar" => GO.Planar(), "spherical" => GO.Spherical())
f, a, b = ops[ARGS[1]], geoms[ARGS[2]], geoms[ARGS[3]]
alg = GO.OverlayNG(manifolds[ARGS[4]])
t1 = @elapsed try f(alg, a, b) catch end
t2 = @elapsed try f(alg, a, b) catch end
println("TTFX_RESULT ", VERSION, " ", t_load, " ", t1, " ", t2)
"""

function probe(op, an, bn, mn)
    cmd = `$JULIA_CMD --startup-file=no --project=$PROJECT -e $CHILD_CODE $op $an $bn $mn`
    buf = IOBuffer()
    ok = success(pipeline(cmd; stdout = buf, stderr = buf))
    out = String(take!(buf))
    m = match(r"TTFX_RESULT (\S+) (\S+) (\S+) (\S+)", out)
    (ok && m !== nothing) || error("child process failed for $op($an, $bn) on $mn:\n$out")
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

const OPS = ("intersection", "union", "difference", "symdifference")
# The four input shapes a first call realistically lands on, on both manifolds
# (the manifold is a type parameter all the way down, so its instances are
# disjoint from the other's). `poly x poly` and `mpoly x mpoly` run the
# arrangement, `line x poly` adds the line builder, `point x poly` takes the
# separate mixed-points path that never reaches the arrangement at all.
const PAIRS = (("poly", "poly"), ("mpoly", "mpoly"), ("line", "poly"), ("point", "poly"))

results = Pair{String, NamedTuple}[]
for op in OPS, (an, bn) in PAIRS
    push!(results, "$op($an, $bn)" => probe(op, an, bn, "planar"))
end
for op in OPS, (an, bn) in PAIRS
    push!(results, "$op($an, $bn) sph" => probe(op, an, bn, "spherical"))
end

println("child: julia $(last(results).second.version) (`$(join(JULIA_CMD.exec, ' '))`)")
println("project: $PROJECT")
println()
printstyled("fresh-process first call (compile) vs second call (steady state)";
    color = :green, bold = true)
println()
@printf("%-32s", "instance")
foreach(c -> @printf(" │ %18s", c), ["package load", "first call", "second call"])
println()
println("─"^(32 + 21 * 3))
for (label, r) in results
    @printf("%-32s", label)
    foreach(t -> @printf(" │ %18s", prettytime(t)), [r.t_load, r.t_first, r.t_second])
    println()
end
println()
@printf("summed first call: %s over %d instances (planar %s, spherical %s)\n\n",
    strip(prettytime(sum(r.t_first for (_, r) in results))), length(results),
    strip(prettytime(sum(r.t_first for (l, r) in results if !endswith(l, "sph")))),
    strip(prettytime(sum(r.t_first for (l, r) in results if endswith(l, "sph")))))

#=
Representative output (2026-07-28, Apple M4 Pro, macOS — Darwin 25.5.0;
Julia 1.12.6, GeometryOps @ the phase-3 stack tip; warm precompile caches;
`julia --project=docs benchmarks/overlayng_ttfx.jl`, ~2 min per table).

## Before — no OverlayNG block in `src/precompile.jl`

fresh-process first call (compile) vs second call (steady state)
instance                         │       package load │         first call │        second call
───────────────────────────────────────────────────────────────────────────────────────────────
intersection(poly, poly)         │           609.8 ms │            1.10 s  │            11.8 μs
intersection(mpoly, mpoly)       │           232.4 ms │            1.21 s  │            53.4 μs
intersection(line, poly)         │           231.5 ms │            1.11 s  │            32.1 μs
intersection(point, poly)        │           234.5 ms │           423.1 ms │             7.2 μs
union(poly, poly)                │           231.7 ms │            1.10 s  │            13.1 μs
union(mpoly, mpoly)              │           236.6 ms │            1.14 s  │            20.2 μs
union(line, poly)                │           235.5 ms │            1.11 s  │            34.5 μs
union(point, poly)               │           232.6 ms │           545.7 ms │            45.0 μs
difference(poly, poly)           │           231.8 ms │            1.10 s  │            13.0 μs
difference(mpoly, mpoly)         │           236.9 ms │            1.13 s  │            20.2 μs
difference(line, poly)           │           235.3 ms │            1.11 s  │            21.2 μs
difference(point, poly)          │           246.4 ms │           499.3 ms │             6.5 μs
symdifference(poly, poly)        │           250.3 ms │            1.16 s  │            20.8 μs
symdifference(mpoly, mpoly)      │           247.8 ms │            1.18 s  │            28.4 μs
symdifference(line, poly)        │           249.5 ms │            1.15 s  │            41.5 μs
symdifference(point, poly)       │           245.8 ms │           556.5 ms │            22.2 μs
intersection(poly, poly) sph     │           243.0 ms │            2.10 s  │            71.5 μs
intersection(mpoly, mpoly) sph   │           264.5 ms │            2.02 s  │           116.1 μs
intersection(line, poly) sph     │           231.3 ms │            2.00 s  │           327.7 μs
intersection(point, poly) sph    │           239.5 ms │            1.15 s  │            50.8 μs
union(poly, poly) sph            │           233.7 ms │            1.98 s  │           120.0 μs
union(mpoly, mpoly) sph          │           230.8 ms │            1.99 s  │           112.1 μs
union(line, poly) sph            │           232.6 ms │            2.00 s  │           312.9 μs
union(point, poly) sph           │           229.4 ms │            1.28 s  │            72.2 μs
difference(poly, poly) sph       │           235.7 ms │            2.03 s  │           187.7 μs
difference(mpoly, mpoly) sph     │           231.0 ms │            2.07 s  │           289.2 μs
difference(line, poly) sph       │           232.7 ms │            2.00 s  │           319.9 μs
difference(point, poly) sph      │           230.5 ms │            1.17 s  │            38.3 μs
symdifference(poly, poly) sph    │           236.3 ms │            2.12 s  │           217.5 μs
symdifference(mpoly, mpoly) sph  │           255.7 ms │            2.09 s  │           309.4 μs
symdifference(line, poly) sph    │           234.3 ms │            2.01 s  │           299.3 μs
symdifference(point, poly) sph   │           240.9 ms │            1.37 s  │           101.8 μs

summed first call: 45.01 s over 32 instances (planar 15.61 s, spherical 29.39 s)

## After — with the OverlayNG block in `src/precompile.jl`

fresh-process first call (compile) vs second call (steady state)
instance                         │       package load │         first call │        second call
───────────────────────────────────────────────────────────────────────────────────────────────
intersection(poly, poly)         │           690.2 ms │           584.7 μs │            11.2 μs
intersection(mpoly, mpoly)       │           287.4 ms │            18.4 ms │            23.5 μs
intersection(line, poly)         │           293.3 ms │           241.6 μs │            20.9 μs
intersection(point, poly)        │           286.7 ms │           106.7 μs │             8.4 μs
union(poly, poly)                │           367.4 ms │             2.0 ms │            12.0 μs
union(mpoly, mpoly)              │           306.5 ms │            16.3 ms │            22.7 μs
union(line, poly)                │           291.1 ms │             2.2 ms │            44.8 μs
union(point, poly)               │           289.1 ms │            47.6 ms │            10.0 μs
difference(poly, poly)           │           285.3 ms │             2.3 ms │            18.7 μs
difference(mpoly, mpoly)         │           318.7 ms │            15.7 ms │            20.1 μs
difference(line, poly)           │           283.6 ms │             2.2 ms │            36.7 μs
difference(point, poly)          │           284.2 ms │             5.8 ms │             4.2 μs
symdifference(poly, poly)        │           287.6 ms │             1.9 ms │             8.5 μs
symdifference(mpoly, mpoly)      │           274.8 ms │            15.2 ms │            17.6 μs
symdifference(line, poly)        │           271.6 ms │             2.0 ms │            24.2 μs
symdifference(point, poly)       │           279.2 ms │            46.7 ms │            12.2 μs
intersection(poly, poly) sph     │           274.1 ms │             2.7 ms │            61.0 μs
intersection(mpoly, mpoly) sph   │           276.3 ms │            18.0 ms │           112.8 μs
intersection(line, poly) sph     │           275.2 ms │             2.7 ms │           277.0 μs
intersection(point, poly) sph    │           277.4 ms │             1.2 ms │            39.8 μs
union(poly, poly) sph            │           274.3 ms │             4.3 ms │            58.5 μs
union(mpoly, mpoly) sph          │           277.1 ms │            17.8 ms │           111.8 μs
union(line, poly) sph            │           273.1 ms │             4.7 ms │           266.6 μs
union(point, poly) sph           │           275.1 ms │            48.8 ms │            53.7 μs
difference(poly, poly) sph       │           277.1 ms │            13.8 ms │           158.1 μs
difference(mpoly, mpoly) sph     │           277.8 ms │            35.7 ms │           227.7 μs
difference(line, poly) sph       │           278.4 ms │             4.7 ms │           256.6 μs
difference(point, poly) sph      │           278.8 ms │             7.3 ms │            34.6 μs
symdifference(poly, poly) sph    │           272.7 ms │            14.0 ms │           141.8 μs
symdifference(mpoly, mpoly) sph  │           277.3 ms │            35.4 ms │           232.1 μs
symdifference(line, poly) sph    │           278.0 ms │             4.6 ms │           274.8 μs
symdifference(point, poly) sph   │           276.3 ms │            48.8 ms │            55.2 μs

summed first call: 443.9 ms over 32 instances (planar 179.3 ms, spherical 264.5 ms)

## The trade

Measured on the same machine, same session, `Pkg.precompile("GeometryOps")`
after invalidating `src/precompile.jl` (median of 3), and the pkgimage
`.dylib` beside the `.ji` in the depot:

                        precompile      pkgimage      summed first call (32)
    before                 14.59 s       23.99 MB                    45.01 s
    after                  18.98 s       32.84 MB                     0.44 s
    delta                  +4.39 s      +8.85 MB                    -44.6 s

Package load also grows with the image, from ~235 ms to ~277 ms (+42 ms) —
paid by every user of the package, including those who never call an overlay.

## What is in, what is out, and why

The workload is four calls per manifold: `intersection` on poly x poly,
mpoly x mpoly, line x poly, point x poly. The marginal ledger behind that
shape was measured by building each candidate as its own workload variant.
Run-to-run precompile *time* is only good to about ±1 s here, so the marginal
rows below are quoted in pkgimage bytes — which are exactly reproducible — and
in TTFX seconds; only the total precompile delta above is quoted in seconds.

    variant                                    pkgimage    Δ vs baseline
    baseline (no OverlayNG block)              23.99 MB               —
    + planar poly x poly                       25.99 MB          +2.00 MB
    + planar mpoly / line / point              28.45 MB          +4.46 MB
    + spherical poly x poly                    31.66 MB          +7.67 MB
    + spherical mpoly / line / point           32.84 MB          +8.85 MB
    (the same set without either point call)   30.07 MB          +6.08 MB

- **The other three ops: OUT.** `_overlay_ng` takes the op as a *value*
  (`op::_OverlayOpCode`), so one call caches the driver for all four. With
  only `intersection` precompiled, first-call `union`/`difference`/
  `symdifference` on poly x poly is ~2 ms rather than ~1.1 s; precompiling
  them too would buy ~6 ms per manifold. Not worth an instance.
- **The four input shapes, planar: IN.** +4.46 MB, and the planar half of the
  matrix goes from 15.61 s to 0.18 s. Note that 2.00 MB of that is the single
  poly x poly call, which alone drops the *whole* planar half except the point
  rows to a few ms: the engine core (`NodedArrangement{P}`, `OverlayGraph`, the
  builders) is typed on the kernel point type only, exactly like RelateNG's, so
  one call caches it for every input geometry type and the extra shapes add
  only their ingest and builder layers.
- **Spherical: IN.** +4.39 MB, and the spherical half goes from 29.39 s to
  0.26 s. The manifold is a type parameter all the way down, so *nothing* is
  shared with the planar instances: without this the first spherical overlay in
  a session costs 2 s. Spherical overlay is the reason this port exists, so
  this is squarely a path real users hit.
- **point x area: IN, on both manifolds.** +2.77 MB of the total (+2.07 planar,
  +0.70 spherical), removing 1.7 s (planar) + 4.9 s (spherical) of first-call
  latency. It is the most expensive megabyte in the block, because
  `_overlay_mixed_points` is a second pipeline that never touches the
  arrangement and therefore shares nothing with the rows above — which is
  exactly why leaving it out leaves `intersection(points, polygon)`, an
  operation the `OverlayNG` docstring advertises, at half a second.
- **mline / mpoint inputs: OUT.** They differ from `line` / `point` only in
  the ingest layer, which `mpoly` already exercises for the multi case; the
  residual is in the tens of ms.
- **`exact = False()`: OUT.** It doubles the engine instances (`exact` is a
  type parameter) for a configuration the `OverlayNG` docstring tells users
  not to select unless they are measuring the cost of exactness.

## Reading notes

- The residual after precompilation is ~18 ms for mpoly x mpoly and ~36-49 ms
  for the spherical multi/point rows — the per-input-type ingest layer, which
  scales with the number of type combinations rather than with engine size.
  Chasing it would mean precompiling a combinatorial matrix; the returns stop
  here.
- The first row of each table shows an inflated package load (~600-690 ms vs
  ~235-290 ms): that is the OS page cache warming on the pkgimage, not a real
  first-instance cost. Ignore it and read the steady rows.
- Steady-state (third column) is unchanged by precompilation, as it must be:
  ~10-50 us planar and ~40-330 us spherical on these tiny inputs. The spherical
  factor of ~5x is the kernel-point conversion and the exact spherical
  predicates, and it is also what `benchmarks/overlayng.jl` measures at scale.
=#
