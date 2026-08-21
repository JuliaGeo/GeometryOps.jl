# # Controlled OverlayNG scaling benchmarks
#
#=
This harness separates three sources of overlay cost:

- `vertices`: two convex polygons with two boundary crossings while their
  vertex counts grow;
- `crossings`: two 4,096-segment polygons while their boundary-crossing count
  grows from 2 to 1,024;
- `serrated`: a radial saw polygon against a rotated copy, so both input size
  and the number of local tooth intersections grow together.

Geometry construction, conversion, validation, and crossing counts happen
outside the timed expression. GeometryOps receives `GO.tuples` geometries and
LibGEOS receives `LibGEOS.Polygon`s.

Run every workload with:

```sh
julia --project=docs benchmarks/overlayng_scaling.jl
```

Pass workload names to select a subset:

```sh
julia --project=docs benchmarks/overlayng_scaling.jl vertices serrated
```

`OVERLAY_SCALING_SECONDS` controls the time spent on each benchmark cell and
defaults to 0.5 seconds. Output is tab-separated so it can later feed the
benchmark plotting/reporting infrastructure directly. The harness writes no
files and is not a CI performance gate. LibGEOS's `bytes` and `allocations`
only describe its Julia result wrapper; GEOS's native heap allocations are not
visible to Chairmarks and should not be compared with OverlayNG's values.
=#

import GeoInterface as GI
import GeometryOps as GO
import LibGEOS as LG
using Chairmarks
using Statistics

const ALG = GO.OverlayNG()
const SECONDS = parse(Float64, get(ENV, "OVERLAY_SCALING_SECONDS", "0.5"))
const OPERATIONS = (:intersection, :union, :difference, :symdifference)
const WORKLOADS = isempty(ARGS) ? ("vertices", "crossings", "serrated") : Tuple(ARGS)

all(workload -> workload in ("vertices", "crossings", "serrated"), WORKLOADS) ||
    error("workloads must be selected from: vertices, crossings, serrated")

function radial_ring(n, radius; center=(0.0, 0.0), amplitude=0.0,
                     frequency=1, phase=0.0)
    cx, cy = center
    ring = Vector{Tuple{Float64,Float64}}(undef, n + 1)
    for i in 0:(n - 1)
        angle = 2pi * i / n
        r = radius * (1 + amplitude * cos(frequency * angle + phase))
        ring[i + 1] = (cx + r * cos(angle), cy + r * sin(angle))
    end
    ring[end] = ring[1]
    return ring
end

function serrated_ring(teeth; inner_radius=10.0, outer_radius=20.0, rotation=0.0)
    n = 2 * teeth
    ring = Vector{Tuple{Float64,Float64}}(undef, n + 1)
    for i in 0:(n - 1)
        angle = 2pi * i / n + rotation
        radius = iseven(i) ? outer_radius : inner_radius
        ring[i + 1] = (radius * cos(angle), radius * sin(angle))
    end
    ring[end] = ring[1]
    return ring
end

polygon(ring) = GI.Polygon([ring])

function vertex_fixtures()
    return [
        (; name="n_$n", parameter=n,
           a=polygon(radial_ring(n, 10.0)),
           b=polygon(radial_ring(n, 10.0; center=(5.0, 0.0), phase=0.017)))
        for n in (16, 64, 256, 1024, 4096, 8192)
    ]
end

function crossing_fixtures()
    n = 4096
    return [
        (; name="k_$(2 * frequency)", parameter=2 * frequency,
           a=polygon(radial_ring(n, 10.0;
                                 amplitude=0.22, frequency, phase=0.31)),
           b=polygon(radial_ring(n, 10.0)))
        for frequency in (1, 4, 16, 64, 256, 512)
    ]
end

function serrated_fixtures()
    return [
        (; name="teeth_$teeth", parameter=teeth,
           a=polygon(serrated_ring(teeth)),
           b=polygon(serrated_ring(teeth; rotation=pi / teeth)))
        for teeth in (8, 32, 128, 512, 2048)
    ]
end

function native_pair(a, b)
    return (; go_a=GO.tuples(a), go_b=GO.tuples(b),
            lg_a=GI.convert(LG, a), lg_b=GI.convert(LG, b))
end

function boundary_crossings(lg_a, lg_b)
    result = LG.intersection(LG.boundary(lg_a), LG.boundary(lg_b))
    return GI.npoint(result)
end

overlayng(::Val{:intersection}, a, b) = GO.intersection(ALG, a, b)
overlayng(::Val{:union}, a, b) = GO.union(ALG, a, b)
overlayng(::Val{:difference}, a, b) = GO.difference(ALG, a, b)
overlayng(::Val{:symdifference}, a, b) = GO.symdifference(ALG, a, b)

libgeos(::Val{:intersection}, a, b) = LG.intersection(a, b)
libgeos(::Val{:union}, a, b) = LG.union(a, b)
libgeos(::Val{:difference}, a, b) = LG.difference(a, b)
libgeos(::Val{:symdifference}, a, b) = LG.symmetricDifference(a, b)

function validate(op, inputs)
    go_result = GI.convert(LG, overlayng(op, inputs.go_a, inputs.go_b))
    lg_result = libgeos(op, inputs.lg_a, inputs.lg_b)
    LG.equals(go_result, lg_result) || error("OverlayNG and LibGEOS results differ")
    return nothing
end

function measure(op, inputs)
    go_a, go_b = inputs.go_a, inputs.go_b
    lg_a, lg_b = inputs.lg_a, inputs.lg_b
    ng_trial = @be overlayng($op, $go_a, $go_b) seconds=SECONDS
    lg_trial = @be libgeos($op, $lg_a, $lg_b) seconds=SECONDS
    return median(ng_trial), median(lg_trial)
end

function print_result(workload, fixture, segments, crossings, operation, engine, sample)
    println(join((workload, fixture.name, fixture.parameter, segments, crossings,
                  operation, engine, sample.time, sample.bytes, sample.allocs), '\t'))
end

println(join(("workload", "fixture", "parameter", "segments_per_polygon",
              "boundary_crossings", "operation", "engine", "median_seconds",
              "bytes", "allocations"), '\t'))

for workload in WORKLOADS
    fixtures = workload == "vertices"  ? vertex_fixtures() :
               workload == "crossings" ? crossing_fixtures() :
                                          serrated_fixtures()
    for fixture in fixtures
        inputs = native_pair(fixture.a, fixture.b)
        segments = GI.npoint(fixture.a) - 1
        crossings = boundary_crossings(inputs.lg_a, inputs.lg_b)
        for operation in OPERATIONS
            op = Val(operation)
            validate(op, inputs)
            ng_sample, lg_sample = measure(op, inputs)
            print_result(workload, fixture, segments, crossings,
                         operation, "OverlayNG", ng_sample)
            print_result(workload, fixture, segments, crossings,
                         operation, "LibGEOS", lg_sample)
        end
    end
end
