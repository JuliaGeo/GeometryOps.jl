# # Controlled OverlayNG scaling benchmarks
#
#=
This harness separates three sources of geometry-operation cost:

- `vertices`: two convex polygons with two boundary crossings while their
  vertex counts grow;
- `crossings`: two 4,096-segment polygons while their boundary-crossing count
  grows from 2 to 1,024;
- `serrated`: a radial saw polygon against a half-tooth-rotated copy, so both
  input size and the number of local tooth intersections grow together.

For every fixture it measures:

- the four constructive overlay operations with OverlayNG and LibGEOS;
- `intersects` with GeometryOps RelateNG and LibGEOS, prepared and unprepared.

Prepared geometries cannot be operands of constructive overlay in either
library. GeometryOps preparation is currently implemented only by `RelateNG`,
and GEOS prepared geometries expose predicates and relate operations, not
intersection/union/difference/symmetric difference. The prepared series is
therefore a predicate comparison rather than a misleading "prepared overlay".

Geometry construction, conversion, validation, preparation, and boundary-
crossing counts happen outside the timed expressions. GeometryOps receives
`GO.tuples` geometries and LibGEOS receives `LibGEOS.Polygon`s.

Run every workload with:

```sh
julia --project=docs benchmarks/overlayng_scaling.jl
```

Pass workload names to select a subset:

```sh
julia --project=docs benchmarks/overlayng_scaling.jl vertices serrated
```

From a Julia session, use `include("benchmarks/overlayng_scaling.jl")` to keep
the resulting `results::DataFrame` and `figures` named tuple available for
further analysis.

`OVERLAY_SCALING_SECONDS` controls the time spent on each benchmark cell and
defaults to 0.5 seconds. Results remain in the `results` DataFrame and the
AlgebraOfGraphics output remains in `figures`; the harness writes no files and
is not a CI performance gate. LibGEOS's `bytes` and `allocations` only describe
its Julia result wrapper: GEOS native heap allocations are invisible to
Chairmarks and should not be compared with OverlayNG's values.
=#

import GeoInterface as GI
import GeometryOps as GO
import LibGEOS as LG
using AlgebraOfGraphics
using CairoMakie
using Chairmarks
using DataFrames
using Statistics

const OVERLAY_ALG = GO.OverlayNG()
const RELATE_ALG = GO.RelateNG()
const SECONDS = parse(Float64, get(ENV, "OVERLAY_SCALING_SECONDS", "0.5"))
const OVERLAY_OPERATIONS = (:intersection, :union, :difference, :symdifference)
const WORKLOADS = isempty(ARGS) ? ("vertices", "crossings", "serrated") : Tuple(ARGS)

all(workload -> workload in ("vertices", "crossings", "serrated"), WORKLOADS) ||
    error("workloads must be selected from: vertices, crossings, serrated")

const ResultRow = NamedTuple{
    (:family, :workload, :fixture, :parameter, :segments_per_polygon,
     :boundary_crossings, :operation, :engine, :preparation, :median_seconds,
     :bytes, :allocations),
    Tuple{String,String,String,Int,Int,Int,String,String,String,Float64,Float64,Float64},
}

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

overlayng(::Val{:intersection}, a, b) = GO.intersection(OVERLAY_ALG, a, b)
overlayng(::Val{:union}, a, b) = GO.union(OVERLAY_ALG, a, b)
overlayng(::Val{:difference}, a, b) = GO.difference(OVERLAY_ALG, a, b)
overlayng(::Val{:symdifference}, a, b) = GO.symdifference(OVERLAY_ALG, a, b)

libgeos(::Val{:intersection}, a, b) = LG.intersection(a, b)
libgeos(::Val{:union}, a, b) = LG.union(a, b)
libgeos(::Val{:difference}, a, b) = LG.difference(a, b)
libgeos(::Val{:symdifference}, a, b) = LG.symmetricDifference(a, b)

function validate_overlay(op, inputs)
    go_result = GI.convert(LG, overlayng(op, inputs.go_a, inputs.go_b))
    lg_result = libgeos(op, inputs.lg_a, inputs.lg_b)
    LG.equals(go_result, lg_result) || error("OverlayNG and LibGEOS results differ")
    return nothing
end

function result_row(family, workload, fixture, segments, crossings,
                    operation, engine, preparation, sample)::ResultRow
    return (; family=String(family), workload=String(workload),
            fixture=String(fixture.name), parameter=Int(fixture.parameter),
            segments_per_polygon=Int(segments), boundary_crossings=Int(crossings),
            operation=String(operation), engine=String(engine),
            preparation=String(preparation), median_seconds=Float64(sample.time),
            bytes=Float64(sample.bytes), allocations=Float64(sample.allocs))
end

function measure_overlay!(rows, workload, fixture, inputs, segments, crossings)
    go_a, go_b = inputs.go_a, inputs.go_b
    lg_a, lg_b = inputs.lg_a, inputs.lg_b
    for operation in OVERLAY_OPERATIONS
        op = Val(operation)
        validate_overlay(op, inputs)
        ng_sample = median(@be overlayng($op, $go_a, $go_b) seconds=SECONDS)
        lg_sample = median(@be libgeos($op, $lg_a, $lg_b) seconds=SECONDS)
        push!(rows, result_row("overlay", workload, fixture, segments, crossings,
                               operation, "OverlayNG", "not applicable", ng_sample))
        push!(rows, result_row("overlay", workload, fixture, segments, crossings,
                               operation, "LibGEOS", "not applicable", lg_sample))
    end
    return nothing
end

function measure_predicates!(rows, workload, fixture, inputs, segments, crossings)
    go_a, go_b = inputs.go_a, inputs.go_b
    lg_a, lg_b = inputs.lg_a, inputs.lg_b
    go_prepared = GO.prepare(RELATE_ALG, go_a; validate=false)
    lg_prepared = LG.prepareGeom(lg_a)
    predicate = GO.pred_intersects()

    answers = (
        GO.intersects(RELATE_ALG, go_a, go_b),
        GO.relate_predicate(go_prepared, predicate, go_b),
        LG.intersects(lg_a, lg_b),
        LG.intersects(lg_prepared, lg_b),
    )
    all(==(first(answers)), answers) || error("intersects implementations disagree")

    samples = (
        ("GeometryOps", "unprepared",
         median(@be GO.intersects($RELATE_ALG, $go_a, $go_b) seconds=SECONDS)),
        ("GeometryOps", "prepared",
         median(@be GO.relate_predicate($go_prepared, $predicate, $go_b) seconds=SECONDS)),
        ("LibGEOS", "unprepared",
         median(@be LG.intersects($lg_a, $lg_b) seconds=SECONDS)),
        ("LibGEOS", "prepared",
         median(@be LG.intersects($lg_prepared, $lg_b) seconds=SECONDS)),
    )
    for (engine, preparation, sample) in samples
        push!(rows, result_row("predicate", workload, fixture, segments, crossings,
                               "intersects", engine, preparation, sample))
    end
    return nothing
end

function run_benchmarks(workloads=WORKLOADS)
    rows = ResultRow[]
    for workload in workloads
        fixtures = workload == "vertices"  ? vertex_fixtures() :
                   workload == "crossings" ? crossing_fixtures() :
                                              serrated_fixtures()
        for fixture in fixtures
            inputs = native_pair(fixture.a, fixture.b)
            segments = GI.npoint(fixture.a) - 1
            crossings = boundary_crossings(inputs.lg_a, inputs.lg_b)
            measure_overlay!(rows, workload, fixture, inputs, segments, crossings)
            measure_predicates!(rows, workload, fixture, inputs, segments, crossings)
        end
    end
    return DataFrame(rows)
end

function make_figures(results)
    overlay_results = subset(results, :family => ByRow(==("overlay")))
    predicate_results = subset(results, :family => ByRow(==("predicate")))

    overlay_plot = data(overlay_results) *
        mapping(:parameter => "Workload parameter",
                :median_seconds => "Median time (s)";
                col=:workload => "Workload",
                row=:operation => "Operation",
                color=:engine => "Engine") *
        visual(ScatterLines)

    predicate_plot = data(predicate_results) *
        mapping(:parameter => "Workload parameter",
                :median_seconds => "Median time (s)";
                col=:workload => "Workload",
                row=:preparation => "Preparation",
                color=:engine => "Engine") *
        visual(ScatterLines)

    axis = (; xscale=log10, yscale=log10)
    return (;
        overlay=draw(overlay_plot; axis,
                     figure=(; size=(1400, 1400), title="Overlay scaling")),
        predicates=draw(predicate_plot; axis,
                        figure=(; size=(1400, 700), title="Prepared predicate scaling")),
    )
end

results = run_benchmarks()
figures = make_figures(results)
println("Collected $(nrow(results)) benchmark rows in `results`.")
display(figures.overlay)
display(figures.predicates)
