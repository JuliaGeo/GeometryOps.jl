# # OverlayNG — the exact-arrangement overlay algorithm

export OverlayNG, symdifference

#=
## What is OverlayNG?

`OverlayNG` is an `Algorithm` for the four set-theoretic overlay
operations on geometries — intersection, union, difference and symmetric
difference — computed as an **exact arrangement**: the two inputs are noded
against each other with exact predicates, the crossings are carried as symbolic
nodes, the arrangement is labelled and assembled as a half-edge graph, and
Float64 coordinates are produced *once*, at emission. No constructed coordinate
ever feeds a decision, so there is no snapping, no tolerance, and no precision
model anywhere in the pipeline.

It is **opt-in**: pass it as the first argument.

```@example overlayng
import GeometryOps as GO
import GeoInterface as GI

a = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
b = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])

GO.intersection(GO.OverlayNG(), a, b)
```

Unlike the default Foster–Hormann engine, an overlay returns a *single*
GeoInterface geometry — the most specific one that fits the result — rather than
a list, and it can return components of several dimensions at once:

```@example overlayng
GO.symdifference(GO.OverlayNG(), a, b)
```

```@example overlayng
#-- two boxes sharing only an edge: the intersection is a line, not an area
c = GI.Polygon([[(2.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 0.0)]])
GO.intersection(GO.OverlayNG(), a, c)
```

If you only ever want one dimension back, say so with `target` and the result
shape stops depending on the data — this pair shares only a line, so the areal
answer is the empty vector rather than a `LineString` you did not ask for:

```@example overlayng
GO.intersection(GO.OverlayNG(), a, c; target = GI.PolygonTrait())
```

Point inputs work too, which makes intersection an efficient "which of these
points are inside" query:

```@example overlayng
GO.intersection(GO.OverlayNG(), GI.MultiPoint([(0.5, 0.5), (5.0, 5.0)]), a)
```

Every operation is available on the sphere by giving the algorithm a manifold,
and the spherical answers are computed with spherical predicates throughout
(great-circle edges, exact point-in-area location), not by overlaying longitude
and latitude as if they were `x` and `y`:

```@example overlayng
sph = GO.OverlayNG(GO.Spherical())
u = GO.union(sph, a, b)
i = GO.intersection(sph, a, b)
#-- area is conserved to machine precision: area(A∪B) + area(A∩B) == area(A) + area(B)
(GO.area(GO.Spherical(), u) + GO.area(GO.Spherical(), i)) /
    (GO.area(GO.Spherical(), a) + GO.area(GO.Spherical(), b)) - 1
```

Spherical results come back in the chart the engine works in — unit-sphere
`xyz`, as `UnitSphericalPoint`s — so a result vertex that is an input vertex is
bit-for-bit the ingested one, with no lon/lat round trip to lose the last ULP
in. Ask for `point_type = Tuple{Float64,Float64}` to get `(lon, lat)` degrees
back instead:

```@example overlayng
GI.coordinates(GO.intersection(GO.OverlayNG(GO.Spherical(); point_type = Tuple{Float64,Float64}), a, b))
```

## Implementation

This file is the public surface of the OverlayNG port; everything it calls is
internal. The engine lives in the sibling files of this directory and is a
file-by-file port of JTS's `operation/overlayng` package: the noding substrate
under `noding/`, then `overlay_label.jl` → `overlay_graph.jl` →
`overlay_labeller.jl` → `maximal_edge_ring.jl`/`polygon_builder.jl` →
`line_builder.jl`/`intersection_point_builder.jl` →
`overlay_points.jl`/`overlay_mixed_points.jl`, tied together by the
`_overlay_ng` driver in `overlay_ng.jl` (the port of `OverlayNG.getResult`).

The four op codes are the internal `_OverlayOpCode` enum; the methods below are
the only place a user-facing name is bound to one. They follow the house
`GO.f(alg, a, b)` idiom (cf. `GO.intersects(GO.RelateNG(), a, b)`). Of the
Foster–Hormann entry points' extra arguments they take `target` — with the same
meaning and the same return shapes, see the docstring — but not
`fix_multipoly` (there is nothing to fix: the arrangement never emits
overlapping components) and not `T` (the result is always emitted at `Float64`).

JTS's precision-model machinery (`SnappingNoder`, `SnapRoundingNoder`,
`OverlayNGRobust`, `PrecisionReducer`) is deliberately *not* ported. It exists
to patch an inexact noding substrate; this one has no such defect, so the retry
ladder would be dead weight.
=#

"""
    OverlayNG(; manifold = Planar(), exact = True(), point_type = ...)
    OverlayNG(manifold::Manifold; kwargs...)

The exact-arrangement overlay algorithm, a port of the JTS OverlayNG engine by
Martin Davis, extended to the sphere.

`OverlayNG` computes [`intersection`](@ref), `union`, [`difference`](@ref) and
[`symdifference`](@ref) of two geometries of any dimension. It is **opt-in**:
the algorithm is the first argument.

```julia
GO.intersection(GO.OverlayNG(), a, b)
GO.union(GO.OverlayNG(GO.Spherical()), a, b)
```

Foster–Hormann clipping remains the default engine for `intersection`, `union`
and [`difference`](@ref) when no algorithm is given; those defaults are
unchanged by the presence of this algorithm. [`symdifference`](@ref) is the one
exception — see its docstring.

## Keyword arguments

- `manifold`: `Planar()` (default) or `Spherical()`. `Spherical(; oriented)`
  selects how a ring denotes its region, and overlay follows that choice
  throughout — nothing extra is needed here.
- `exact`: `True()` (default) to decide every uncertain filter with an exact
  predicate, `False()` to stay in Float64. Leave it at the default unless you
  are measuring the cost of exactness.
- `point_type`: the type of the coordinates in the result. Defaults to the
  manifold's own working point type — `Tuple{Float64,Float64}` on `Planar()`,
  `UnitSphericalPoint{Float64}` (3D unit-sphere `xyz`) on `Spherical()`. On the
  sphere `Tuple{Float64,Float64}` is also accepted and gives `(lon, lat)`
  degrees; see "Output coordinates" below for what that costs.

## Output coordinates

The arrangement is exact and symbolic: rounding to Float64 happens exactly once,
when a node's coordinate is realized for output. `point_type` chooses the chart
that rounding lands in, and on the sphere the choice is not neutral.

`Spherical()` works in unit-sphere `xyz` throughout — an input vertex is
converted to a `UnitSphericalPoint` at ingest and every predicate reads it there
— so emitting `UnitSphericalPoint{Float64}` is a pass-through for every result
vertex that is an input vertex: the coordinate that comes out is *bit-for-bit*
the one that went in, and `GO.area` of a clipped cell agrees with the uncut
original exactly rather than in all but the last ULP.

Emitting `(lon, lat)` instead sends that vertex back out through `atan`/`asin`,
which is a rounding of a value that had an exact image in the format it was
already in. Measured over 200 000 uniformly random directions, the round trip
displaces a point by up to 3.2 ULPs of the unit sphere below 60° latitude, 7.9
at 75°, and 126 above 89.5°; swept along single parallels the worst displacement
is 7.8 ULPs at 89°, 1 364 at 89.99° and 2 915 (6.5e-13 rad) at 89.999°. The
growth is the chart's, not the arithmetic's — a degree of longitude is
`cos φ` of an arc — and it is why a polar grid cell survives an overlay in
`xyz` and does not in lon/lat.

For a crossing node there is no exact Float64 answer in either chart: the
position is a `Rational{BigInt}` direction with no finite decimal form. What
`xyz` buys there is rounding once (normalize the direction) rather than twice
(normalize, then trigonometry).

`Tuple{Float64,Float64}` therefore exists for callers that need lon/lat
coordinates back and are willing to pay for them, not as a lossless alternative.

## Result shape

The result is a single GeoInterface geometry, the most specific one that fits:

- one component of one dimension → `Point` / `LineString` / `Polygon`;
- several components of one dimension → the corresponding `Multi` geometry;
- components of several dimensions → a `GeometryCollection`, ordered areas,
  then lines, then points;
- nothing → an empty `MultiPoint`/`MultiLineString`/`MultiPolygon`, whose
  dimension follows the OGC rule for the operation (`min` of the inputs for
  intersection, `max` for union and symmetric difference, the left input's for
  difference).

Lower-dimension components are *included*: the intersection of two polygons
that share a boundary segment and also overlap is a `GeometryCollection` of the
overlap polygon and the shared line. This is JTS's original (non-strict)
overlay semantics.

## `target`: asking for one dimension

Each operation takes a `target` keyword that narrows the result to a single
dimension, chosen up front rather than by the data. It is the `target` of the
Foster–Hormann entry points, with the same return shapes — a singular trait
gives the `Vector` of atomic components, a `Multi` trait gives the one
multi-geometry:

| `target`                    | returns                |
|:----------------------------|:-----------------------|
| `nothing` (default)         | as above — most specific over every dimension |
| `GI.PolygonTrait()`         | `Vector{<:Polygon}`    |
| `GI.MultiPolygonTrait()`    | `MultiPolygon`         |
| `GI.LineStringTrait()`      | `Vector{<:LineString}` |
| `GI.MultiLineStringTrait()` | `MultiLineString`      |
| `GI.PointTrait()`           | `Vector{<:Point}`      |
| `GI.MultiPointTrait()`      | `MultiPoint`           |

```julia
#-- always a MultiPolygon, whatever `a` and `b` turn out to share
GO.intersection(GO.OverlayNG(), a, b; target = GI.MultiPolygonTrait())
```

An empty targeted result is the empty `Vector` or empty multi-geometry of that
same concrete type, so the return type no longer depends on the inputs — which
is the point: without `target`, code that only wants areas has to handle a
`GeometryCollection` that appears only when the inputs happen to touch.

`target` also removes work. A target above the result's OGC dimension is
answered from the input dimensions alone, with no noding at all (an areal target
on a line ∩ area is empty for every possible input), and the builds the target
cannot want are skipped. The saving is one-directional, because the three builds
form a dependency chain: an areal target skips both the line and point builds, a
line target skips the point build, and a point target skips neither. Result
polygons are built even for a line or point target — whether a line lies inside
the result area is part of deciding it, and there is no cheaper equivalent test.

## Inputs

`Point`, `MultiPoint`, `LineString`, `LinearRing`, `MultiLineString`, `Polygon`
and `MultiPolygon`, in any A × B combination. `GeometryCollection` inputs raise
an `ArgumentError`.

Inputs must be **valid**: rings may not self-cross, and the components of a
multi-geometry may not overlap. The engine nodes A against B but never A against
itself, so an invalid input yields an undefined result rather than an error. Fix
inputs first (e.g. with `CrossingEdgeSplit` / `AntipodalEdgeSplit` on the
sphere) if you are not sure.

## Robustness

All topological decisions are made by exact predicates on the *input*
coordinates and on symbolic crossing keys. Float64 values are used only as
filters with certified error bounds, and only appear in the output. There is
therefore no snapping, no tolerance and no precision model to configure, and no
robustness failure mode to retry around.

## Limitation: the full sphere

On `Spherical()`, an overlay whose result covers the entire sphere and has no
boundary at all cannot be returned, and raises an `ArgumentError`. A polygon
denotes the region bounded by its rings, and a polygon with no rings already
means the empty geometry, so GeometryOps has no spelling for the full sphere.
Reformulate such an operation as a `difference` from the covering region.
"""
struct OverlayNG{M <: Manifold, E, P} <: GeometryOpsCore.Algorithm{M}
    manifold::M
    exact::E
    point_type::Type{P}
end

# The manifold is POSITIONAL here and the keyword form forwards to it, not the
# other way round: `point_type`'s default is `_kernel_point_type(manifold)`, and
# a keyword argument is not specialized on, so computing it in the keyword-only
# method leaves the algorithm's third parameter — and with it the return type of
# every targeted overlay — inferred as `Type` rather than as a constant.
function OverlayNG(m::Manifold; exact = True(), point_type = nothing)
    #-- `nothing` rather than `_kernel_point_type(m)` as the default: a keyword
    #-- default is evaluated before the body, and an unsupported manifold has no
    #-- kernel point type — spelling the default here would replace the message
    #-- below with a `MethodError` from RelateNG's ingest.
    m isa Union{Planar, Spherical} || throw(ArgumentError(
        "OverlayNG supports the `Planar()` and `Spherical()` manifolds; got " *
        "$(typeof(m))"))
    pt = point_type === nothing ? _kernel_point_type(m) : point_type
    _overlay_supports_point_type(m, pt) ||
        throw(ArgumentError(_bad_point_type_msg(m, pt)))
    return OverlayNG(m, exact, pt)
end
OverlayNG(; manifold::Manifold = Planar(), kw...) = OverlayNG(manifold; kw...)

#=
The output point types each manifold can emit. Both lists are closed, and
deliberately so: the emitter has one method per (kernel point, output point)
pair, and every one of them is a rounding argument that had to be made
explicitly (noding/emit.jl). A manifold/point-type pair that is not listed here
has no such argument behind it, so it is rejected rather than approximated.
=#
_overlay_supports_point_type(::Planar, ::Type{Tuple{Float64, Float64}}) = true
_overlay_supports_point_type(::Spherical, ::Type{Tuple{Float64, Float64}}) = true
_overlay_supports_point_type(::Spherical, ::Type{UnitSphericalPoint{Float64}}) = true
_overlay_supports_point_type(::Manifold, ::Any) = false

_bad_point_type_msg(m::Manifold, pt) =
    "OverlayNG: `point_type = $(pt)` is not a supported output coordinate type " *
    "on $(typeof(m).name.name)(). $(m isa Planar ?
        "`Planar()` emits `Tuple{Float64,Float64}` (x, y) only." :
        "`Spherical()` emits `UnitSphericalPoint{Float64}` (unit-sphere xyz, the " *
        "default and the chart the engine works in) or `Tuple{Float64,Float64}` " *
        "((lon, lat) degrees, one extra rounding).")"

GeometryOpsCore.manifold(alg::OverlayNG) = alg.manifold

# The output point type is manifold-derived by DEFAULT, so a rebuild onto another
# manifold re-derives it; only a choice that deviates from the old manifold's
# default is a choice the caller made, and that one is carried across whenever
# the new manifold can honour it. (The single deviation available today —
# `Spherical()` emitting lon/lat tuples — is the planar default, so it survives
# every rebuild.)
GeometryOpsCore.rebuild(alg::OverlayNG, m::Manifold) =
    OverlayNG(m; exact = alg.exact, point_type = _rebuilt_point_type(alg, m))

_rebuilt_point_type(alg::OverlayNG, m::Manifold) =
    (alg.point_type === _kernel_point_type(alg.manifold) ||
     !_overlay_supports_point_type(m, alg.point_type)) ?
        _kernel_point_type(m) : alg.point_type

# ## The four operations
#
# Each binds a user-facing name to one `_OverlayOpCode` and hands off to the
# driver. `intersection`, `union` and `difference` gain an `OverlayNG` method
# beside their existing Foster–Hormann ones; `symdifference` is new.

intersection(alg::OverlayNG, geom_a, geom_b; target = nothing) =
    _overlay_ng(alg.manifold, OVERLAY_INTERSECTION, geom_a, geom_b;
                exact = alg.exact, point_type = alg.point_type, target)

union(alg::OverlayNG, geom_a, geom_b; target = nothing) =
    _overlay_ng(alg.manifold, OVERLAY_UNION, geom_a, geom_b;
                exact = alg.exact, point_type = alg.point_type, target)

difference(alg::OverlayNG, geom_a, geom_b; target = nothing) =
    _overlay_ng(alg.manifold, OVERLAY_DIFFERENCE, geom_a, geom_b;
                exact = alg.exact, point_type = alg.point_type, target)

"""
    symdifference([alg::OverlayNG], geom_a, geom_b)
    symdifference(manifold::Manifold, geom_a, geom_b)

The symmetric difference of two geometries: everything that lies in exactly one
of them, i.e. `union(difference(a, b), difference(b, a))`.

```jldoctest
import GeoInterface as GI, GeometryOps as GO

a = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
b = GI.Polygon([[(1.0, 0.0), (3.0, 0.0), (3.0, 2.0), (1.0, 2.0), (1.0, 0.0)]])
GO.area(GO.symdifference(a, b))

# output
4.0
```

## Engine

`symdifference` is **OverlayNG-only**: there is no Foster–Hormann symmetric
difference, and none is planned — symmetric difference is the operation the
exact arrangement gets for free and the ent/exit tracer does not.

That makes it the one member of the overlay family whose algorithm-free form
does not run Foster–Hormann: `symdifference(a, b)` is
`symdifference(OverlayNG(Planar()), a, b)`, and `symdifference(m, a, b)` is
`symdifference(OverlayNG(m), a, b)`. The defaults of [`intersection`](@ref),
`union` and [`difference`](@ref) are unaffected.

See [`OverlayNG`](@ref) for the result shape, the input contract, and the
spherical full-sphere limitation.
"""
symdifference(alg::OverlayNG, geom_a, geom_b; target = nothing) =
    _overlay_ng(alg.manifold, OVERLAY_SYMDIFFERENCE, geom_a, geom_b;
                exact = alg.exact, point_type = alg.point_type, target)

symdifference(m::Manifold, geom_a, geom_b; kwargs...) =
    symdifference(OverlayNG(m), geom_a, geom_b; kwargs...)
symdifference(geom_a, geom_b; kwargs...) =
    symdifference(OverlayNG(Planar()), geom_a, geom_b; kwargs...)
