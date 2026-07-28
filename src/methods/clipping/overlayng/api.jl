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
`GO.f(alg, a, b)` idiom (cf. `GO.intersects(GO.RelateNG(), a, b)`), and they
deliberately do *not* accept the `target`/`fix_multipoly`/`T` arguments of the
Foster–Hormann entry points: OverlayNG has one result for a given pair of
inputs and an operation, and it is always emitted at `Float64`.

JTS's precision-model machinery (`SnappingNoder`, `SnapRoundingNoder`,
`OverlayNGRobust`, `PrecisionReducer`) is deliberately *not* ported. It exists
to patch an inexact noding substrate; this one has no such defect, so the retry
ladder would be dead weight.
=#

"""
    OverlayNG(; manifold = Planar(), exact = True())
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
struct OverlayNG{M <: Manifold, E} <: GeometryOpsCore.Algorithm{M}
    manifold::M
    exact::E
end

function OverlayNG(; manifold::Manifold = Planar(), exact = True())
    manifold isa Union{Planar, Spherical} || throw(ArgumentError(
        "OverlayNG supports the `Planar()` and `Spherical()` manifolds; got " *
        "$(typeof(manifold))"))
    return OverlayNG(manifold, exact)
end
OverlayNG(m::Manifold; kw...) = OverlayNG(; manifold = m, kw...)

GeometryOpsCore.manifold(alg::OverlayNG) = alg.manifold
GeometryOpsCore.rebuild(alg::OverlayNG, m::Manifold) = OverlayNG(m; exact = alg.exact)

# ## The four operations
#
# Each binds a user-facing name to one `_OverlayOpCode` and hands off to the
# driver. `intersection`, `union` and `difference` gain an `OverlayNG` method
# beside their existing Foster–Hormann ones; `symdifference` is new.

intersection(alg::OverlayNG, geom_a, geom_b) =
    _overlay_ng(alg.manifold, OVERLAY_INTERSECTION, geom_a, geom_b; exact = alg.exact)

union(alg::OverlayNG, geom_a, geom_b) =
    _overlay_ng(alg.manifold, OVERLAY_UNION, geom_a, geom_b; exact = alg.exact)

difference(alg::OverlayNG, geom_a, geom_b) =
    _overlay_ng(alg.manifold, OVERLAY_DIFFERENCE, geom_a, geom_b; exact = alg.exact)

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
symdifference(alg::OverlayNG, geom_a, geom_b) =
    _overlay_ng(alg.manifold, OVERLAY_SYMDIFFERENCE, geom_a, geom_b; exact = alg.exact)

symdifference(m::Manifold, geom_a, geom_b) = symdifference(OverlayNG(m), geom_a, geom_b)
symdifference(geom_a, geom_b) = symdifference(OverlayNG(Planar()), geom_a, geom_b)
