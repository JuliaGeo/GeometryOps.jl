# # Buffer
export buffer, ChenMcMains

#=
## What is a buffer?

The buffer of a geometry at distance `d` is the set of points within `d` of it.
For `d > 0` that grows the geometry; for `d < 0` it erodes an areal one, and
empties anything of lower dimension. The result is always polygonal, whatever
the input dimension was — a buffered point is a disc, a buffered line is a
sausage.

```@example buffer
import GeometryOps as GO
import GeoInterface as GI
using CairoMakie

line = GI.LineString([(0.0, 0.0), (3.0, 0.5), (4.0, 3.0), (1.5, 2.0)])
fig, ax, _ = poly(GO.buffer(line, 0.6); color = (:dodgerblue, 0.4),
                  strokecolor = :dodgerblue, strokewidth = 1.5,
                  axis = (; aspect = DataAspect(), title = "buffer(line, 0.6)"))
lines!(ax, [(GI.x(p), GI.y(p)) for p in GI.getpoint(line)]; color = :black, linewidth = 2)
fig
```

Eroding a polygon is the same call with a negative distance, and the result may
split into several parts or vanish entirely:

```@example buffer
poly_with_hole = GI.Polygon([
    GI.LinearRing([(0.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0), (0.0, 0.0)]),
    GI.LinearRing([(2.5, 0.8), (3.5, 0.8), (3.5, 3.2), (2.5, 3.2), (2.5, 0.8)]),
])
fig, ax, _ = poly(poly_with_hole; color = (:grey, 0.3),
                  axis = (; aspect = DataAspect(), title = "erosion splits the polygon"))
poly!(ax, GO.buffer(poly_with_hole, -0.5); color = (:orangered, 0.5))
fig
```

## Engines

| Call | Engine |
|:--|:--|
| `buffer(geom, d)` | [`ChenMcMains`](@ref), native Julia, planar (the default) |
| `buffer(ChenMcMains(; quadsegs = 16), geom, d)` | the same, with a finer fillet |
| `buffer(GEOS(; joinStyle = :mitre), geom, d)` | GEOS, via LibGEOS |

The native engine is round-cap, round-join only. `endCapStyle`, `joinStyle` and
`mitreLimit` are accepted at their round defaults and otherwise *rejected*, with
a message naming `GEOS()` — see [`ChenMcMains`](@ref) for why that is an error
rather than a silent forward.

## Implementation

Two stages, one file each:

1. `buffer_offset_curve.jl` builds the **raw offset curves** — one closed,
   freely self-intersecting ring per input ring, line and point, transcribed
   from JTS `OffsetSegmentGenerator` so the fillet vertices match GEOS's bit for
   bit.
2. `clipping/overlayng/winding_overlay.jl` runs the **winding-number overlay**
   over all of those rings at once and keeps the faces with winding number
   `>= 1`.

That second stage is what makes the first one simple: the loops a raw offset
curve traces at inside turns cancel themselves in the winding number, so none of
JTS's trimming, closing-segment, input-simplification or curve-inversion
heuristics are needed. The two files' headers carry the details and the places
where this deliberately disagrees with GEOS.
=#

"""
    ChenMcMains(; manifold = Planar(), quadsegs = 8, prune_eroded_rings = true,
                  single_sided = false, exact = True())
    ChenMcMains(manifold::Manifold; quadsegs = 8, prune_eroded_rings = true,
                  single_sided = false, exact = True())

The native planar [`buffer`](@ref) algorithm, and the default engine for
`buffer(geom, distance)`.

Named for Chen & McMains, *Polygon offsetting by computing winding numbers*
(ASME IDETC/CIE 2005), whose selection rule it implements: build the raw offset
curve of the input and keep the region its linework winds around at least once.
The offset curve itself is a transcription of the JTS `OffsetSegmentGenerator`
primitives, and the winding overlay runs on the exact
[`OverlayNG`](@ref) arrangement, so the arrangement's topology is exact and only
the output coordinates are rounded.

## Keywords

- `manifold`: `Planar()` only today. `Spherical()` throws; the design for it is
  recorded in `winding_overlay.jl`.
- `quadsegs = 8`: segments per quadrant of a fillet, i.e. per 90° of turn. The
  count for one corner is `round(turn / (90° / quadsegs))`, matching JTS and
  GEOS 3.14.1 exactly (both round rather than ceil, so turns below half a
  quantum get no intermediate vertex at all).
- `prune_eroded_rings = true`: drop a ring whose own erosion is already empty —
  a shell at `distance < 0`, a hole at `distance > 0` — instead of offsetting it
  (JTS `BufferCurveSetBuilder.isRingFullyEroded`). This is a speed guard, not a
  correctness one: the winding identity returns the same answer either way, and
  did so on every one of the 2275 measured collapse cases, but reaching that
  answer through the engine costs Θ(n²) self-crossings. A 1000-gon eroded past
  its apothem is 0.04 ms with the check and 69 s without. Turn it off to read
  the engine's own answer at collapse — auditing the identity, or ruling the
  check out when a result at the collapse threshold looks wrong.
- `exact = True()`: use exact predicates in the arrangement. `False()` is
  available for experiments and is not recommended.
- `endCapStyle`, `joinStyle`, `mitreLimit`: accepted **only** at their round
  defaults (`:round`, `:round`, and any `mitreLimit`, which round joins ignore).
  Any other value throws an `ArgumentError` naming the value and `GEOS()`.
- `single_sided = false`: the native buffer is two-sided. `true` throws, and
  there is no route to a single-sided buffer in GeometryOps today — the
  `GEOS()` path does not expose one either.

## Why the style keywords error

Before the native engine existed, `buffer(geom, d; kwargs...)` forwarded every
keyword to `GEOS`, so `endCapStyle = :flat` worked. The three ways to keep such
a call alive are all worse than an error:

- *silently ignore it* — a flat-capped buffer and a round-capped one are
  materially different geometries, and a caller who asked for flat caps and got
  round ones has no way to tell;
- *warn and forward to GEOS* — makes the engine, the return type, and the
  dependency on LibGEOS all depend on which keywords were passed, and fails with
  a confusing `MethodError` when LibGEOS is not loaded;
- *accept and error later* — hides the problem from `ChenMcMains`'s own
  constructor, which is where it is cheapest to report.

So this errors at construction, and the message names `GEOS()` as today's route
to caps and joins. The fields exist, and the generator has the hooks
(`_offset_end_cap!`, `_offset_corner_fillet!`), so implementing them is additive.
"""
struct ChenMcMains{M <: Manifold, E} <: GeometryOpsCore.Algorithm{M}
    manifold::M
    quadsegs::Int
    prune_eroded_rings::Bool
    endCapStyle::Symbol
    joinStyle::Symbol
    mitreLimit::Float64
    exact::E
end

# The manifold is POSITIONAL, as in `OverlayNG`: the keyword form forwards here,
# so a manifold passed positionally stays inferrable as a type parameter.
function ChenMcMains(m::Manifold; quadsegs::Integer = 8, prune_eroded_rings::Bool = true,
        endCapStyle = :round, joinStyle = :round, mitreLimit::Real = 5.0,
        single_sided::Bool = false, exact = True())
    m isa Planar || throw(ArgumentError(
        "ChenMcMains: the native buffer is implemented on the `Planar()` manifold " *
        "only; got $(typeof(m)). There is no spherical buffer in GeometryOps yet."))
    quadsegs >= 1 || throw(ArgumentError(
        "ChenMcMains: `quadsegs` must be at least 1; got $quadsegs"))
    Symbol(endCapStyle) === :round || throw(ArgumentError(
        _buffer_style_message(:endCapStyle, endCapStyle, (:flat, :square))))
    Symbol(joinStyle) === :round || throw(ArgumentError(
        _buffer_style_message(:joinStyle, joinStyle, (:mitre, :bevel))))
    single_sided === false || throw(ArgumentError(
        "ChenMcMains: `single_sided = $(repr(single_sided))` is not supported — the " *
        "native buffer is two-sided. `GEOS()` is the escape hatch for what the native " *
        "engine does not implement:\n    buffer(GEOS(; joinStyle = :mitre), geom, distance)\n" *
        "(requires `using LibGEOS`.) Single-sided buffers are not wired through that " *
        "path either, so there is no route to one in GeometryOps today."))
    return ChenMcMains(m, Int(quadsegs), prune_eroded_rings, :round, :round,
                       Float64(mitreLimit), exact)
end
ChenMcMains(; manifold::Manifold = Planar(), kw...) = ChenMcMains(manifold; kw...)

_buffer_style_message(key::Symbol, got, others) =
    "ChenMcMains: `$key = $(repr(got))` is not supported — the native buffer is " *
    "round-cap, round-join only. Use the GEOS engine for $(join(map(repr, others), " or ")) " *
    "styles:\n    buffer(GEOS(; $key = $(repr(got))), geom, distance)\n" *
    "(requires `using LibGEOS`.)"

GeometryOpsCore.manifold(alg::ChenMcMains) = alg.manifold
GeometryOpsCore.rebuild(alg::ChenMcMains, m::Manifold) =
    ChenMcMains(m; alg.quadsegs, alg.prune_eroded_rings, alg.endCapStyle, alg.joinStyle,
                alg.mitreLimit, alg.exact)

# Buffer changes the trait of what it touches — a point becomes a polygon — so it
# must match at the top-level geometry and handle Multi*/GeometryCollection
# itself. A deeper target would have `apply` rebuild e.g. a MultiPoint around
# polygons; `apply`'s own docstring warns about exactly this.
const _BUFFER_TARGETS = TraitTarget{GI.AbstractGeometryTrait}()

"""
    buffer([alg::Algorithm], geom, distance; kwargs...)
    buffer(manifold::Manifold, geom, distance; kwargs...)

The region within `distance` of `geom`: a polygonal geometry, whatever the
dimension of the input.

`distance > 0` grows `geom`; `distance < 0` erodes an areal `geom` and empties
anything of lower dimension. Distances are in the units of the input
coordinates.

A non-finite `distance` is an `ArgumentError`. A `distance` whose magnitude is
below half an ulp of the geometry's own largest coordinate cannot displace any
vertex, so it is treated as zero: an areal geometry comes back unchanged and
anything of lower dimension comes back empty, as at `distance == 0` and as in
GEOS.

## Returns

The most specific polygonal geometry that fits the result, mirroring
[`OverlayNG`](@ref):

| Result | Returned |
|:--|:--|
| one part | `GI.Polygon` |
| several parts | `GI.MultiPolygon` |
| nothing (total erosion, or a negative distance on a line or point) | an empty `GI.MultiPolygon` |

Coordinates are always `Tuple{Float64, Float64}` — the arrangement's topology is
exact, but it emits `Float64`, exactly as the overlay engine does. `GI.crs(geom)`
is carried onto the result; an `Extent` is attached only when
`calc_extent = true`.

Test an empty result with `GI.ngeom(result) == 0`, not `GI.isempty`: GeoInterface
reports `isempty(::MultiPolygon) == false` for a multi-geometry with no parts,
and every OverlayNG empty result has the same shape.

Given a vector of geometries, a Tables.jl table, or a `FeatureCollection`, the
same container comes back with each geometry buffered (this is `apply` at
`TraitTarget{GI.AbstractGeometryTrait}()`). Note that `calc_extent` then applies
to the buffered geometries, not to the rebuilt container.

## Keywords

Algorithm keywords (`quadsegs`, `prune_eroded_rings`, `exact`, and the rejected
`endCapStyle` / `joinStyle` / `mitreLimit` / `single_sided` — see [`ChenMcMains`](@ref)) may be passed directly to
`buffer(geom, distance; ...)`, which forwards them to the default algorithm's
constructor. `threaded` and any other keyword goes to `apply`; `calc_extent`
defaults to `false`.

## Examples

```jldoctest buffer
import GeoInterface as GI, GeometryOps as GO

GO.area(GO.buffer(GI.Point(0.0, 0.0), 1.0)) ≈ 32 * sin(π / 32) * cos(π / 32)

# output
true
```

A negative distance on a line is empty rather than an error:

```jldoctest buffer
import GeoInterface as GI, GeometryOps as GO

GI.ngeom(GO.buffer(GI.LineString([(0.0, 0.0), (2.0, 0.0)]), -1.0))

# output
0
```

## Engines

`buffer(geom, distance)` is the native planar [`ChenMcMains`](@ref) engine.
`buffer(GEOS(; kwargs...), geom, distance)` calls GEOS through LibGEOS and is
the route to flat or square caps and mitre or bevel joins.
"""
function buffer(geom, distance; quadsegs::Integer = 8, prune_eroded_rings::Bool = true,
        endCapStyle = :round, joinStyle = :round, mitreLimit::Real = 5.0,
        single_sided::Bool = false, exact = True(), kwargs...)
    alg = ChenMcMains(Planar(); quadsegs, prune_eroded_rings, endCapStyle, joinStyle,
                      mitreLimit, single_sided, exact)
    return buffer(alg, geom, distance; kwargs...)
end

function buffer(m::Manifold, geom, distance; quadsegs::Integer = 8,
        prune_eroded_rings::Bool = true, endCapStyle = :round, joinStyle = :round,
        mitreLimit::Real = 5.0, single_sided::Bool = false, exact = True(), kwargs...)
    alg = ChenMcMains(m; quadsegs, prune_eroded_rings, endCapStyle, joinStyle,
                      mitreLimit, single_sided, exact)
    return buffer(alg, geom, distance; kwargs...)
end

function buffer(alg::ChenMcMains, geom, distance; calc_extent = false, kwargs...)
    d = Float64(distance)
    isfinite(d) || throw(ArgumentError(
        "buffer: `distance` must be finite; got $distance. An infinite buffer has no " *
        "polygonal boundary, and a NaN one has no meaning."))
    ce = booltype(calc_extent)
    return apply(_BUFFER_TARGETS, geom; kwargs...) do g
        _buffer(alg, g, d, ce)
    end
end

function _buffer(alg::ChenMcMains{<:Planar}, geom, d::Float64, calc_extent::BoolsAsTypes)
    T = _kernel_point_type(alg.manifold)
    rings = _raw_offset_curves(T, geom, d, alg.quadsegs; alg.prune_eroded_rings)
    polys = isempty(rings) ? _result_poly_type(T)[] :
        _winding_overlay(alg.manifold, T, rings;
                         keep = _winding_at_least_one, exact = alg.exact)
    return _buffer_result(T, polys, GI.crs(geom), calc_extent)
end

#=
The most specific geometry over the result polygons, with the input's CRS and an
optional extent. `_build_polygons` emits CRS-less, extent-less wrappers (the
overlay engine attaches neither), so a single-polygon result is re-wrapped around
the same ring vector rather than copied.
=#
function _buffer_result(::Type{T}, polys, crs, calc_extent::BoolsAsTypes) where {T}
    if isempty(polys)
        #-- an empty result is a polygonal geometry with no parts, NOT an error and
        #-- not a zero-ring `Polygon`: `rebuild` reads Z/M off `first(child_geoms)`
        return GI.MultiPolygon{false, false, Vector{_result_poly_type(T)}, typeof(crs),
                               Nothing}(_result_poly_type(T)[], crs, nothing)
    elseif length(polys) == 1
        p = polys[1]
        return GI.Polygon(GI.getgeom(p); crs, extent = _buffer_extent(p, calc_extent))
    end
    return GI.MultiPolygon(polys; crs,
        extent = _buffer_extent(GI.MultiPolygon(polys), calc_extent))
end

_buffer_extent(geom, ::True) = GI.extent(geom)
_buffer_extent(geom, ::False) = nothing

# Add an error hint for `buffer` if LibGEOS is not loaded!
function _buffer_error_hinter(io, exc, argtypes, kwargs)
    if isnothing(Base.get_extension(GeometryOps, :GeometryOpsLibGEOSExt)) && exc.f == buffer && first(argtypes) == GEOS
        print(io, "\n\nThe `buffer` method requires the LibGEOS.jl package to be explicitly loaded.\n")
        print(io, "You can do this by simply typing ")
        printstyled(io, "using LibGEOS"; color = :cyan, bold = true)
        println(io, " in your REPL, \nor otherwise loading LibGEOS.jl via using or import.")
    end
end
