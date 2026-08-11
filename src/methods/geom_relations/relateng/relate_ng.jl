# # Relate

export relate, RelateNG, prepare

#=
## What is relate?

The `relate` function computes the full topological relationship between two
geometries — not just a single yes/no question like [`intersects`](@ref) or
[`within`](@ref), but the complete description from which *every* such
question can be answered. That description is the
**dimensionally extended nine-intersection model (DE-9IM)** matrix.

To provide an example, consider these two overlapping polygons:
```@example relateng
import GeometryOps as GO
import GeoInterface as GI
using Makie
using CairoMakie

p1 = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
p2 = GI.Polygon([[(2.0, 2.0), (5.0, 2.0), (5.0, 5.0), (2.0, 5.0), (2.0, 2.0)]])
f, a, p = poly(collect(GI.getpoint(p1)); color = (:blue, 0.5), axis = (; aspect = DataAspect()))
poly!(collect(GI.getpoint(p2)); color = (:orange, 0.5))
f
```
Their relationship is captured by a single matrix:
```@example relateng
GO.relate(p1, p2)
```

## What is the DE-9IM?

Every geometry partitions the plane into three point sets: its *interior*,
its *boundary*, and its *exterior*. The DE-9IM matrix records, for each of
the nine pairings of those sets between geometry A (rows) and geometry B
(columns), the dimension of the pair's intersection: `F` for empty, `0` for
points, `1` for lines, `2` for areas. The matrix above, `"212101212"`, reads
row-major:

|                | B interior | B boundary | B exterior |
|----------------|------------|------------|------------|
| **A interior** | `2`        | `1`        | `2`        |
| **A boundary** | `1`        | `0`        | `1`        |
| **A exterior** | `2`        | `1`        | `2`        |

The interiors meet in an area (`2`), the boundaries cross at points (`0`),
and each geometry has interior outside the other (`2` in the exterior
column/row) — exactly the *overlaps* relationship. Named predicates are
just patterns over this matrix, where `T` means "non-empty (any dimension)"
and `*` means "don't care". You can match a pattern directly, which
short-circuits as soon as the answer is known instead of computing the full
matrix:
```@example relateng
GO.relate(p1, p2, "T*T***T**") # the `overlaps` pattern for two areas
```
or evaluate a named predicate through the [`RelateNG`](@ref) algorithm:
```@example relateng
GO.overlaps(GO.RelateNG(), p1, p2)
```
When one geometry is tested against many others, [`prepare`](@ref) it once
and reuse the cached indexes:
```@example relateng
prep = GO.prepare(GO.RelateNG(), p1)
GO.relate(prep, p2)
```

## Implementation

This file is the port of JTS `RelateNG.java` — the evaluation engine that
drives all the RelateNG machinery (the kernel, point locators, topology
computer, and node analysis in the surrounding files) to compute the value
of a topological predicate between two geometries, based on the DE-9IM
model.

Method order parallels the Java file (`evaluate`, `hasRequiredEnvelopeInteraction`,
`finishValue`, `computePP`, `computeAtPoints`, `computePoints`, `computePoint`,
`computeLineEnds`, `computeLineEnd`, `computeAreaVertex` ×2, `computeAtEdges`,
`computeEdgesAll`, `computeEdgesMutual`), so this file diffs against its
Java counterpart. Idiom changes:

- The Java class holds `geomA` (for prepared mode); here the unprepared
  entry points build both `RelateGeometry`s per call, while prepared mode
  carries the A side — with its lazy caches forced and the segment
  strings/segment tree prebuilt — in a [`PreparedRelate`](@ref), threaded
  through the evaluation as the optional `prep` argument.
- The algorithm configuration (manifold, accelerator, exactness flag,
  boundary node rule) travels in the `RelateNG` algorithm struct, the
  house `Algorithm{M}` idiom (cf. `FosterHormannClipping`).
- Java's `computeEdgesAll` feeds the combined A∪B edge set through one
  `EdgeSetIntersector` (every unordered chain pair once). Here that is
  phased as A×B (`process_edge_intersections!`), then A×A and B×B
  (`process_self_intersections!`) — the same pair set, possibly visited
  in a different order, which only affects *when* a short-circuiting
  predicate exits, never the final value (dimension updates are
  monotone).
- Java's null `Envelope` is `nothing` here: empty geometries have a
  `nothing` extent, `ext_intersects`/`ext_covers` return `false` for it
  (as Java's null envelope does), and `computeAtEdges` early-returns
  before ever forwarding an extent filter that could be `nothing`.
=#

"""
    RelateNG{M <: Manifold, A <: IntersectionAccelerator, E, BR <: BoundaryNodeRule}

The next-generation DE-9IM topological-relationship algorithm, a port of the
JTS RelateNG algorithm by Martin Davis. Capabilities:

1. Efficient short-circuited evaluation of topological predicates (including
   matching custom DE-9IM matrix patterns).
2. Robust evaluation: all answers are computed from exact predicates on the
   input coordinates (no constructed intersection points), so invalid
   topology does not cause failures.
3. `GeometryCollection` inputs containing mixed types and overlapping
   polygons are supported, using *union semantics*.
4. Zero-length LineStrings are treated as being topologically identical to
   Points.
5. Support for [`BoundaryNodeRule`](@ref)s.

All coordinates are evaluated as `Float64`: input coordinates are converted
on extraction, and the exact-predicate machinery (adaptive orientation
predicates, rational-arithmetic node coincidence) assumes `Float64` inputs.
Non-`Float64` geometries are accepted but evaluated at `Float64` precision.

Keyword arguments (all optional): `manifold` (default `Planar()`),
`accelerator` (default [`AutoAccelerator`](@ref)), `exact` (default
`True()`), `boundary_rule` (default [`Mod2Boundary`](@ref), the OGC SFS
rule).

## Unprepared performance

Every unprepared evaluation rebuilds an extent-annotated view of both
inputs (the stand-in for the envelope cache JTS keeps on each Geometry),
one coordinate pass per call. Inputs that already carry extents at every
level skip that pass — stamp them once with

```julia
geom = GO.tuples(geom; calc_extent = true)
```

and repeated unprepared calls read the stored extents instead. When one
geometry is queried many times, [`prepare`](@ref) it instead: the prepared
form also caches the point locators and edge index.

See [`relate`](@ref) and [`relate_predicate`](@ref) for the entry points.
"""
struct RelateNG{M <: Manifold, A <: IntersectionAccelerator, E, BR <: BoundaryNodeRule} <: GeometryOpsCore.Algorithm{M}
    manifold::M
    accelerator::A
    exact::E
    boundary_rule::BR
end
RelateNG(; manifold::Manifold = Planar(), accelerator = AutoAccelerator(),
        exact = True(), boundary_rule = Mod2Boundary()) =
    RelateNG(manifold, accelerator, exact, boundary_rule)
RelateNG(m::Manifold; kw...) = RelateNG(; manifold = m, kw...)
GeometryOpsCore.manifold(alg::RelateNG) = alg.manifold

#==========================================================================
## Entry points (the static RelateNG.relate(...) overloads)
==========================================================================#

"""
    relate([alg::RelateNG], a, b)::DE9IM

Computes the [`DE9IM`](@ref) matrix for the topological relationship between
geometries `a` and `b`.

Port of `RelateNG.relate(Geometry, Geometry)`.
"""
function relate(alg::RelateNG, a, b)
    pred = RelateMatrixPredicate()
    relate_predicate(alg, pred, a, b)
    return result_im(pred)
end
relate(a, b) = relate(RelateNG(), a, b)

"""
    relate([alg::RelateNG], a, b, im_pattern::AbstractString)::Bool

Tests whether the topological relationship between geometries `a` and `b`
matches the DE-9IM matrix pattern `im_pattern` (9 characters over
`012TF*`).

Port of `RelateNG.relate(Geometry, Geometry, String)`.
"""
relate(alg::RelateNG, a, b, im_pattern::AbstractString) =
    relate_predicate(alg, pred_matches(im_pattern), a, b)
relate(a, b, im_pattern::AbstractString) = relate(RelateNG(), a, b, im_pattern)

# Method-ambiguity guard between `relate(alg, a, b)` and
# `relate(a, b, pattern)`: a String third argument with a RelateNG first
# argument is missing its B geometry.
relate(alg::RelateNG, a, im_pattern::AbstractString) =
    throw(ArgumentError("`relate(alg, a, b, pattern)` requires both geometries; got a pattern in place of the B geometry"))

"""
    relate_predicate(alg::RelateNG, predicate::TopologyPredicate, a, b)::Bool

Tests whether the topological relationship between geometries `a` and `b`
satisfies the given [`TopologyPredicate`](@ref). This is the core evaluation
entry point (the port of `RelateNG.evaluate(Geometry, TopologyPredicate)`,
via the static `RelateNG.relate(a, b, pred)`).

!!! note
    Predicates are mutable accumulators — pass a freshly constructed one
    (e.g. `pred_intersects()`) per evaluation.
"""
function relate_predicate(alg::RelateNG, predicate::TopologyPredicate, a, b)
    m = GeometryOpsCore.manifold(alg)
    geom_a = RelateGeometry(m, a; exact = alg.exact, boundary_rule = alg.boundary_rule)
    return evaluate!(alg, geom_a, b, predicate)
end

#==========================================================================
## Named-predicate methods

The ports of the JTS RelateNG static predicate overloads. These add
`RelateNG` algorithm methods to the existing GO predicate functions,
following the house `GO.f(alg::Algorithm, a, b)` idiom (cf.
`GO.intersects(GO.GEOS(), a, b)` in the LibGEOS extension). They are
opt-in: the two-argument forms `GO.intersects(a, b)` etc. keep
dispatching to the old engines (design D4).

`equals` maps to `pred_equalstopo`, i.e. *topological* equality (the
DE-9IM `T*F**FFF*` sense), which can differ from the structural equality
the two-argument `GO.equals` implements only in exotic cases (both
treat rotated/reversed rings and repeated points as equal).
==========================================================================#
intersects(alg::RelateNG, g1, g2) = relate_predicate(alg, pred_intersects(), g1, g2)
disjoint(alg::RelateNG, g1, g2)   = relate_predicate(alg, pred_disjoint(), g1, g2)
contains(alg::RelateNG, g1, g2)   = relate_predicate(alg, pred_contains(), g1, g2)
within(alg::RelateNG, g1, g2)     = relate_predicate(alg, pred_within(), g1, g2)
covers(alg::RelateNG, g1, g2)     = relate_predicate(alg, pred_covers(), g1, g2)
coveredby(alg::RelateNG, g1, g2)  = relate_predicate(alg, pred_coveredby(), g1, g2)
crosses(alg::RelateNG, g1, g2)    = relate_predicate(alg, pred_crosses(), g1, g2)
overlaps(alg::RelateNG, g1, g2)   = relate_predicate(alg, pred_overlaps(), g1, g2)
touches(alg::RelateNG, g1, g2)    = relate_predicate(alg, pred_touches(), g1, g2)
equals(alg::RelateNG, g1, g2)     = relate_predicate(alg, pred_equalstopo(), g1, g2)

#==========================================================================
## Evaluation (port of RelateNG.evaluate and helpers)
==========================================================================#

# Port of RelateNG.evaluate(Geometry b, TopologyPredicate predicate):
# the phased evaluation against a prebuilt A-side RelateGeometry. In
# prepared mode `prep` is the `PreparedRelate` carrying the cached A-side
# edges/index; otherwise `nothing`.
function evaluate!(alg::RelateNG, geom_a::RelateGeometry, b, predicate::TopologyPredicate,
        prep = nothing)
    #-- Java performs the envelope fast-exit before building the B
    #-- RelateGeometry, reading the envelope cached on the Geometry. Here the
    #-- RelateGeometry constructor is what caches the extents (one coordinate
    #-- pass either way), so B is built first and the check reads its extent.
    geom_b = RelateGeometry(geom_a.m, b; exact = geom_a.exact,
        boundary_rule = geom_a.boundary_rule)
    if !has_required_envelope_interaction(geom_a, get_extent(geom_b), predicate)
        return false
    end

    dim_a = get_dimension_real(geom_a)
    dim_b = get_dimension_real(geom_b)

    #-- check if predicate is determined by dimension or envelope
    init_dims!(predicate, dim_a, dim_b)
    is_known(predicate) && return finish_value!(predicate)

    init_bounds!(predicate, get_extent(geom_a), get_extent(geom_b))
    is_known(predicate) && return finish_value!(predicate)

    tc = TopologyComputer(predicate, geom_a, geom_b)

    #-- optimized P/P evaluation
    if dim_a == DIM_P && dim_b == DIM_P
        compute_pp!(tc, geom_a, geom_b)
        finish!(tc)
        return get_result(tc)
    end

    #-- test points against (potentially) indexed geometry first
    compute_at_points!(tc, geom_b, GEOM_B, geom_a)
    is_result_known(tc) && return get_result(tc)
    compute_at_points!(tc, geom_a, GEOM_A, geom_b)
    is_result_known(tc) && return get_result(tc)

    if has_edges(geom_a) && has_edges(geom_b)
        compute_at_edges!(alg, tc, geom_a, geom_b, prep)
    end

    #-- after all processing, set remaining unknown values in IM
    finish!(tc)
    return get_result(tc)
end

# Port of RelateNG.hasRequiredEnvelopeInteraction (private). `env_b` is the
# B extent computed directly from the raw geometry by the caller, as Java
# reads `b.getEnvelopeInternal()` before constructing the B RelateGeometry;
# an empty geometry yields a `nothing` (null) extent, for which `ext_covers`/
# `ext_intersects` return false, exactly like Java's null Envelope.
function has_required_envelope_interaction(geom_a::RelateGeometry, env_b, predicate::TopologyPredicate)
    is_interacts = false
    if require_covers(predicate, GEOM_A)
        if !ext_covers(get_extent(geom_a), env_b)
            return false
        end
        is_interacts = true
    elseif require_covers(predicate, GEOM_B)
        if !ext_covers(env_b, get_extent(geom_a))
            return false
        end
        is_interacts = true
    end
    if !is_interacts &&
            require_interaction(predicate) &&
            !ext_intersects(get_extent(geom_a), env_b)
        return false
    end
    return true
end

# Port of RelateNG.finishValue (private).
function finish_value!(predicate::TopologyPredicate)
    finish!(predicate)
    return predicate_value(predicate)
end

#=
Port of RelateNG.computePP (private): an optimized algorithm for evaluating
P/P cases. It tests one point set against the other.
=#
function compute_pp!(tc::TopologyComputer, geom_a::RelateGeometry, geom_b::RelateGeometry)
    pts_a = get_unique_points(geom_a)
    #TODO: only query points in interaction extent?
    pts_b = get_unique_points(geom_b)

    num_b_in_a = 0
    for pt_b in pts_b
        if pt_b in pts_a
            num_b_in_a += 1
            add_point_on_point_interior!(tc, pt_b)
        else
            add_point_on_point_exterior!(tc, GEOM_B, pt_b)
        end
        is_result_known(tc) && return nothing
    end
    #=
    If number of matched B points is less than size of A,
    there must be at least one A point in the exterior of B
    =#
    if num_b_in_a < length(pts_a)
        #TODO: determine actual exterior point?
        add_point_on_point_exterior!(tc, GEOM_A, nothing)
    end
    return nothing
end

# Port of RelateNG.computeAtPoints (private).
function compute_at_points!(tc::TopologyComputer, geom::RelateGeometry, is_a::Bool,
        geom_target::RelateGeometry)
    is_result_known_ = compute_points!(tc, geom, is_a, geom_target)
    is_result_known_ && return nothing

    #=
    Performance optimization: only check points against target
    if it has areas OR if the predicate requires checking for
    exterior interaction.
    In particular, this avoids testing line ends against lines
    for the intersects predicate (since these are checked
    during segment/segment intersection checking anyway).
    Checking points against areas is necessary, since the input
    linework is disjoint if one input lies wholly inside an area,
    so segment intersection checking is not sufficient.
    =#
    check_disjoint_points = has_dimension(geom_target, DIM_A) ||
        is_exterior_check_required(tc, is_a)
    check_disjoint_points || return nothing

    is_result_known_ = compute_line_ends!(tc, geom, is_a, geom_target)
    is_result_known_ && return nothing

    compute_area_vertex!(tc, geom, is_a, geom_target)
    return nothing
end

# Port of RelateNG.computePoints (private).
function compute_points!(tc::TopologyComputer, geom::RelateGeometry, is_a::Bool,
        geom_target::RelateGeometry)
    has_dimension(geom, DIM_P) || return false

    points = get_effective_points(geom)
    for point in points
        #TODO: exit when all possible target locations (E,I,B) have been found?
        GI.isempty(point) && continue

        pt = _to_kernel_point(geom.m, point)
        compute_point!(tc, is_a, pt, geom_target)
        is_result_known(tc) && return true
    end
    return false
end

# Port of RelateNG.computePoint (private).
function compute_point!(tc::TopologyComputer, is_a::Bool, pt, geom_target::RelateGeometry)
    loc_dim_target = locate_with_dim(geom_target, pt)
    loc_target = dimloc_location(loc_dim_target)
    dim_target = dimloc_dimension(loc_dim_target, get_dimension(tc, !is_a))
    add_point_on_geometry!(tc, is_a, loc_target, dim_target, pt)
    return nothing
end

# Port of RelateNG.computeLineEnds (private). The Java
# GeometryCollectionIterator + `instanceof LineString` filter becomes a
# recursive walk over atomic curve elements; the walk threads the
# `hasExteriorIntersection` flag and reports early exit through its first
# return value.
function compute_line_ends!(tc::TopologyComputer, geom::RelateGeometry, is_a::Bool,
        geom_target::RelateGeometry)
    has_dimension(geom, DIM_L) || return false
    done, _ = _compute_line_ends_walk!(tc, geom.geom, geom, is_a, geom_target, false)
    return done
end

function _compute_line_ends_walk!(tc::TopologyComputer, elem, geom::RelateGeometry,
        is_a::Bool, geom_target::RelateGeometry, has_exterior_intersection::Bool)
    trait = GI.trait(elem)
    if trait isa GI.AbstractGeometryCollectionTrait
        for g in GI.getgeom(elem)
            done, has_exterior_intersection = _compute_line_ends_walk!(
                tc, g, geom, is_a, geom_target, has_exterior_intersection)
            done && return (true, has_exterior_intersection)
        end
        return (false, has_exterior_intersection)
    end
    trait isa GI.AbstractCurveTrait || return (false, has_exterior_intersection)
    GI.isempty(elem) && return (false, has_exterior_intersection)

    #-- once an intersection with target exterior is recorded, skip further known-exterior points
    if has_exterior_intersection && _elem_env_disjoint(geom.m, elem, get_extent(geom_target))
        return (false, has_exterior_intersection)
    end

    e0 = _to_kernel_point(geom.m, GI.getpoint(elem, 1))
    has_exterior_intersection |= compute_line_end!(tc, geom, is_a, e0, geom_target)
    is_result_known(tc) && return (true, has_exterior_intersection)

    if !_line_is_closed(elem)
        e1 = _to_kernel_point(geom.m, GI.getpoint(elem, GI.npoint(elem)))
        has_exterior_intersection |= compute_line_end!(tc, geom, is_a, e1, geom_target)
        is_result_known(tc) && return (true, has_exterior_intersection)
    end
    #TODO: break when all possible locations have been found?
    return (false, has_exterior_intersection)
end

#=
Port of RelateNG.computeLineEnd (private): compute the topology of a line
endpoint. Also reports if the line end is in the exterior of the target
geometry, to optimize testing multiple exterior endpoints.
=#
function compute_line_end!(tc::TopologyComputer, geom::RelateGeometry, is_a::Bool, pt,
        geom_target::RelateGeometry)
    loc_dim_line_end = locate_line_end_with_dim(geom, pt)
    dim_line_end = dimloc_dimension(loc_dim_line_end, get_dimension(tc, is_a))
    #-- skip line ends which are in a GC area
    dim_line_end != DIM_L && return false
    loc_line_end = dimloc_location(loc_dim_line_end)

    loc_dim_target = locate_with_dim(geom_target, pt)
    loc_target = dimloc_location(loc_dim_target)
    dim_target = dimloc_dimension(loc_dim_target, get_dimension(tc, !is_a))
    add_line_end_on_geometry!(tc, is_a, loc_line_end, loc_target, dim_target, pt)
    return loc_target == LOC_EXTERIOR
end

# Port of RelateNG.computeAreaVertex(geom, isA, geomTarget, topoComputer)
# (private): the recursive walk over atomic polygon elements.
function compute_area_vertex!(tc::TopologyComputer, geom::RelateGeometry, is_a::Bool,
        geom_target::RelateGeometry)
    has_dimension(geom, DIM_A) || return false
    #-- evaluate for line and area targets only, since points are handled in the reverse direction
    get_dimension(geom_target) < DIM_L && return false
    done, _ = _compute_area_vertex_walk!(tc, geom.geom, geom, is_a, geom_target, false)
    return done
end

function _compute_area_vertex_walk!(tc::TopologyComputer, elem, geom::RelateGeometry,
        is_a::Bool, geom_target::RelateGeometry, has_exterior_intersection::Bool)
    trait = GI.trait(elem)
    if trait isa GI.AbstractGeometryCollectionTrait
        for g in GI.getgeom(elem)
            done, has_exterior_intersection = _compute_area_vertex_walk!(
                tc, g, geom, is_a, geom_target, has_exterior_intersection)
            done && return (true, has_exterior_intersection)
        end
        return (false, has_exterior_intersection)
    end
    trait isa GI.AbstractPolygonTrait || return (false, has_exterior_intersection)
    GI.isempty(elem) && return (false, has_exterior_intersection)

    #-- once an intersection with target exterior is recorded, skip further known-exterior points
    if has_exterior_intersection && _elem_env_disjoint(geom.m, elem, get_extent(geom_target))
        return (false, has_exterior_intersection)
    end

    has_exterior_intersection |=
        compute_area_vertex_on_ring!(tc, geom, is_a, GI.getexterior(elem), geom_target)
    is_result_known(tc) && return (true, has_exterior_intersection)

    for hole in GI.gethole(elem)
        has_exterior_intersection |=
            compute_area_vertex_on_ring!(tc, geom, is_a, hole, geom_target)
        is_result_known(tc) && return (true, has_exterior_intersection)
    end
    return (false, has_exterior_intersection)
end

# Port of RelateNG.computeAreaVertex(geom, isA, ring, geomTarget, topoComputer)
# (private).
function compute_area_vertex_on_ring!(tc::TopologyComputer, geom::RelateGeometry,
        is_a::Bool, ring, geom_target::RelateGeometry)
    #TODO: use extremal (highest) point to ensure one is on boundary of polygon cluster
    pt = _to_kernel_point(geom.m, GI.getpoint(ring, 1))

    loc_area = locate_area_vertex(geom, pt)
    loc_dim_target = locate_with_dim(geom_target, pt)
    loc_target = dimloc_location(loc_dim_target)
    dim_target = dimloc_dimension(loc_dim_target, get_dimension(tc, !is_a))
    add_area_vertex!(tc, is_a, loc_area, loc_target, dim_target, pt)
    return loc_target == LOC_EXTERIOR
end

# Port of RelateNG.computeAtEdges (private). The interaction envelope
# replaces Java's `Envelope.intersection` + `isNull` with
# `Extents.intersection`, which returns `nothing` for disjoint extents; an
# empty input's `nothing` extent also short-circuits here, so the extent
# filter passed to `extract_segment_strings` is never `nothing`
# (see the warning on that function).
function compute_at_edges!(alg::RelateNG, tc::TopologyComputer,
        geom_a::RelateGeometry, geom_b::RelateGeometry, prep = nothing)
    ext_a = get_extent(geom_a)
    ext_b = get_extent(geom_b)
    (ext_a === nothing || ext_b === nothing) && return nothing
    env_int = Extents.intersection(ext_a, ext_b)
    env_int === nothing && return nothing

    edges_b = extract_segment_strings(geom_b, GEOM_B, env_int)

    if is_self_noding_required(tc)
        #-- predicates requiring self-noding bypass the prepared cache: as in
        #-- Java, computeEdgesAll re-extracts the A edges per evaluation,
        #-- filtered by the interaction envelope (`prep` is not forwarded)
        compute_edges_all!(alg, tc, geom_a, edges_b, env_int)
    else
        compute_edges_mutual!(alg, tc, geom_a, edges_b, env_int, prep)
    end
    is_result_known(tc) && return nothing

    evaluate_nodes!(tc)
    return nothing
end

# Port of RelateNG.computeEdgesAll (private): the self-noding path. Java
# feeds A∪B through one EdgeSetIntersector (each unordered chain pair once,
# including self pairs); here the same pair set is phased as A×B, A×A, B×B.
function compute_edges_all!(alg::RelateNG, tc::TopologyComputer,
        geom_a::RelateGeometry, edges_b, env_int)
    #TODO: find a way to reuse prepared index?
    edges_a = extract_segment_strings(geom_a, GEOM_A, env_int)

    #-- mutual A×B pairs
    process_edge_intersections!(tc, edges_a, edges_b, alg.accelerator)
    is_result_known(tc) && return nothing
    #-- guarded A×A self pairs
    process_self_intersections!(tc, edges_a, alg.accelerator)
    is_result_known(tc) && return nothing
    #-- guarded B×B self pairs
    process_self_intersections!(tc, edges_b, alg.accelerator)
    return nothing
end

# Port of RelateNG.computeEdgesMutual (private). In prepared mode (`prep`
# a `PreparedRelate`) the cached A-side segment strings and prebuilt segment
# tree are reused — the port of Java's cached
# `MCIndexSegmentSetMutualIntersector`, which is built over A edges extracted
# with a *null* filter (`envExtract = geomA.isPrepared() ? null : envInt`):
# the cached strings are unfiltered so they serve any future B.
function compute_edges_mutual!(alg::RelateNG, tc::TopologyComputer,
        geom_a::RelateGeometry, edges_b, env_int, prep = nothing)
    if prep === nothing
        edges_a = extract_segment_strings(geom_a, GEOM_A, env_int)
        process_edge_intersections!(tc, edges_a, edges_b, alg.accelerator)
    else
        _process_prepared_edges!(tc, prep.segs_a, prep.edge_tree, edges_b)
    end
    return nothing
end

#==========================================================================
## Prepared mode

The port of the RelateNG.prepare entry points and the prepared-mode
branches.
==========================================================================#

"""
    PreparedRelate{ALG, RG, SS, T}

A prepared RelateNG instance for optimized repeated evaluation of
topological relationships against a single geometry `a` (the "prepared
mode" of JTS `RelateNG.prepare`). Holds:

- `alg`: the [`RelateNG`](@ref) algorithm configuration,
- `geom_a`: the A-side [`RelateGeometry`](@ref), constructed with
  `is_prepared = true` and with its lazy locator/unique-points caches
  forced,
- `segs_a`: the A segment strings, extracted once *without* an
  interaction-envelope filter so they serve any B geometry,
- `edge_tree`: the prebuilt segment index over `segs_a` (`_relate_edge_index`,
  edge_intersector.jl — the stand-in for Java's cached
  `MCIndexSegmentSetMutualIntersector`), or `nothing` below the accelerator
  size threshold (where the nested loop wins).

Construct with [`prepare`](@ref); evaluate with [`relate`](@ref) /
[`relate_predicate`](@ref).

!!! warning
    Not safe for concurrent use: self-noding evaluations mutate the held
    `RelateGeometry` (edge re-extraction, element-id counter). Use one
    `PreparedRelate` per thread.
"""
struct PreparedRelate{ALG <: RelateNG, RG <: RelateGeometry,
        SS <: AbstractVector{<:RelateSegmentString}, T <: Union{Nothing, RTree}}
    alg::ALG
    geom_a::RG
    segs_a::SS
    edge_tree::T
end

"""
    prepare(alg::RelateNG, a; validate = <manifold-dependent>)::PreparedRelate

`prepare` is the generic entry point for prepared-geometry optimizations in
GeometryOps; `RelateNG` is currently the only algorithm implementing it.

Creates a prepared relate instance to optimize the repeated evaluation of
relationships against the single geometry `a`.

Port of `RelateNG.prepare(Geometry)` (the algorithm's `boundary_rule` plays
the role of the `prepare(Geometry, BoundaryNodeRule)` overload). The A-side
`RelateGeometry` is constructed with `is_prepared = true`, and the lazy
caches that the Java instance accumulates across evaluations are forced
eagerly:

- the [`RelatePointLocator`](@ref) — the prepared flag selects the
  per-polygonal-element [`IndexedPointInAreaLocator`](@ref) caches, as in
  Java (see `point_locator.jl` / `indexed_point_in_area.jl`);
- the unique-point set, when `a` has effective dimension P (the only case
  the P/P fast path consults it);
- the A segment strings, extracted ONCE **without** an interaction-envelope
  filter (Java's prepared-mode `envExtract = null` in `computeEdgesMutual`)
  so the cache serves any future B, plus the prebuilt segment tree over
  them.

!!! note
    Predicates whose evaluation requires self-noding
    (`is_self_noding_required`) bypass the cached edges entirely: as in the
    Java prepared branch, `computeEdgesAll` re-extracts the A edges per
    evaluation, filtered by the A/B interaction envelope.

## Validation

`validate` controls a ring self-crossing check over the prepared geometry:
a self-join of the segment set (reusing the prepared edge index) that
detects PROPER crossings — transversal, interior to both edges — between
non-adjacent edges of the same polygonal element. Shared-endpoint adjacency
is excluded and vertex touches are NOT flagged: the scope is exactly the
crossing class that breaks the engine's containment parity, not full OGC
validity. On the first crossing found an `ArgumentError` is thrown naming
the ring and the edge pair (\"edge i crosses edge j\"); the documented
remedy is the [`CrossingEdgeSplit`](@ref) correction.

The default is manifold-dependent, and deliberately so:

- `Spherical` → `validate = true`. A planar-valid ring whose edges cross
  when reinterpreted as great-circle arcs is *undetectable by standard
  planar tooling* (planar validity checks pass it), and undetected it can
  invert containment globally — the figure-eight's lobes cancel the
  curvature the interior bootstrap reads, so every query lands on the
  wrong side (Natural Earth 110m Sudan is a real instance). The check is
  a small fraction of the ~100 ms spherical prepare build.
- `Planar` → `validate = false`. Planar invalidity of the same class is
  visible to ordinary planar tools, and the planar engine degrades
  gracefully under even-odd ray-crossing parity instead of inverting;
  JTS/GEOS prepared geometries do not validate either. A planar prepare
  costs ~600 µs, which a validation join would dominate, destroying the
  build-cost amortization. Pass `validate = true` to opt in.
"""
function prepare(alg::RelateNG, a;
        validate::Bool = _prepare_validate_default(GeometryOpsCore.manifold(alg)))
    m = GeometryOpsCore.manifold(alg)
    geom_a = RelateGeometry(m, a; exact = alg.exact, is_prepared = true,
        boundary_rule = alg.boundary_rule)
    #-- cached A edges: extracted once, unfiltered (Java's null envExtract).
    #-- Extracted (and validated) before the locator build below, so an
    #-- invalid ring fails fast instead of after the expensive locator pass.
    segs_a = extract_segment_strings(geom_a, GEOM_A, nothing)
    edge_tree = _build_prepared_edge_index(m, alg.accelerator, segs_a)
    validate && _validate_ring_crossings(m, alg.exact, segs_a, edge_tree)
    #-- force the lazy caches that repeated evaluations reuse
    _get_locator(geom_a)
    get_dimension_real(geom_a) == DIM_P && get_unique_points(geom_a)
    return PreparedRelate(alg, geom_a, segs_a, edge_tree)
end

# The manifold-dependent `validate` default of `prepare` (see its docstring
# for the rationale): only `Spherical`, where this invalidity class is both
# invisible to planar tooling and containment-inverting, validates by
# default.
_prepare_validate_default(::Spherical) = true
_prepare_validate_default(::Manifold) = false

#=
The validation join: enumerate extent-interacting segment pairs within the
A segment strings — through the prepared edge tree when one was built
(`dual_depth_first_search` of the tree against itself, exactly like the
self-noding tree path), otherwise a nested loop with the same per-pair
extent prune — and throw on the first PROPER crossing between non-adjacent
edges of the same polygonal element. Only ring edges are checked
(`dim == DIM_A` on both strings, same element id): lines may legitimately
self-cross, and rings of different elements (e.g. two polygons of a
MultiPolygon) may overlap without breaking per-element containment parity.
Shared-endpoint pairs are skipped by coordinate equality before any
predicate runs: adjacency and vertex touches are out of scope (see the
`prepare` docstring), and an exactly-zero orient would otherwise force the
adaptive exact stage on every adjacent pair.
=#
function _validate_ring_crossings(m::Manifold, exact, segs_a, edge_tree)
    if edge_tree === nothing
        _validate_ring_crossings_nested(m, exact, segs_a)
    else
        SpatialTreeInterface.dual_depth_first_search(Extents.intersects, edge_tree, edge_tree) do ia, ib
            ia < ib || return nothing
            (sa, ka) = edge_tree.data[ia]
            (sb, kb) = edge_tree.data[ib]
            _check_ring_crossing(m, exact, segs_a[sa], ka, segs_a[sb], kb)
            return nothing
        end
    end
    return nothing
end

function _validate_ring_crossings_nested(m::Manifold, exact, segs_a)
    for si in eachindex(segs_a)
        ssa = segs_a[si]
        ssa.dim == DIM_A || continue
        for sj in si:lastindex(segs_a)
            ssb = segs_a[sj]
            (ssb.dim == DIM_A && ssb.id == ssa.id) || continue
            for ka in 1:(length(ssa.pts) - 1)
                a0 = ssa.pts[ka]
                a1 = ssa.pts[ka + 1]
                kb0 = si == sj ? ka + 1 : 1
                for kb in kb0:(length(ssb.pts) - 1)
                    _segment_envs_disjoint(m, a0, a1, ssb.pts[kb], ssb.pts[kb + 1]) && continue
                    _check_ring_crossing(m, exact, ssa, ka, ssb, kb)
                end
            end
        end
    end
    return nothing
end

# One candidate pair of the validation join: filter (same-element rings,
# no shared endpoint), confirm (`_edges_cross_properly`), throw.
function _check_ring_crossing(m::Manifold, exact,
        ssa::RelateSegmentString, ka::Integer, ssb::RelateSegmentString, kb::Integer)
    (ssa.dim == DIM_A && ssb.dim == DIM_A && ssa.id == ssb.id) || return nothing
    a0 = ssa.pts[ka]
    a1 = ssa.pts[ka + 1]
    b0 = ssb.pts[kb]
    b1 = ssb.pts[kb + 1]
    #-- kernel-point `==` (all coordinates): adjacency and vertex touches
    #-- are out of scope, and skipping them here keeps the exactly-zero
    #-- orients of shared endpoints out of the adaptive exact stage
    (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1) && return nothing
    _edges_cross_properly(m, a0, a1, b0, b1; exact) &&
        _throw_ring_crossing(m, ssa, ka, ssb, kb)
    return nothing
end

@noinline function _throw_ring_crossing(m::Manifold, ssa, ka, ssb, kb)
    ia = _input_edge_index(m, ssa, ka)
    ib = _input_edge_index(m, ssb, kb)
    ssa === ssb && ib < ia && ((ia, ib) = (ib, ia))
    ring_desc(ss) = ss.ring_id == 0 ? "the shell" : "hole $(ss.ring_id)"
    place = ssa === ssb ?
        "edge $ia crosses edge $ib in $(ring_desc(ssa)) of polygonal element $(ssa.id)" :
        "edge $ia of $(ring_desc(ssa)) crosses edge $ib of $(ring_desc(ssb)) in polygonal element $(ssa.id)"
    throw(ArgumentError(
        "prepare: geometry is invalid on the `$(nameof(typeof(m)))` manifold: $place " *
        "(a proper crossing between non-adjacent ring edges). Left undetected, such a " *
        "crossing can invert containment globally; repair the geometry first with the " *
        "`CrossingEdgeSplit` correction (it splits each ring at its crossing points " *
        "into separate loops), or pass `validate = false` to skip this check"))
end

#=
Input-order edge index for the error message: segment strings store ring
vertices in engine orientation (shells and holes reversed as needed, repeated
points removed), so a stored segment index does not generally match an edge
of the ring as the user wrote it. Re-derive the input index by locating the
crossing edge's kernel endpoints among the input ring's consecutive vertex
pairs, in either direction. Falls back to the stored index if the ring or
pair cannot be found (it always should be). Error-path only.
=#
function _input_edge_index(m::Manifold, ss::RelateSegmentString, k::Integer)
    ring = _extracted_ring(ss)
    ring === nothing && return Int(k)
    e0 = ss.pts[k]
    e1 = ss.pts[k + 1]
    n = GI.npoint(ring)
    prev = _to_kernel_point(m, GI.getpoint(ring, 1))
    for i in 2:n
        cur = _to_kernel_point(m, GI.getpoint(ring, i))
        ((prev == e0 && cur == e1) || (prev == e1 && cur == e0)) && return i - 1
        prev = cur
    end
    return Int(k)
end

# The input ring a ring segment string was extracted from: the `ring_id`-th
# ring of the `id`-th non-empty atomic element of the input geometry.
function _extracted_ring(ss::RelateSegmentString)
    elem, _ = _nth_atomic_element(ss.input_geom.geom, Int(ss.id))
    (elem === nothing || !(GI.trait(elem) isa GI.AbstractPolygonTrait)) && return nothing
    ss.ring_id == 0 && return GI.getexterior(elem)
    for (i, hole) in enumerate(GI.gethole(elem))
        i == ss.ring_id && return hole
    end
    return nothing
end

# The `remaining`-th non-empty atomic element of `geom` in extraction order —
# the walk of `_extract_segment_strings!`, whose `element_id` counter ticks
# once per non-empty atomic element (collections, including Multi* types,
# recurse). Returns `(element_or_nothing, remaining_after)`.
function _nth_atomic_element(geom, remaining::Int)
    if GI.trait(geom) isa GI.AbstractGeometryCollectionTrait
        for g in GI.getgeom(geom)
            elem, remaining = _nth_atomic_element(g, remaining)
            elem === nothing || return (elem, remaining)
        end
        return (nothing, remaining)
    end
    GI.isempty(geom) && return (nothing, remaining)
    remaining -= 1
    return (remaining == 0 ? geom : nothing, remaining)
end

"""
    relate(p::PreparedRelate, b)::DE9IM

Computes the DE-9IM matrix for the topological relationship of the prepared
geometry to `b`. Port of the instance method `RelateNG.evaluate(Geometry)`.
"""
function relate(p::PreparedRelate, b)
    pred = RelateMatrixPredicate()
    relate_predicate(p, pred, b)
    return result_im(pred)
end

"""
    relate(p::PreparedRelate, b, im_pattern::AbstractString)::Bool

Tests whether the topological relationship of the prepared geometry to `b`
matches the DE-9IM pattern. Port of `RelateNG.evaluate(Geometry, String)`.
"""
relate(p::PreparedRelate, b, im_pattern::AbstractString) =
    relate_predicate(p, pred_matches(im_pattern), b)

"""
    relate_predicate(p::PreparedRelate, predicate::TopologyPredicate, b)::Bool

Tests whether the topological relationship of the prepared geometry to `b`
satisfies the predicate. Port of the instance method
`RelateNG.evaluate(Geometry, TopologyPredicate)` in prepared mode.
"""
relate_predicate(p::PreparedRelate, predicate::TopologyPredicate, b) =
    evaluate!(p.alg, p.geom_a, b, predicate, p)

# Whether to prebuild the A-side segment tree, mirroring the dispatch of
# `process_edge_intersections!` + `_select_edge_set_accelerator`: an explicit
# `NestedLoop` accelerator never uses a tree; `AutoAccelerator` uses one on
# the manifolds with a segment-extent kernel (`Planar`/`Spherical`) above
# the clipping size threshold only (B is unknown at prepare time, so the
# decision is made on A's segment count alone); any other explicit
# accelerator always takes the tree path.
_build_prepared_edge_index(m::Manifold, ::IntersectionAccelerator, segs_a) =
    _relate_edge_index(m, segs_a)
_build_prepared_edge_index(::Manifold, ::NestedLoop, segs_a) = nothing
_build_prepared_edge_index(::Manifold, ::AutoAccelerator, segs_a) = nothing
function _build_prepared_edge_index(m::Union{Planar, Spherical}, ::AutoAccelerator, segs_a)
    _total_segment_count(segs_a) >= GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS ||
        return nothing
    return _relate_edge_index(m, segs_a)
end

# The prepared counterpart of the mutual-pair enumeration: no prebuilt tree
# means the cached strings go through the plain nested loop.
_process_prepared_edges!(tc::TopologyComputer, segs_a, ::Nothing, edges_b) =
    process_edge_intersections!(tc, segs_a, edges_b, NestedLoop())

# Tree path: the prebuilt A tree is dual-traversed against a per-call tree
# over B's (envelope-filtered) segment extents — cf. the unprepared tree path
# in `process_edge_intersections!`, which builds both trees per call.
function _process_prepared_edges!(tc::TopologyComputer, segs_a,
        tree_a::RTree, edges_b;
        m::Manifold = _manifold(tc), exact = _exact(tc))
    tree_b = _relate_edge_index(m, edges_b)
    tree_b === nothing && return nothing
    SpatialTreeInterface.dual_depth_first_search(Extents.intersects, tree_a, tree_b) do ia, ib
        (sa, ka) = tree_a.data[ia]
        (sb, kb) = tree_b.data[ib]
        process_intersections!(tc, segs_a[sa], ka, edges_b[sb], kb; m, exact)
        #-- the Java noder's isDone() early-exit hook
        is_result_known(tc) && return Action(:full_return, nothing)
        return nothing
    end
    return nothing
end

#==========================================================================
## Small geometry helpers
==========================================================================#

# Java `elem.getEnvelopeInternal().disjoint(geomTarget.getEnvelope())`. A
# null (empty-geometry) target extent intersects nothing, hence is disjoint
# (only reachable when the target is empty but the predicate still requires
# exterior checks).
_elem_env_disjoint(m::Manifold, elem, target_ext) =
    target_ext === nothing || !Extents.intersects(rk_interaction_bounds(m, elem), target_ext)

# Java LineString.isClosed: false for an empty line, otherwise exact 2D
# coordinate equality of the endpoints.
function _line_is_closed(line)
    n = GI.npoint(line)
    n == 0 && return false
    return _equals2(_node_point(GI.getpoint(line, 1)), _node_point(GI.getpoint(line, n)))
end
