# # RelateNG engine
#
# Port of JTS `RelateNG.java` — the evaluation engine that drives all the
# RelateNG machinery (Tasks 1–20) to compute the value of a topological
# predicate between two geometries, based on the DE-9IM model.
#
# Method order parallels the Java file (`evaluate`, `hasRequiredEnvelopeInteraction`,
# `finishValue`, `computePP`, `computeAtPoints`, `computePoints`, `computePoint`,
# `computeLineEnds`, `computeLineEnd`, `computeAreaVertex` ×2, `computeAtEdges`,
# `computeEdgesAll`, `computeEdgesMutual`), so this file diffs against its
# Java counterpart. Idiom changes:
#
# - The Java class holds `geomA` (for prepared mode); here the unprepared
#   entry points build both `RelateGeometry`s per call. Prepared mode is
#   Task 22 (`PreparedRelate`).
# - The algorithm configuration (manifold, accelerator, exactness flag,
#   boundary node rule) travels in the `RelateNG` algorithm struct, the
#   house `Algorithm{M}` idiom (cf. `FosterHormannClipping`).
# - Java's `computeEdgesAll` feeds the combined A∪B edge set through one
#   `EdgeSetIntersector` (every unordered chain pair once). Here that is
#   phased as A×B (`process_edge_intersections!`), then A×A and B×B
#   (`process_self_intersections!`) — the same pair set, possibly visited
#   in a different order, which only affects *when* a short-circuiting
#   predicate exits, never the final value (dimension updates are
#   monotone).
# - Java's null `Envelope` is `nothing` here: empty geometries have a
#   `nothing` extent, `ext_intersects`/`ext_covers` return `false` for it
#   (as Java's null envelope does), and `computeAtEdges` early-returns
#   before ever forwarding an extent filter that could be `nothing`.

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

Keyword arguments (all optional): `manifold` (default `Planar()`),
`accelerator` (default [`AutoAccelerator`](@ref)), `exact` (default
`True()`), `boundary_rule` (default [`Mod2Boundary`](@ref), the OGC SFS
rule).

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
# Entry points (the static RelateNG.relate(...) overloads)
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
# Evaluation (port of RelateNG.evaluate and helpers)
==========================================================================#

# Port of RelateNG.evaluate(Geometry b, TopologyPredicate predicate):
# the phased evaluation against a prebuilt A-side RelateGeometry.
function evaluate!(alg::RelateNG, geom_a::RelateGeometry, b, predicate::TopologyPredicate)
    #-- fast envelope checks
    if !has_required_envelope_interaction(geom_a, b, predicate)
        return false
    end

    geom_b = RelateGeometry(geom_a.m, b; exact = geom_a.exact,
        boundary_rule = geom_a.boundary_rule)

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
        compute_at_edges!(alg, tc, geom_a, geom_b)
    end

    #-- after all processing, set remaining unknown values in IM
    finish!(tc)
    return get_result(tc)
end

# Port of RelateNG.hasRequiredEnvelopeInteraction (private). The B extent is
# computed directly from the raw geometry, as Java reads
# `b.getEnvelopeInternal()` before constructing the B RelateGeometry; an
# empty geometry yields a `nothing` (null) extent, for which `ext_covers`/
# `ext_intersects` return false, exactly like Java's null Envelope.
function has_required_envelope_interaction(geom_a::RelateGeometry, b, predicate::TopologyPredicate)
    env_b = _relate_extent(geom_a.m, b)
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

        pt = _node_point(point)
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

    e0 = _node_point(GI.getpoint(elem, 1))
    has_exterior_intersection |= compute_line_end!(tc, geom, is_a, e0, geom_target)
    is_result_known(tc) && return (true, has_exterior_intersection)

    if !_line_is_closed(elem)
        e1 = _node_point(GI.getpoint(elem, GI.npoint(elem)))
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
    pt = _node_point(GI.getpoint(ring, 1))

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
        geom_a::RelateGeometry, geom_b::RelateGeometry)
    ext_a = get_extent(geom_a)
    ext_b = get_extent(geom_b)
    (ext_a === nothing || ext_b === nothing) && return nothing
    env_int = Extents.intersection(ext_a, ext_b)
    env_int === nothing && return nothing

    edges_b = extract_segment_strings(geom_b, GEOM_B, env_int)

    if is_self_noding_required(tc)
        compute_edges_all!(alg, tc, geom_a, edges_b, env_int)
    else
        compute_edges_mutual!(alg, tc, geom_a, edges_b, env_int)
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

# Port of RelateNG.computeEdgesMutual (private). (The Java prepared-mode
# index reuse — null extract filter + cached MCIndexSegmentSetMutualIntersector
# — is Task 22.)
function compute_edges_mutual!(alg::RelateNG, tc::TopologyComputer,
        geom_a::RelateGeometry, edges_b, env_int)
    edges_a = extract_segment_strings(geom_a, GEOM_A, env_int)
    process_edge_intersections!(tc, edges_a, edges_b, alg.accelerator)
    return nothing
end

#==========================================================================
# Small geometry helpers
==========================================================================#

# Java `elem.getEnvelopeInternal().disjoint(geomTarget.getEnvelope())`. A
# null (empty-geometry) target extent intersects nothing, hence is disjoint
# (only reachable when the target is empty but the predicate still requires
# exterior checks).
_elem_env_disjoint(m::Manifold, elem, target_ext) =
    target_ext === nothing || rk_bounds_disjoint(rk_interaction_bounds(m, elem), target_ext)

# Java LineString.isClosed: false for an empty line, otherwise exact 2D
# coordinate equality of the endpoints.
function _line_is_closed(line)
    n = GI.npoint(line)
    n == 0 && return false
    return _equals2(_node_point(GI.getpoint(line, 1)), _node_point(GI.getpoint(line, n)))
end
