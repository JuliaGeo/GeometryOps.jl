# # Executable RelateNG paths

"""
    relate_evaluate_predicate(alg, a, b, predicate)

Evaluate a RelateNG predicate for the currently ported RelateNG paths.
"""
function relate_evaluate_predicate(alg::RelateNG, a, b, predicate::TopologyPredicate)
    computer = RelateTopologyComputer(alg, predicate, a, b)
    relate_is_result_known(computer) && return _relate_finish_predicate!(computer)

    _relate_can_evaluate_current_paths(computer) ||
        throw(ArgumentError("RelateNG edge-vs-edge evaluation is not implemented yet."))

    dim_a = relate_dimension_real(computer.geom_a)
    dim_b = relate_dimension_real(computer.geom_b)
    if dim_a == dim_point && dim_b == dim_point
        relate_compute_point_point!(computer)
        return _relate_finish_predicate!(computer)
    end

    relate_compute_at_points!(computer, computer.geom_b, input_b, computer.geom_a)
    relate_is_result_known(computer) && return computer.predicate

    relate_compute_at_points!(computer, computer.geom_a, input_a, computer.geom_b)
    relate_is_result_known(computer) && return computer.predicate

    if relate_has_edges(computer.geom_a) && relate_has_edges(computer.geom_b)
        relate_compute_edge_intersections!(computer)
        relate_is_result_known(computer) && return computer.predicate
        relate_evaluate_nodes!(computer)
        relate_is_result_known(computer) && return computer.predicate
    end

    return _relate_finish_predicate!(computer)
end

function _relate_finish_predicate!(computer::RelateTopologyComputer)
    relate_finish!(computer)
    return computer.predicate
end

function _relate_can_evaluate_current_paths(computer::RelateTopologyComputer)
    dim_a = relate_dimension_real(computer.geom_a)
    dim_b = relate_dimension_real(computer.geom_b)
    dim_a == dim_false && return true
    dim_b == dim_false && return true
    dim_a == dim_point && return true
    dim_b == dim_point && return true
    computer.predicate isa RelateInteractionPredicate && return true
    !(relate_has_area_and_line(computer.geom_a) || relate_has_area_and_line(computer.geom_b)) &&
        return true
    return false
end

"""
    relate_compute_point_point!(computer)

Optimized point-set evaluation for point/multipoint operands.
"""
function relate_compute_point_point!(computer::RelateTopologyComputer)
    points_a = relate_unique_points(computer.geom_a)
    points_b = relate_unique_points(computer.geom_b)

    matched_b_in_a = 0
    for point_b in points_b
        if point_b in points_a
            matched_b_in_a += 1
            relate_add_point_on_point_interior!(computer, point_b)
        else
            relate_add_point_on_point_exterior!(computer, input_b, point_b)
        end
        relate_is_result_known(computer) && return computer
    end

    if matched_b_in_a < length(points_a)
        relate_add_point_on_point_exterior!(computer, input_a)
    end
    return computer
end

function relate_compute_at_points!(
    computer::RelateTopologyComputer,
    source::RelateGeometry,
    source_side::NGInputSide,
    target::RelateGeometry,
)
    relate_compute_points!(computer, source, source_side, target)
    relate_is_result_known(computer) && return computer

    check_disjoint_points =
        relate_has_dimension(target, dim_area) ||
        relate_is_exterior_check_required(computer, source_side)
    check_disjoint_points || return computer

    relate_compute_line_ends!(computer, source, source_side, target)
    relate_is_result_known(computer) && return computer

    relate_compute_area_vertices!(computer, source, source_side, target)
    return computer
end

function relate_compute_points!(
    computer::RelateTopologyComputer,
    source::RelateGeometry,
    source_side::NGInputSide,
    target::RelateGeometry,
)
    relate_has_dimension(source, dim_point) || return computer

    for point in relate_effective_points(source)
        relate_compute_point!(computer, source_side, point.point, target)
        relate_is_result_known(computer) && return computer
    end
    return computer
end

function relate_compute_point!(
    computer::RelateTopologyComputer,
    point_side::NGInputSide,
    point,
    target::RelateGeometry,
)
    target_location = relate_locate_with_dim(target, point)
    relate_add_point_on_geometry!(
        computer,
        point_side,
        target_location.location,
        _relate_dimension_or(target_location, target.dimension),
        point,
    )
    return computer
end

function relate_compute_line_ends!(
    computer::RelateTopologyComputer,
    source::RelateGeometry,
    source_side::NGInputSide,
    target::RelateGeometry,
)
    relate_has_dimension(source, dim_line) || return computer
    _relate_compute_line_ends!(
        computer,
        GI.trait(source.geom),
        source.geom,
        source,
        source_side,
        target,
    )
    return computer
end

_relate_compute_line_ends!(computer, ::GI.PointTrait, geom, source, source_side, target) = computer
_relate_compute_line_ends!(computer, ::GI.MultiPointTrait, geom, source, source_side, target) = computer
_relate_compute_line_ends!(computer, ::GI.PolygonTrait, geom, source, source_side, target) = computer
_relate_compute_line_ends!(computer, ::GI.MultiPolygonTrait, geom, source, source_side, target) = computer

function _relate_compute_line_ends!(
    computer,
    ::GI.AbstractCurveTrait,
    curve,
    source,
    source_side,
    target,
)
    GI.isempty(curve) && return computer
    GI.npoint(curve) == 0 && return computer

    relate_compute_line_end!(computer, source, source_side, GI.getpoint(curve, 1), target)
    relate_is_result_known(computer) && return computer

    if !_relate_curve_is_closed(curve)
        relate_compute_line_end!(
            computer,
            source,
            source_side,
            GI.getpoint(curve, GI.npoint(curve)),
            target,
        )
    end
    return computer
end

function _relate_compute_line_ends!(
    computer,
    ::GI.AbstractGeometryTrait,
    geom,
    source,
    source_side,
    target,
)
    GI.isempty(geom) && return computer
    for child in GI.getgeom(geom)
        _relate_compute_line_ends!(
            computer,
            GI.trait(child),
            child,
            source,
            source_side,
            target,
        )
        relate_is_result_known(computer) && return computer
    end
    return computer
end

function relate_compute_line_end!(
    computer::RelateTopologyComputer,
    source::RelateGeometry,
    source_side::NGInputSide,
    point,
    target::RelateGeometry,
)
    line_end = relate_locate_line_end_with_dim(source, point)
    line_end.dimension == dim_line || return computer

    target_location = relate_locate_with_dim(target, point)
    relate_add_line_end_on_geometry!(
        computer,
        source_side,
        line_end.location,
        target_location.location,
        _relate_dimension_or(target_location, target.dimension),
        point,
    )
    return computer
end

function relate_compute_area_vertices!(
    computer::RelateTopologyComputer,
    source::RelateGeometry,
    source_side::NGInputSide,
    target::RelateGeometry,
)
    relate_has_dimension(source, dim_area) || return computer
    dimension_value(target.dimension) < dimension_value(dim_line) && return computer
    _relate_compute_area_vertices!(
        computer,
        GI.trait(source.geom),
        source.geom,
        source,
        source_side,
        target,
    )
    return computer
end

_relate_compute_area_vertices!(computer, ::GI.PointTrait, geom, source, source_side, target) = computer
_relate_compute_area_vertices!(computer, ::GI.MultiPointTrait, geom, source, source_side, target) = computer
_relate_compute_area_vertices!(computer, ::GI.AbstractCurveTrait, geom, source, source_side, target) = computer

function _relate_compute_area_vertices!(
    computer,
    ::GI.PolygonTrait,
    polygon,
    source,
    source_side,
    target,
)
    GI.isempty(polygon) && return computer

    relate_compute_area_vertex!(
        computer,
        source,
        source_side,
        GI.getexterior(polygon),
        target,
    )
    relate_is_result_known(computer) && return computer

    for hole in GI.gethole(polygon)
        relate_compute_area_vertex!(computer, source, source_side, hole, target)
        relate_is_result_known(computer) && return computer
    end
    return computer
end

function _relate_compute_area_vertices!(
    computer,
    ::GI.MultiPolygonTrait,
    multipolygon,
    source,
    source_side,
    target,
)
    GI.isempty(multipolygon) && return computer
    for polygon in GI.getgeom(multipolygon)
        _relate_compute_area_vertices!(
            computer,
            GI.trait(polygon),
            polygon,
            source,
            source_side,
            target,
        )
        relate_is_result_known(computer) && return computer
    end
    return computer
end

function _relate_compute_area_vertices!(
    computer,
    ::GI.AbstractGeometryTrait,
    geom,
    source,
    source_side,
    target,
)
    GI.isempty(geom) && return computer
    for child in GI.getgeom(geom)
        _relate_compute_area_vertices!(
            computer,
            GI.trait(child),
            child,
            source,
            source_side,
            target,
        )
        relate_is_result_known(computer) && return computer
    end
    return computer
end

function relate_compute_area_vertex!(
    computer::RelateTopologyComputer,
    source::RelateGeometry,
    source_side::NGInputSide,
    ring,
    target::RelateGeometry,
)
    GI.isempty(ring) && return computer
    GI.npoint(ring) == 0 && return computer

    point = _tuple_point(GI.getpoint(ring, 1), Float64)
    area_location = relate_locate_area_vertex(source, point).location
    target_location = relate_locate_with_dim(target, point)
    relate_add_area_vertex!(
        computer,
        source_side,
        area_location,
        target_location.location,
        _relate_dimension_or(target_location, target.dimension),
        point,
    )
    return computer
end

function _relate_dimension_or(dimloc::DimensionLocation, fallback::TopologicalDimension)
    dimloc.location == loc_exterior && dimloc.dimension == dim_false && return fallback
    return dimloc.dimension
end

function _relate_curve_is_closed(curve)
    GI.npoint(curve) > 1 || return false
    return _tuple_point(GI.getpoint(curve, 1)) == _tuple_point(GI.getpoint(curve, GI.npoint(curve)))
end

relate(alg::RelateNG, a, b, predicate::TopologyPredicate) =
    predicate_value(relate_evaluate_predicate(alg, a, b, predicate))

relate(alg::RelateNG, a, b, pattern::AbstractString) =
    relate(alg, a, b, relate_matches_predicate(pattern))

relate_matrix(alg::RelateNG, a, b) =
    predicate_matrix(relate_evaluate_predicate(alg, a, b, relate_matrix_predicate()))

intersects(alg::RelateNG, a, b) = relate(alg, a, b, relate_intersects_predicate())
disjoint(alg::RelateNG, a, b) = relate(alg, a, b, relate_disjoint_predicate())
contains(alg::RelateNG, a, b) = relate(alg, a, b, relate_contains_predicate())
within(alg::RelateNG, a, b) = relate(alg, a, b, relate_within_predicate())
covers(alg::RelateNG, a, b) = relate(alg, a, b, relate_covers_predicate())
coveredby(alg::RelateNG, a, b) = relate(alg, a, b, relate_coveredby_predicate())
crosses(alg::RelateNG, a, b) = relate(alg, a, b, relate_crosses_predicate())
equals(alg::RelateNG, a, b) = relate(alg, a, b, relate_equals_predicate())
overlaps(alg::RelateNG, a, b) = relate(alg, a, b, relate_overlaps_predicate())
touches(alg::RelateNG, a, b) = relate(alg, a, b, relate_touches_predicate())
