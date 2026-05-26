using Test
import GeometryOps as GO
import GeoInterface as GI
import GeoInterface.Extents: Extents

@testset "RelatePointLocator point and line locations" begin
    point_locator = GO.RelatePointLocator(GI.Point(1.0, 1.0))
    point_hit = GO.relate_locate_with_dim(point_locator, (1.0, 1.0))
    point_miss = GO.relate_locate_with_dim(point_locator, (2.0, 2.0))

    @test point_hit == GO.DimensionLocation(GO.dim_point, GO.loc_interior)
    @test point_miss == GO.RELATE_EXTERIOR
    @test GO.relate_locate(point_locator, (1.0, 1.0)) == GO.loc_interior

    line_locator = GO.RelatePointLocator(GI.LineString([(0.0, 0.0), (2.0, 0.0)]))
    @test GO.relate_locate_with_dim(line_locator, (1.0, 0.0)) ==
        GO.DimensionLocation(GO.dim_line, GO.loc_interior)
    @test GO.relate_locate_with_dim(line_locator, (0.0, 0.0)) ==
        GO.DimensionLocation(GO.dim_line, GO.loc_boundary)
    @test GO.relate_locate_with_dim(line_locator, (0.0, 1.0)) == GO.RELATE_EXTERIOR

    closed_line = GI.LineString([(0.0, 0.0), (1.0, 0.0), (0.0, 0.0)])
    mod2_locator = GO.RelatePointLocator(closed_line)
    endpoint_locator = GO.RelatePointLocator(
        closed_line;
        boundary_node_rule = GO.EndpointBoundaryNodeRule(),
    )
    @test GO.relate_locate_with_dim(mod2_locator, (0.0, 0.0)) ==
        GO.DimensionLocation(GO.dim_line, GO.loc_interior)
    @test GO.relate_locate_with_dim(endpoint_locator, (0.0, 0.0)) ==
        GO.DimensionLocation(GO.dim_line, GO.loc_boundary)
end

@testset "RelatePointLocator polygon and collection precedence" begin
    polygon = GI.Polygon([[
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.0, 0.0),
    ]])
    polygon_locator = GO.RelatePointLocator(polygon)

    @test GO.relate_locate_with_dim(polygon_locator, (1.0, 1.0)) ==
        GO.DimensionLocation(GO.dim_area, GO.loc_interior)
    @test GO.relate_locate_with_dim(polygon_locator, (0.0, 1.0)) ==
        GO.DimensionLocation(GO.dim_area, GO.loc_boundary)
    @test GO.relate_locate_with_dim(polygon_locator, (3.0, 3.0)) == GO.RELATE_EXTERIOR
    @test GO.relate_locate_node_with_dim(polygon_locator, (0.0, 0.0), polygon) ==
        GO.DimensionLocation(GO.dim_area, GO.loc_boundary)

    collection = GI.GeometryCollection([
        GI.Point(1.0, 1.0),
        GI.LineString([(0.0, 1.0), (2.0, 1.0)]),
        polygon,
    ])
    collection_locator = GO.RelatePointLocator(collection)

    @test GO.relate_locate_with_dim(collection_locator, (1.0, 1.0)) ==
        GO.DimensionLocation(GO.dim_area, GO.loc_interior)
    @test GO.relate_locate_with_dim(collection_locator, (1.0, 1.0); is_node = true) ==
        GO.DimensionLocation(GO.dim_area, GO.loc_interior)

    line_point_collection = GI.GeometryCollection([
        GI.Point(1.0, 0.0),
        GI.LineString([(0.0, 0.0), (2.0, 0.0)]),
    ])
    line_point_locator = GO.RelatePointLocator(line_point_collection)
    @test GO.relate_locate_with_dim(line_point_locator, (1.0, 0.0)) ==
        GO.DimensionLocation(GO.dim_line, GO.loc_interior)
end

@testset "RelatePointLocator polygonal union semantics" begin
    left = GI.Polygon([[
        (0.0, 0.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (0.0, 0.0),
    ]])
    right = GI.Polygon([[
        (1.0, 0.0),
        (2.0, 0.0),
        (2.0, 1.0),
        (1.0, 1.0),
        (1.0, 0.0),
    ]])
    multipolygon_locator = GO.RelatePointLocator(GI.MultiPolygon([left, right]))

    @test GO.relate_locate_with_dim(multipolygon_locator, (1.0, 0.5)) ==
        GO.DimensionLocation(GO.dim_area, GO.loc_interior)
    @test GO.relate_locate_with_dim(multipolygon_locator, (0.0, 0.5)) ==
        GO.DimensionLocation(GO.dim_area, GO.loc_boundary)

    corner = GI.Polygon([[
        (1.0, 1.0),
        (2.0, 1.0),
        (2.0, 2.0),
        (1.0, 2.0),
        (1.0, 1.0),
    ]])
    corner_touch_locator = GO.RelatePointLocator(GI.MultiPolygon([left, corner]))
    @test GO.relate_locate_with_dim(corner_touch_locator, (1.0, 1.0)) ==
        GO.DimensionLocation(GO.dim_area, GO.loc_boundary)

    area_line_collection = GI.GeometryCollection([
        GI.LineString([(0.5, 0.5), (1.5, 0.5)]),
        left,
    ])
    area_line_locator = GO.RelatePointLocator(area_line_collection)
    @test GO.relate_locate_line_end_with_dim(area_line_locator, (0.5, 0.5)) ==
        GO.DimensionLocation(GO.dim_area, GO.loc_interior)
end

@testset "RelateGeometry metadata and lazy caches" begin
    zero_line = GI.LineString([(2.0, 2.0), (2.0, 2.0), (2.0, 2.0)])
    zero_relate = GO.RelateGeometry(zero_line)

    @test zero_relate.dimension == GO.dim_line
    @test zero_relate.has_lines
    @test !zero_relate.has_points
    @test zero_relate.all_linework_zero_length
    @test GO.relate_dimension_real(zero_relate) == GO.dim_point
    @test GO.relate_has_dimension(zero_relate, GO.dim_line)
    @test GO.relate_has_edges(zero_relate)
    @test !GO.relate_has_area_and_line(zero_relate)
    @test GO.relate_is_self_noding_required(zero_relate)

    line = GI.LineString([(0.0, 0.0), (1.0, 0.0)])
    polygon = GI.Polygon([[
        (0.0, -1.0),
        (2.0, -1.0),
        (2.0, 1.0),
        (0.0, 1.0),
        (0.0, -1.0),
    ]])
    collection = GI.GeometryCollection([
        GI.Point(0.5, 0.0),
        GI.Point(3.0, 3.0),
        line,
        polygon,
    ])
    relate = GO.RelateGeometry(collection; prepared = true)

    @test relate.prepared
    @test relate.dimension == GO.dim_area
    @test GO.relate_dimension_real(relate) == GO.dim_area
    @test relate.has_points
    @test relate.has_lines
    @test relate.has_areas
    @test GO.relate_has_area_and_line(relate)
    @test GO.relate_has_edges(relate)
    @test GO.relate_is_self_noding_required(relate)

    locator = GO.relate_point_locator(relate)
    @test locator === GO.relate_point_locator(relate)
    @test GO.relate_locate_with_dim(relate, (0.5, 0.0)) ==
        GO.DimensionLocation(GO.dim_area, GO.loc_interior)
    @test GO.relate_locate_with_dim(relate, (3.0, 3.0)) ==
        GO.DimensionLocation(GO.dim_point, GO.loc_interior)

    unique_points = GO.relate_unique_points(relate)
    @test unique_points === GO.relate_unique_points(relate)
    @test unique_points == Set([(0.5, 0.0), (3.0, 3.0)])

    effective_points = GO.relate_effective_points(relate)
    @test getproperty.(effective_points, :point) == [(3.0, 3.0)]
end

@testset "RelateGeometry segment strings and prepared index" begin
    polygon = GI.Polygon([[
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 4.0),
        (0.0, 4.0),
        (0.0, 0.0),
    ], [
        (1.0, 1.0),
        (1.0, 3.0),
        (3.0, 3.0),
        (3.0, 1.0),
        (1.0, 1.0),
    ]])
    relate = GO.RelateGeometry(polygon)
    segments = GO.relate_segment_strings(relate; input_side = GO.input_b)

    @test segments === GO.relate_segment_strings(relate; input_side = GO.input_b)
    @test length(segments) == 2
    @test all(segment -> segment.source.input_side == GO.input_b, segments)
    @test segments[1].source.ring_role == GO.ring_shell
    @test segments[2].source.ring_role == GO.ring_hole
    @test GO.ng_is_clockwise(segments[1].points)
    @test !GO.ng_is_clockwise(segments[2].points)

    prepared_index = GO.relate_prepared_edge_index(relate)
    @test prepared_index === GO.relate_prepared_edge_index(relate)
    @test !isnothing(prepared_index)
    @test length(prepared_index.lines) == 8
    @test Extents.intersects(Extents.extent(prepared_index.index), GI.extent(prepared_index.lines[1]))
end

@testset "Relate topology interaction predicates" begin
    a = Extents.Extent(X = (0.0, 1.0), Y = (0.0, 1.0))
    b = Extents.Extent(X = (2.0, 3.0), Y = (2.0, 3.0))

    intersects_pred = GO.relate_intersects_predicate()
    @test GO.predicate_name(intersects_pred) == :intersects
    @test !GO.require_self_noding(intersects_pred)
    @test GO.require_interaction(intersects_pred)
    @test !GO.require_exterior_check(intersects_pred, GO.input_a)
    GO.relate_init_extents!(intersects_pred, a, b)
    @test GO.predicate_is_known(intersects_pred)
    @test !GO.predicate_value(intersects_pred)

    disjoint_pred = GO.relate_disjoint_predicate()
    @test !GO.require_self_noding(disjoint_pred)
    @test !GO.require_interaction(disjoint_pred)
    GO.relate_init_extents!(disjoint_pred, a, b)
    @test GO.predicate_is_known(disjoint_pred)
    @test GO.predicate_value(disjoint_pred)

    intersects_update = GO.relate_intersects_predicate()
    GO.relate_update_dimension!(
        intersects_update,
        GO.loc_boundary,
        GO.loc_interior,
        GO.dim_point,
    )
    @test GO.predicate_is_known(intersects_update)
    @test GO.predicate_value(intersects_update)
end

@testset "Relate topology IM predicates" begin
    contains_pred = GO.relate_contains_predicate()
    @test GO.require_covers(contains_pred, GO.input_a)
    @test !GO.require_covers(contains_pred, GO.input_b)
    @test !GO.require_exterior_check(contains_pred, GO.input_a)
    @test GO.require_exterior_check(contains_pred, GO.input_b)

    GO.relate_init_dimensions!(contains_pred, GO.dim_area, GO.dim_point)
    GO.relate_update_dimension!(contains_pred, GO.loc_interior, GO.loc_interior, GO.dim_point)
    GO.relate_finish!(contains_pred)
    @test GO.predicate_is_known(contains_pred)
    @test GO.predicate_value(contains_pred)

    bad_contains = GO.relate_contains_predicate()
    GO.relate_init_dimensions!(bad_contains, GO.dim_point, GO.dim_area)
    @test GO.predicate_is_known(bad_contains)
    @test !GO.predicate_value(bad_contains)

    covers_pred = GO.relate_covers_predicate()
    GO.relate_init_dimensions!(covers_pred, GO.dim_area, GO.dim_line)
    GO.relate_update_dimension!(covers_pred, GO.loc_boundary, GO.loc_interior, GO.dim_line)
    GO.relate_finish!(covers_pred)
    @test GO.predicate_value(covers_pred)

    touches_pred = GO.relate_touches_predicate()
    GO.relate_init_dimensions!(touches_pred, GO.dim_line, GO.dim_line)
    GO.relate_update_dimension!(touches_pred, GO.loc_boundary, GO.loc_interior, GO.dim_point)
    GO.relate_finish!(touches_pred)
    @test GO.predicate_value(touches_pred)

    point_touches = GO.relate_touches_predicate()
    GO.relate_init_dimensions!(point_touches, GO.dim_point, GO.dim_point)
    @test GO.predicate_is_known(point_touches)
    @test !GO.predicate_value(point_touches)
end

@testset "Relate topology pattern and matrix predicates" begin
    pattern_pred = GO.relate_matches_predicate("T********")
    @test GO.predicate_name(pattern_pred) == :matches
    @test GO.require_interaction(pattern_pred)
    GO.relate_update_dimension!(
        pattern_pred,
        GO.loc_interior,
        GO.loc_interior,
        GO.dim_point,
    )
    GO.relate_finish!(pattern_pred)
    @test GO.predicate_value(pattern_pred)

    false_pattern = GO.relate_matches_predicate("F********")
    GO.relate_update_dimension!(
        false_pattern,
        GO.loc_interior,
        GO.loc_interior,
        GO.dim_point,
    )
    @test GO.predicate_is_known(false_pattern)
    @test !GO.predicate_value(false_pattern)

    matrix_pred = GO.relate_matrix_predicate()
    @test !GO.require_interaction(matrix_pred)
    @test !GO.predicate_is_known(matrix_pred)
    @test GO.predicate_matrix(matrix_pred)[GO.loc_exterior, GO.loc_exterior] == GO.dim_area
    GO.relate_update_dimension!(
        matrix_pred,
        GO.loc_interior,
        GO.loc_boundary,
        GO.dim_line,
    )
    @test GO.predicate_matrix(matrix_pred)[GO.loc_interior, GO.loc_boundary] == GO.dim_line
end
