using Test
import GeometryOps as GO
import GeoInterface as GI
import GeoInterface.Extents: Extents
using GeometryOpsTestHelpers

const _JTS_RELATE_MISC_FIXTURE_DIR = normpath(joinpath(
    @__DIR__,
    "..",
    "..",
    "..",
    "..",
    "jts",
    "modules",
    "tests",
    "src",
    "test",
    "resources",
    "testxml",
    "misc",
))

const _JTS_RELATE_GENERAL_FIXTURE_DIR = normpath(joinpath(
    @__DIR__,
    "..",
    "..",
    "..",
    "..",
    "jts",
    "modules",
    "tests",
    "src",
    "test",
    "resources",
    "testxml",
    "general",
))

const _JTS_RELATE_ROBUST_FIXTURE_DIR = normpath(joinpath(
    @__DIR__,
    "..",
    "..",
    "..",
    "..",
    "jts",
    "modules",
    "tests",
    "src",
    "test",
    "resources",
    "testxml",
    "robust",
))

_relateng_jts_misc_fixtures_available() =
    isfile(joinpath(_JTS_RELATE_MISC_FIXTURE_DIR, "TestRelateGC.xml"))

_relateng_jts_empty_fixtures_available() =
    isfile(joinpath(_JTS_RELATE_MISC_FIXTURE_DIR, "TestRelateEmpty.xml"))

_relateng_jts_general_fixtures_available() =
    isfile(joinpath(_JTS_RELATE_GENERAL_FIXTURE_DIR, "TestRelatePP.xml"))

_relateng_jts_robust_fixtures_available() =
    isfile(joinpath(_JTS_RELATE_ROBUST_FIXTURE_DIR, "TestRobustRelate.xml"))

function _relateng_fixture_value(alg::GO.RelateNG, op::JTSOperation)
    name = lowercase(op.name)
    a = _relateng_fixture_argument(op.arguments[1])
    b = _relateng_fixture_argument(op.arguments[2])
    if name == "relate"
        return GO.relate(alg, a, b, op.arguments[3])
    elseif name == "intersects"
        return GO.intersects(alg, a, b)
    elseif name == "contains"
        return GO.contains(alg, a, b)
    elseif name == "covers"
        return GO.covers(alg, a, b)
    elseif name == "coveredby"
        return GO.coveredby(alg, a, b)
    elseif name == "within"
        return GO.within(alg, a, b)
    elseif name == "touches"
        return GO.touches(alg, a, b)
    elseif name == "crosses"
        return GO.crosses(alg, a, b)
    elseif name == "overlaps"
        return GO.overlaps(alg, a, b)
    elseif name == "disjoint"
        return GO.disjoint(alg, a, b)
    elseif name == "equalstopo"
        return GO.equals(alg, a, b)
    end
    error("Unsupported RelateNG fixture operation: $(op.name)")
end

_relateng_fixture_argument(geom) = geom
_relateng_fixture_argument(::JTSEmptyGeometry) = GI.FeatureCollection(Any[])

@testset "RelatePointLocator point and line locations" begin
    point_locator = GO.RelatePointLocator(GI.Point(1.0, 1.0))
    point_hit = GO.relate_locate_with_dim(point_locator, (1.0, 1.0))
    point_miss = GO.relate_locate_with_dim(point_locator, (2.0, 2.0))

    @test point_hit == GO.DimensionLocation(GO.dim_point, GO.loc_interior)
    @test point_miss == GO.RELATE_EXTERIOR
    @test GO.relate_locate(point_locator, (1.0, 1.0)) == GO.loc_interior

    empty_locator = GO.RelatePointLocator(GI.FeatureCollection(Any[]))
    @test GO.relate_locate_with_dim(empty_locator, (0.0, 0.0)) == GO.RELATE_EXTERIOR

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
    @test length(prepared_index.records) == 8
    @test length(prepared_index.lines) == 8
    @test Extents.intersects(Extents.extent(prepared_index.index), GI.extent(prepared_index.lines[1]))
end

@testset "RelateNG prepared geometry queries" begin
    alg = GO.RelateNG(; prepared = true)
    polygon = GI.Polygon([[
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 4.0),
        (0.0, 4.0),
        (0.0, 0.0),
    ]])
    prepared_polygon = GO.relate_prepare(alg, polygon)

    @test prepared_polygon.prepared
    @test GO.relate_prepare(alg, prepared_polygon) === prepared_polygon

    unprepared_polygon = GO.RelateGeometry(polygon)
    @test GO.relate_prepare(unprepared_polygon) === unprepared_polygon
    @test unprepared_polygon.prepared

    locator = GO.relate_point_locator(prepared_polygon)
    segments = GO.relate_segment_strings(prepared_polygon; input_side = GO.input_a)
    edge_index = GO.relate_prepared_edge_index(prepared_polygon)

    crossing_line = GI.LineString([(-1.0, 2.0), (5.0, 2.0)])
    interior_line = GI.LineString([(1.0, 1.0), (3.0, 1.0)])

    @test GO.intersects(alg, prepared_polygon, crossing_line)
    @test GO.crosses(alg, prepared_polygon, crossing_line)
    @test GO.contains(alg, prepared_polygon, interior_line)
    @test GO.covers(alg, prepared_polygon, interior_line)

    @test GO.relate_point_locator(prepared_polygon) === locator
    @test GO.relate_segment_strings(prepared_polygon; input_side = GO.input_a) === segments
    @test GO.relate_prepared_edge_index(prepared_polygon) === edge_index

    endpoint_alg = GO.RelateNG(; boundary_node_rule = GO.EndpointBoundaryNodeRule())
    endpoint_prepared = GO.relate_prepare(endpoint_alg, polygon)
    @test_throws ArgumentError GO.intersects(alg, endpoint_prepared, interior_line)
end

@testset "RelateNG JTS XML conformance smoke" begin
    xml = """
    <run>
      <desc>RelateNG staged conformance smoke</desc>
      <case>
        <desc>PP equal</desc>
        <a>POINT (0 0)</a>
        <b>POINT (0 0)</b>
        <test><op name="relate" arg1="A" arg2="B" arg3="0FFFFFFF2">true</op></test>
        <test><op name="intersects" arg1="A" arg2="B">true</op></test>
        <test><op name="disjoint" arg1="A" arg2="B">false</op></test>
        <test><op name="equalsTopo" arg1="A" arg2="B">true</op></test>
      </case>
      <case>
        <desc>PL endpoint</desc>
        <a>POINT (0 0)</a>
        <b>LINESTRING (0 0, 2 0)</b>
        <test><op name="intersects" arg1="A" arg2="B">true</op></test>
        <test><op name="touches" arg1="A" arg2="B">true</op></test>
        <test><op name="within" arg1="A" arg2="B">false</op></test>
        <test><op name="coveredBy" arg1="A" arg2="B">true</op></test>
      </case>
      <case>
        <desc>PA interior</desc>
        <a>POINT (1 1)</a>
        <b>POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))</b>
        <test><op name="intersects" arg1="A" arg2="B">true</op></test>
        <test><op name="within" arg1="A" arg2="B">true</op></test>
        <test><op name="coveredBy" arg1="A" arg2="B">true</op></test>
        <test><op name="touches" arg1="A" arg2="B">false</op></test>
      </case>
      <case>
        <desc>LL crossing</desc>
        <a>LINESTRING (0 0, 2 2)</a>
        <b>LINESTRING (0 2, 2 0)</b>
        <test><op name="intersects" arg1="A" arg2="B">true</op></test>
        <test><op name="crosses" arg1="A" arg2="B">true</op></test>
        <test><op name="touches" arg1="A" arg2="B">false</op></test>
        <test><op name="overlaps" arg1="A" arg2="B">false</op></test>
      </case>
      <case>
        <desc>LA boundary line</desc>
        <a>LINESTRING (0 0, 2 0)</a>
        <b>POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))</b>
        <test><op name="intersects" arg1="A" arg2="B">true</op></test>
        <test><op name="touches" arg1="A" arg2="B">true</op></test>
        <test><op name="within" arg1="A" arg2="B">false</op></test>
        <test><op name="coveredBy" arg1="A" arg2="B">true</op></test>
      </case>
      <case>
        <desc>AA overlap</desc>
        <a>POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))</a>
        <b>POLYGON ((1 0, 3 0, 3 2, 1 2, 1 0))</b>
        <test><op name="intersects" arg1="A" arg2="B">true</op></test>
        <test><op name="overlaps" arg1="A" arg2="B">true</op></test>
        <test><op name="touches" arg1="A" arg2="B">false</op></test>
      </case>
      <case>
        <desc>AA adjacent</desc>
        <a>POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))</a>
        <b>POLYGON ((2 0, 4 0, 4 2, 2 2, 2 0))</b>
        <test><op name="intersects" arg1="A" arg2="B">true</op></test>
        <test><op name="touches" arg1="A" arg2="B">true</op></test>
        <test><op name="overlaps" arg1="A" arg2="B">false</op></test>
      </case>
    </run>
    """

    mktemp() do path, io
        write(io, xml)
        close(io)

        test_set = load_test_set(path)
        inventory = conformance_inventory(test_set)
        @test inventory[:point_point] == 1
        @test inventory[:point_line] == 1
        @test inventory[:point_area] == 1
        @test inventory[:line_line] == 1
        @test inventory[:line_area] == 1
        @test inventory[:area_area] == 2

        for alg in (GO.RelateNG(), GO.RelateNG(; prepared = true))
            for case in test_set.cases
                for op in case.operations
                    @test is_relate_operation(op)
                    @test is_runnable(op)
                    @test _relateng_fixture_value(alg, op) === op.expected
                end
            end
        end
    end
end

@testset "RelateNG JTS general fixtures" begin
    if !_relateng_jts_general_fixtures_available()
        @test_skip _relateng_jts_general_fixtures_available()
    else
        fixtures = (
            ("TestRelatePP.xml", 1:4),
            ("TestRelatePL.xml", 1:8),
            ("TestRelatePA.xml", 1:11),
            ("TestRelateAA.xml", 1:14),
            ("TestRelateLA.xml", 1:11),
            ("TestRelateLL.xml", [1:20; 22:25]),
        )

        matched_operations = 0
        for alg in (GO.RelateNG(), GO.RelateNG(; prepared = true))
            for (filename, case_indices) in fixtures
                test_set = load_test_set(joinpath(_JTS_RELATE_GENERAL_FIXTURE_DIR, filename))
                for case_index in case_indices
                    case = test_set.cases[case_index]
                    for op in case.operations
                        is_relate_operation(op) || continue
                        matched_operations += 1
                        @test _relateng_fixture_value(alg, op) === op.expected
                    end
                end
            end
        end
        @test matched_operations == 458
    end
end

@testset "RelateNG JTS empty fixtures" begin
    if !_relateng_jts_empty_fixtures_available()
        @test_skip _relateng_jts_empty_fixtures_available()
    else
        test_set = load_test_set(joinpath(_JTS_RELATE_MISC_FIXTURE_DIR, "TestRelateEmpty.xml"))
        matched_operations = 0
        for alg in (GO.RelateNG(), GO.RelateNG(; prepared = true))
            for case in test_set.cases
                for op in case.operations
                    is_relate_operation(op) || continue
                    matched_operations += 1
                    @test _relateng_fixture_value(alg, op) === op.expected
                end
            end
        end
        @test matched_operations == 1144
    end
end

@testset "RelateNG JTS robust fixtures" begin
    if !_relateng_jts_robust_fixtures_available()
        @test_skip _relateng_jts_robust_fixtures_available()
    else
        fixtures = (
            ("TestRobustRelate.xml", 1:1),
            ("TestRobustRelateFloat.xml", 2:2),
        )

        matched_operations = 0
        for alg in (GO.RelateNG(), GO.RelateNG(; prepared = true))
            for (filename, case_indices) in fixtures
                test_set = load_test_set(joinpath(_JTS_RELATE_ROBUST_FIXTURE_DIR, filename))
                for case_index in case_indices
                    case = test_set.cases[case_index]
                    for op in case.operations
                        is_relate_operation(op) || continue
                        matched_operations += 1
                        @test _relateng_fixture_value(alg, op) === op.expected
                    end
                end
            end
        end
        @test matched_operations == 4
    end
end

@testset "RelateNG JTS geometry collection fixtures" begin
    if !_relateng_jts_misc_fixtures_available()
        @test_skip _relateng_jts_misc_fixtures_available()
    else
        test_set = load_test_set(joinpath(_JTS_RELATE_MISC_FIXTURE_DIR, "TestRelateGC.xml"))
        matched_operations = 0
        for alg in (GO.RelateNG(), GO.RelateNG(; prepared = true))
            for case in test_set.cases
                for op in case.operations
                    matched_operations += 1
                    @test _relateng_fixture_value(alg, op) === op.expected
                end
            end
        end
        @test matched_operations == 656
    end
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

@testset "RelateTopologyComputer exterior seeding" begin
    point = GO.RelateGeometry(GI.Point(0.0, 0.0))
    line = GO.RelateGeometry(GI.LineString([(0.0, 0.0), (1.0, 0.0)]))
    area = GO.RelateGeometry(GI.Polygon([[
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.0, 0.0),
    ]]))
    empty = GO.RelateGeometry(GI.FeatureCollection(Any[]))

    point_line = GO.RelateTopologyComputer(GO.relate_matrix_predicate(), point, line)
    @test GO.predicate_matrix(point_line.predicate)[GO.loc_exterior, GO.loc_interior] ==
        GO.dim_line

    area_point = GO.RelateTopologyComputer(GO.relate_matrix_predicate(), area, point)
    area_point_matrix = GO.predicate_matrix(area_point.predicate)
    @test area_point_matrix[GO.loc_interior, GO.loc_exterior] == GO.dim_area
    @test area_point_matrix[GO.loc_boundary, GO.loc_exterior] == GO.dim_line

    line_area = GO.RelateTopologyComputer(GO.relate_matrix_predicate(), line, area)
    @test GO.predicate_matrix(line_area.predicate)[GO.loc_exterior, GO.loc_interior] ==
        GO.dim_area

    line_empty = GO.RelateTopologyComputer(GO.relate_matrix_predicate(), line, empty)
    line_empty_matrix = GO.predicate_matrix(line_empty.predicate)
    @test line_empty_matrix[GO.loc_boundary, GO.loc_exterior] == GO.dim_point
    @test line_empty_matrix[GO.loc_interior, GO.loc_exterior] == GO.dim_line

    area_empty = GO.RelateTopologyComputer(GO.relate_matrix_predicate(), empty, area)
    area_empty_matrix = GO.predicate_matrix(area_empty.predicate)
    @test area_empty_matrix[GO.loc_exterior, GO.loc_boundary] == GO.dim_line
    @test area_empty_matrix[GO.loc_exterior, GO.loc_interior] == GO.dim_area
end

@testset "RelateTopologyComputer events and short circuiting" begin
    point = GO.RelateGeometry(GI.Point(1.0, 1.0))
    area = GO.RelateGeometry(GI.Polygon([[
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.0, 0.0),
    ]]))
    line = GO.RelateGeometry(GI.LineString([(-1.0, 1.0), (3.0, 1.0)]))

    intersects_computer = GO.RelateTopologyComputer(
        GO.relate_intersects_predicate(),
        point,
        GO.RelateGeometry(GI.Point(1.0, 1.0)),
    )
    @test !GO.relate_is_result_known(intersects_computer)
    GO.relate_add_point_on_point_interior!(intersects_computer, (1.0, 1.0))
    @test GO.relate_is_result_known(intersects_computer)
    @test GO.relate_result(intersects_computer)

    point_area = GO.RelateTopologyComputer(GO.relate_matrix_predicate(), point, area)
    GO.relate_add_point_on_geometry!(
        point_area,
        GO.input_a,
        GO.DimensionLocation(GO.dim_area, GO.loc_interior),
        (1.0, 1.0),
    )
    point_area_matrix = GO.predicate_matrix(point_area.predicate)
    @test point_area_matrix[GO.loc_interior, GO.loc_interior] == GO.dim_point
    @test point_area_matrix[GO.loc_exterior, GO.loc_boundary] == GO.dim_line

    line_area = GO.RelateTopologyComputer(GO.relate_matrix_predicate(), line, area)
    GO.relate_add_line_end_on_geometry!(
        line_area,
        GO.input_a,
        GO.loc_boundary,
        GO.DimensionLocation(GO.dim_area, GO.loc_exterior),
        (-1.0, 1.0),
    )
    line_area_matrix = GO.predicate_matrix(line_area.predicate)
    @test line_area_matrix[GO.loc_boundary, GO.loc_exterior] == GO.dim_point
    @test line_area_matrix[GO.loc_interior, GO.loc_exterior] == GO.dim_line
    @test line_area_matrix[GO.loc_exterior, GO.loc_exterior] == GO.dim_area

    area_point = GO.RelateTopologyComputer(GO.relate_matrix_predicate(), area, point)
    GO.relate_add_area_vertex!(
        area_point,
        GO.input_a,
        GO.loc_boundary,
        GO.DimensionLocation(GO.dim_point, GO.loc_interior),
        (0.0, 0.0),
    )
    area_point_matrix = GO.predicate_matrix(area_point.predicate)
    @test area_point_matrix[GO.loc_boundary, GO.loc_interior] == GO.dim_point
    @test area_point_matrix[GO.loc_boundary, GO.loc_exterior] == GO.dim_line
    @test area_point_matrix[GO.loc_exterior, GO.loc_exterior] == GO.dim_area
end

@testset "RelateTopologyComputer self-noding and node sections" begin
    line_a_geom = GI.LineString([(0.0, 0.0), (2.0, 0.0)])
    line_b_geom = GI.LineString([(1.0, -1.0), (1.0, 1.0)])
    line_a = GO.RelateGeometry(line_a_geom)
    line_b = GO.RelateGeometry(line_b_geom)

    matrix_computer = GO.RelateTopologyComputer(GO.relate_matrix_predicate(), line_a, line_b)
    @test GO.relate_is_self_noding_required(matrix_computer)

    interaction_computer =
        GO.RelateTopologyComputer(GO.relate_intersects_predicate(), line_a, line_b)
    @test !GO.relate_is_self_noding_required(interaction_computer)

    a_segment = only(GO.relate_segment_strings(line_a; input_side = GO.input_a))
    b_segment = only(GO.relate_segment_strings(line_b; input_side = GO.input_b))
    section_a = GO.RelateNodeSection(a_segment, 1, (1.0, 0.0))
    section_b = GO.RelateNodeSection(b_segment, 1, (1.0, 0.0))

    @test !section_a.is_node_at_vertex
    @test section_a.previous_vertex == (0.0, 0.0)
    @test section_a.next_vertex == (2.0, 0.0)

    GO.relate_add_intersection!(matrix_computer, section_a, section_b)
    @test length(GO.relate_node_sections(matrix_computer, (1.0, 0.0))) == 2
    @test GO.predicate_matrix(matrix_computer.predicate)[GO.loc_interior, GO.loc_interior] ==
        GO.dim_point

    mixed_b = GO.RelateGeometry(GI.GeometryCollection([
        line_b_geom,
        GI.Polygon([[
            (10.0, 10.0),
            (11.0, 10.0),
            (11.0, 11.0),
            (10.0, 11.0),
            (10.0, 10.0),
        ]]),
    ]))
    mixed_computer = GO.RelateTopologyComputer(GO.relate_matrix_predicate(), GO.RelateGeometry(
        GI.Polygon([[
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.0, 0.0),
        ]]),
    ), mixed_b)
    @test GO.relate_is_self_noding_required(mixed_computer)
end

@testset "RelateNG polygon node section conversion" begin
    node_point = (0.0, 0.0)
    @test GO._relate_compare_angle(node_point, (1.0, 0.0), (0.0, 1.0)) < 0
    @test GO._relate_compare_angle(node_point, (2.0, 0.0), (1.0, 0.0)) == 0
    @test GO.relate_polygon_node_is_crossing(
        node_point,
        (-1.0, -1.0),
        (1.0, 1.0),
        (-1.0, 1.0),
        (1.0, -1.0),
    )
    @test !GO.relate_polygon_node_is_crossing(
        node_point,
        (-1.0, 0.0),
        (1.0, 0.0),
        (0.0, 1.0),
        (2.0, 0.0),
    )

    node = GO.RelateNode(node_point)
    GO.relate_add_line_edge!(node, GO.input_a, (0.0, 1.0))
    GO.relate_add_line_edge!(node, GO.input_a, (0.0, -1.0))
    GO.relate_add_line_edge!(node, GO.input_a, (1.0, 0.0))
    GO.relate_add_line_edge!(node, GO.input_a, (-1.0, 0.0))
    @test getproperty.(node.edges, :direction_point) ==
        [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)]

    shell = GO.RelateNodeSection(
        GO.input_a,
        GO.dim_area,
        7,
        0,
        nothing,
        true,
        node_point,
        (-1.0, 0.0),
        (0.0, 1.0),
    )
    hole = GO.RelateNodeSection(
        GO.input_a,
        GO.dim_area,
        7,
        1,
        nothing,
        true,
        node_point,
        (1.0, 0.0),
        (0.0, -1.0),
    )

    converted = GO.relate_convert_polygon_sections([shell, hole])
    @test length(converted) == 2
    @test all(section -> section.ring_id == 0, converted)
    @test getproperty.(converted, :previous_vertex) == [(-1.0, 0.0), (1.0, 0.0)]
    @test getproperty.(converted, :next_vertex) == [(0.0, -1.0), (0.0, 1.0)]

    hole_a = GO.RelateNodeSection(
        GO.input_a,
        GO.dim_area,
        8,
        1,
        nothing,
        true,
        node_point,
        (1.0, 0.0),
        (0.0, 1.0),
    )
    hole_b = GO.RelateNodeSection(
        GO.input_a,
        GO.dim_area,
        8,
        2,
        nothing,
        true,
        node_point,
        (-1.0, 0.0),
        (0.0, -1.0),
    )

    converted_holes = GO.relate_convert_polygon_sections([hole_a, hole_b])
    @test length(converted_holes) == 2
    @test all(section -> section.ring_id == 0, converted_holes)
    @test getproperty.(converted_holes, :previous_vertex) == [(1.0, 0.0), (-1.0, 0.0)]
    @test getproperty.(converted_holes, :next_vertex) == [(0.0, -1.0), (0.0, 1.0)]
end

@testset "RelateNG point path evaluation" begin
    alg = GO.RelateNG()
    point = GI.Point(1.0, 1.0)
    same_point = GI.Point(1.0, 1.0)
    other_point = GI.Point(2.0, 2.0)
    line = GI.LineString([(0.0, 1.0), (2.0, 1.0)])
    polygon = GI.Polygon([[
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.0, 0.0),
    ]])

    @test GO.intersects(alg, point, same_point)
    @test !GO.disjoint(alg, point, same_point)
    @test GO.equals(alg, point, same_point)
    @test GO.contains(alg, point, same_point)
    @test GO.de9im_string(GO.relate_matrix(alg, point, same_point)) == "0FFFFFFF2"

    @test GO.disjoint(alg, point, other_point)
    @test !GO.intersects(alg, point, other_point)
    @test !GO.equals(alg, point, other_point)
    @test GO.de9im_string(GO.relate_matrix(alg, point, other_point)) == "FF0FFF0F2"

    @test GO.within(alg, point, line)
    @test GO.coveredby(alg, point, line)
    @test GO.contains(alg, line, point)
    @test GO.covers(alg, line, point)
    @test GO.de9im_string(GO.relate_matrix(alg, point, line)) == "0FFFFF102"

    line_endpoint = GI.Point(0.0, 1.0)
    @test !GO.within(alg, line_endpoint, line)
    @test GO.coveredby(alg, line_endpoint, line)
    @test GO.touches(alg, line_endpoint, line)

    @test GO.within(alg, point, polygon)
    @test GO.coveredby(alg, point, polygon)
    @test GO.contains(alg, polygon, point)
    @test GO.covers(alg, polygon, point)
    @test GO.relate(alg, point, polygon, "T*F**F***")

    boundary_point = GI.Point(0.0, 1.0)
    outside_point = GI.Point(3.0, 3.0)
    @test !GO.within(alg, boundary_point, polygon)
    @test GO.coveredby(alg, boundary_point, polygon)
    @test GO.touches(alg, boundary_point, polygon)
    @test GO.disjoint(alg, outside_point, polygon)
    @test !GO.intersects(alg, outside_point, polygon)

    multipoint_a = GI.MultiPoint([(0.0, 0.0), (1.0, 1.0)])
    multipoint_b = GI.MultiPoint([(1.0, 1.0), (2.0, 2.0)])
    @test GO.overlaps(alg, multipoint_a, multipoint_b)
    @test !GO.equals(alg, multipoint_a, multipoint_b)

    crossing_points = GI.MultiPoint([(1.0, 1.0), (3.0, 3.0)])
    @test GO.crosses(alg, crossing_points, line)

    collection = GI.GeometryCollection([
        GI.Point(1.0, 1.0),
        GI.LineString([(0.0, 1.0), (2.0, 1.0)]),
        polygon,
    ])
    @test GO.within(alg, point, collection)
    @test GO.contains(alg, collection, point)
end

@testset "RelateNG mixed area-line collection evaluation" begin
    crossing_line_a = GI.LineString([(0.0, 0.0), (2.0, 2.0)])
    mixed_collection = GI.GeometryCollection([
        crossing_line_a,
        GI.Polygon([[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (0.0, 2.0),
            (0.0, 0.0),
        ]]),
    ])
    polygon = GI.Polygon([[
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.0, 0.0),
    ]])

    for alg in (GO.RelateNG(), GO.RelateNG(; prepared = true))
        @test GO.de9im_string(GO.relate_matrix(alg, mixed_collection, polygon)) ==
            "2FFF1FFF2"
        @test GO.contains(alg, polygon, mixed_collection)
        @test GO.contains(alg, mixed_collection, polygon)
        @test GO.equals(alg, mixed_collection, polygon)
    end
end

@testset "RelateNG mutual line edge events" begin
    alg = GO.RelateNG()
    crossing_line_a = GI.LineString([(0.0, 0.0), (2.0, 2.0)])
    crossing_line_b = GI.LineString([(0.0, 2.0), (2.0, 0.0)])

    @test GO.intersects(alg, crossing_line_a, crossing_line_b)
    @test GO.crosses(alg, crossing_line_a, crossing_line_b)
    @test !GO.touches(alg, crossing_line_a, crossing_line_b)
    @test GO.de9im_string(GO.relate_matrix(alg, crossing_line_a, crossing_line_b)) ==
        "0F1FF0102"

    prepared_alg = GO.RelateNG(; prepared = true)
    @test GO.de9im_string(GO.relate_matrix(prepared_alg, crossing_line_a, crossing_line_b)) ==
        GO.de9im_string(GO.relate_matrix(alg, crossing_line_a, crossing_line_b))
    prepared_a = GO.RelateGeometry(crossing_line_a; prepared = true)
    prepared_index = GO.relate_prepared_edge_index(prepared_a)
    @test prepared_index === GO.relate_prepared_edge_index(prepared_a)
    @test length(prepared_index.records) == 1

    touching_line = GI.LineString([(2.0, 2.0), (3.0, 2.0)])
    @test GO.intersects(alg, crossing_line_a, touching_line)
    @test GO.touches(alg, crossing_line_a, touching_line)
    @test !GO.crosses(alg, crossing_line_a, touching_line)
    @test GO.de9im_string(GO.relate_matrix(alg, crossing_line_a, touching_line)) ==
        "FF1F00102"

    overlapping_line = GI.LineString([(1.0, 1.0), (3.0, 3.0)])
    @test GO.intersects(alg, crossing_line_a, overlapping_line)
    @test GO.overlaps(alg, crossing_line_a, overlapping_line)
    @test !GO.equals(alg, crossing_line_a, overlapping_line)
    @test GO.de9im_string(GO.relate_matrix(alg, crossing_line_a, overlapping_line)) ==
        "1010F0102"

    equal_line = GI.LineString([(0.0, 0.0), (2.0, 2.0)])
    @test GO.equals(alg, crossing_line_a, equal_line)
    @test GO.contains(alg, crossing_line_a, equal_line)
    @test GO.covers(alg, crossing_line_a, equal_line)
    @test GO.de9im_string(GO.relate_matrix(alg, crossing_line_a, equal_line)) ==
        "1FFF0FFF2"

    self_crossing_line = GI.LineString([
        (0.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (2.0, 0.0),
    ])
    crossing_query = GI.LineString([(0.0, 1.0), (2.0, 1.0)])
    @test GO.crosses(alg, self_crossing_line, crossing_query)
    @test GO.de9im_string(GO.relate_matrix(alg, self_crossing_line, crossing_query)) ==
        "0F1FF0102"
    @test GO.de9im_string(GO.relate_matrix(prepared_alg, self_crossing_line, crossing_query)) ==
        "0F1FF0102"
end

@testset "RelateNG local node topology for line and area edges" begin
    alg = GO.RelateNG()
    polygon = GI.Polygon([[
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.0, 0.0),
    ]])
    crossing_line = GI.LineString([(-1.0, 1.0), (3.0, 1.0)])
    interior_line = GI.LineString([(0.5, 1.0), (1.5, 1.0)])
    boundary_line = GI.LineString([(2.0, 0.0), (2.0, 1.0)])

    @test GO.de9im_string(GO.relate_matrix(alg, crossing_line, polygon)) == "101FF0212"
    @test GO.crosses(alg, crossing_line, polygon)
    @test !GO.within(alg, crossing_line, polygon)

    @test GO.de9im_string(GO.relate_matrix(alg, interior_line, polygon)) == "1FF0FF212"
    @test GO.within(alg, interior_line, polygon)
    @test GO.contains(alg, polygon, interior_line)

    @test GO.de9im_string(GO.relate_matrix(alg, boundary_line, polygon)) == "F1FF0F212"
    @test GO.touches(alg, boundary_line, polygon)
    @test GO.coveredby(alg, boundary_line, polygon)

    overlapping_polygon = GI.Polygon([[
        (1.0, 1.0),
        (3.0, 1.0),
        (3.0, 3.0),
        (1.0, 3.0),
        (1.0, 1.0),
    ]])
    adjacent_polygon = GI.Polygon([[
        (2.0, 0.0),
        (3.0, 0.0),
        (3.0, 1.0),
        (2.0, 1.0),
        (2.0, 0.0),
    ]])

    @test GO.de9im_string(GO.relate_matrix(alg, polygon, overlapping_polygon)) == "212101212"
    @test GO.overlaps(alg, polygon, overlapping_polygon)
    @test !GO.touches(alg, polygon, overlapping_polygon)

    @test GO.de9im_string(GO.relate_matrix(alg, polygon, adjacent_polygon)) == "FF2F11212"
    @test GO.touches(alg, polygon, adjacent_polygon)
    @test !GO.overlaps(alg, polygon, adjacent_polygon)
end

@testset "RelateNG interaction predicates with area edges" begin
    alg = GO.RelateNG()
    polygon = GI.Polygon([[
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.0, 0.0),
    ]])
    crossing_line = GI.LineString([(-1.0, 1.0), (3.0, 1.0)])
    interior_line = GI.LineString([(0.5, 1.0), (1.5, 1.0)])
    exterior_line = GI.LineString([(3.0, 3.0), (4.0, 4.0)])

    @test GO.intersects(alg, crossing_line, polygon)
    @test !GO.disjoint(alg, crossing_line, polygon)
    @test GO.intersects(alg, interior_line, polygon)
    @test GO.disjoint(alg, exterior_line, polygon)

    overlapping_polygon = GI.Polygon([[
        (1.0, 1.0),
        (3.0, 1.0),
        (3.0, 3.0),
        (1.0, 3.0),
        (1.0, 1.0),
    ]])
    exterior_polygon = GI.Polygon([[
        (3.0, 3.0),
        (4.0, 3.0),
        (4.0, 4.0),
        (3.0, 4.0),
        (3.0, 3.0),
    ]])

    @test GO.intersects(alg, polygon, overlapping_polygon)
    @test !GO.disjoint(alg, polygon, overlapping_polygon)
    @test !GO.intersects(alg, polygon, exterior_polygon)
    @test GO.disjoint(alg, polygon, exterior_polygon)
end
