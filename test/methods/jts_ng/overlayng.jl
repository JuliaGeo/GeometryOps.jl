using Test
import GeometryOps as GO
import GeoInterface as GI
using GeometryOpsTestHelpers

_overlay_tuples(geoms) = map(geom -> GO.tuples(geom), geoms)
_ring_tuples(poly, i = 1) = [GO.tuples(point) for point in GI.getpoint(GI.getring(poly, i))]
_overlay_area(geoms) = sum(geom -> GO.area(geom), geoms; init = 0.0)
_result_points(geoms) = [GO.tuples(geom) for geom in geoms if GI.trait(geom) isa GI.PointTrait]
_result_lines(geoms) = [
    [GO.tuples(point) for point in GI.getpoint(geom)] for geom in geoms if GI.trait(geom) isa GI.LineStringTrait
]

const _JTS_OVERLAY_FIXTURE_DIR = normpath(joinpath(
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

function _overlayng_fixture_value(alg::GO.OverlayNG, op::JTSOperation)
    name = lowercase(op.name)
    a = _overlayng_fixture_argument(op.arguments[1])
    b = _overlayng_fixture_argument(op.arguments[2])
    if name == "intersection" || name == "intersectionng"
        return GO.intersection(alg, a, b)
    elseif name == "union" || name == "unionng"
        return GO.union(alg, a, b)
    elseif name == "difference" || name == "differenceng"
        return GO.difference(alg, a, b)
    elseif name == "symdifference" || name == "symdifferenceng"
        return GO.symdifference(alg, a, b)
    end
    error("Unsupported OverlayNG fixture operation: $(op.name)")
end

_overlayng_fixture_argument(geom) = geom
_overlayng_fixture_argument(::JTSEmptyGeometry) = GI.FeatureCollection(Any[])

function _overlayng_flatten_components(value)
    value isa AbstractVector && return reduce(vcat, map(_overlayng_flatten_components, value); init = Any[])
    value isa JTSEmptyGeometry && return Any[]
    isnothing(value) && return Any[]

    trait = GI.trait(value)
    if trait isa GI.GeometryCollectionTrait ||
       trait isa GI.MultiPointTrait ||
       trait isa GI.MultiLineStringTrait ||
       trait isa GI.MultiPolygonTrait
        return reduce(vcat, map(_overlayng_flatten_components, GI.getgeom(value)); init = Any[])
    end
    return Any[value]
end

_overlayng_segment_key(a, b) = min((a, b), (b, a))

function _overlayng_line_segments(line)
    points = [GO.tuples(point) for point in GI.getpoint(line)]
    segments = Tuple[]
    for index in 1:(length(points) - 1)
        points[index] == points[index + 1] && continue
        push!(segments, _overlayng_segment_key(points[index], points[index + 1]))
    end
    return segments
end

function _overlayng_fixture_summary(value)
    points = Set{Any}()
    line_segments = Set{Any}()
    areas = Float64[]
    ring_counts = Int[]
    for geom in _overlayng_flatten_components(value)
        trait = GI.trait(geom)
        if trait isa GI.PointTrait
            push!(points, GO.tuples(geom))
        elseif trait isa GI.LineTrait ||
               trait isa GI.LineStringTrait ||
               trait isa GI.LinearRingTrait
            union!(line_segments, _overlayng_line_segments(geom))
        elseif trait isa GI.PolygonTrait
            push!(areas, round(GO.area(geom); digits = 8))
            push!(ring_counts, GI.nring(geom))
        else
            error("Unsupported OverlayNG fixture result component: $(typeof(geom))")
        end
    end
    return (; points, line_segments, areas = sort(areas), ring_counts = sort(ring_counts))
end

_overlayng_fixture_operation_allowed(::Nothing, op::JTSOperation) = true
_overlayng_fixture_operation_allowed(operations, op::JTSOperation) =
    lowercase(op.name) in operations

_overlayng_fixture_algorithm(test_set::JTSTestSet) =
    GO.OverlayNG(; precision_model = _overlayng_fixture_precision_model(test_set.precision_model))

_overlayng_fixture_precision_model(::Nothing) = GO.NoPrecisionModel()

function _overlayng_fixture_precision_model(model::JTSPrecisionModel)
    if uppercase(model.model_type) == "FLOATING" && isnothing(model.scale)
        return GO.NoPrecisionModel()
    end
    scale = isnothing(model.scale) ? 1.0 : model.scale
    offset = (something(model.offsetx, 0.0), something(model.offsety, 0.0))
    return GO.FixedPrecisionModel(scale; offset)
end

@testset "OverlayNG input wrappers" begin
    alg = GO.OverlayNG()
    point_input = GO.OverlayInputGeometry(alg, GI.MultiPoint([(0.0, 0.0), (1.0, 1.0)]))
    line_input = GO.OverlayInputGeometry(alg, GI.LineString([(0.0, 0.0), (1.0, 1.0)]))

    @test point_input.dimension == GO.dim_point
    @test isnothing(point_input.locator)
    @test line_input.dimension == GO.dim_line
    @test line_input.locator isa GO.RelatePointLocator
    @test GO.overlay_has_point_dispatch(point_input, line_input)
end

@testset "OverlayNG edge source extraction" begin
    alg = GO.OverlayNG()
    line = GI.LineString([(0.0, 0.0), (1.0, 1.0)])
    line_input = GO.OverlayInputGeometry(alg, line)
    line_segments = GO.overlay_segment_strings(line_input; input_side = GO.input_b)

    @test line_segments === GO.overlay_segment_strings(line_input; input_side = GO.input_b)
    @test length(line_segments) == 1
    line_segment = only(line_segments)
    @test line_segment isa GO.OverlaySegmentString
    @test line_segment.points == [(0.0, 0.0), (1.0, 1.0)]
    @test line_segment.source isa GO.OverlayEdgeSourceInfo
    @test line_segment.source.input_side == GO.input_b
    @test line_segment.source.source_dimension == GO.dim_line
    @test line_segment.source.ring_role == GO.ring_none
    @test line_segment.source.depth_delta == 0
    @test !line_segment.source.is_collapsed
    @test line_segment.source.geometry === line

    polygon = GI.Polygon([
        [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)],
        [(1.0, 1.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0), (1.0, 1.0)],
    ])
    polygon_input = GO.OverlayInputGeometry(alg, polygon)
    shell, hole = GO.overlay_segment_strings(polygon_input; input_side = GO.input_a)

    @test shell.source.source_dimension == GO.dim_area
    @test shell.source.ring_role == GO.ring_shell
    @test hole.source.ring_role == GO.ring_hole
    @test shell.source.source_orientation == GO.ring_counterclockwise
    @test hole.source.source_orientation == GO.ring_clockwise
    @test !shell.source.coordinates_reversed
    @test !hole.source.coordinates_reversed
    @test shell.source.depth_delta == -1
    @test hole.source.depth_delta == -1
    @test shell.source.parent_polygonal === polygon
    @test hole.source.parent_polygonal === polygon
end

@testset "OverlayNG simple noder" begin
    alg = GO.OverlayNG()
    crossing_a = GI.LineString([(0.0, 0.0), (2.0, 2.0)])
    crossing_b = GI.LineString([(0.0, 2.0), (2.0, 0.0)])

    raw_input_a = GO.OverlayInputGeometry(alg, crossing_a)
    raw_input_b = GO.OverlayInputGeometry(alg, crossing_b)
    raw_segments = Any[
        GO.overlay_segment_strings(raw_input_a; input_side = GO.input_a)...,
        GO.overlay_segment_strings(raw_input_b; input_side = GO.input_b)...,
    ]
    @test !GO.overlay_is_fully_noded(raw_segments)

    noded = GO.overlay_node_segment_strings(alg, raw_input_a, raw_input_b)
    @test length(noded) == 4
    @test GO.overlay_is_fully_noded(noded)
    @test getproperty.(noded, :points) == [
        [(0.0, 0.0), (1.0, 1.0)],
        [(1.0, 1.0), (2.0, 2.0)],
        [(0.0, 2.0), (1.0, 1.0)],
        [(1.0, 1.0), (2.0, 0.0)],
    ]
    @test count(segment -> segment.source.input_side == GO.input_a, noded) == 2
    @test count(segment -> segment.source.input_side == GO.input_b, noded) == 2

    overlap_a = GI.LineString([(0.0, 0.0), (3.0, 0.0)])
    overlap_b = GI.LineString([(1.0, 0.0), (2.0, 0.0)])
    overlap_noded = GO.overlay_node_segment_strings(alg, overlap_a, overlap_b)
    @test GO.overlay_is_fully_noded(overlap_noded)
    @test getproperty.(overlap_noded, :points) == [
        [(0.0, 0.0), (1.0, 0.0)],
        [(1.0, 0.0), (2.0, 0.0)],
        [(2.0, 0.0), (3.0, 0.0)],
        [(1.0, 0.0), (2.0, 0.0)],
    ]

    repeated = GI.LineString([(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (2.0, 0.0)])
    crossing = GI.LineString([(1.0, -1.0), (1.0, 1.0)])
    repeated_noded = GO.overlay_node_segment_strings(alg, repeated, crossing)
    @test all(segment -> !segment.is_zero_length, repeated_noded)
    @test any(segment -> segment.had_repeated_coordinates, repeated_noded)
end

@testset "OverlayNG edge merger" begin
    alg = GO.OverlayNG()
    line_a = GI.LineString([(0.0, 0.0), (2.0, 0.0)])
    line_b = GI.LineString([(2.0, 0.0), (0.0, 0.0)])

    merged_line_edges = GO.overlay_merge_edges(alg, line_a, line_b)
    @test length(merged_line_edges) == 1
    line_edge = only(merged_line_edges)
    @test line_edge.key == GO.OverlayEdgeKey((0.0, 0.0), (2.0, 0.0))
    @test line_edge.points == [(0.0, 0.0), (2.0, 0.0)]
    @test length(line_edge.sources) == 2
    @test line_edge.source_directions == [true, false]
    @test line_edge.depth_delta == 0
    @test GO.overlay_is_line_edge(line_edge)
    @test !GO.overlay_is_boundary_edge(line_edge)
    @test GO.overlay_primary_ring_role(line_edge) == GO.ring_none

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
    merged_area_edges = GO.overlay_merge_edges(alg, left, right)
    shared_key = GO.OverlayEdgeKey((1.0, 0.0), (1.0, 1.0))
    shared_edge = only(filter(edge -> edge.key == shared_key, merged_area_edges))

    @test length(shared_edge.sources) == 2
    @test shared_edge.source_directions == [true, false]
    @test shared_edge.depth_delta == 0
    @test GO.overlay_is_boundary_edge(shared_edge)
    @test !GO.overlay_is_line_edge(shared_edge)
    @test GO.overlay_primary_ring_role(shared_edge) == GO.ring_shell
end

@testset "OverlayNG edge labels" begin
    alg = GO.OverlayNG()
    line_a = GI.LineString([(0.0, 0.0), (2.0, 0.0)])
    line_b = GI.LineString([(2.0, 0.0), (0.0, 0.0)])
    line_edge = only(GO.overlay_merge_edges(alg, line_a, line_b))
    line_label = GO.overlay_label(line_edge)

    @test line_label.input_a.dimension == GO.dim_line
    @test line_label.input_a.on_location == GO.loc_interior
    @test line_label.input_a.left_location == GO.loc_exterior
    @test line_label.input_a.right_location == GO.loc_exterior
    @test line_label.input_a.line_state == GO.overlay_line_part
    @test line_label.input_a.collapse_role == GO.overlay_not_collapsed
    @test line_label.input_b.line_state == GO.overlay_line_part

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
    shared_key = GO.OverlayEdgeKey((1.0, 0.0), (1.0, 1.0))
    shared_edge = only(filter(edge -> edge.key == shared_key, GO.overlay_merge_edges(alg, left, right)))
    label = GO.overlay_label(shared_edge)

    @test label.input_a.dimension == GO.dim_area
    @test label.input_a.on_location == GO.loc_interior
    @test label.input_a.left_location == GO.loc_interior
    @test label.input_a.right_location == GO.loc_exterior
    @test label.input_a.line_state == GO.overlay_boundary_part
    @test label.input_a.collapse_role == GO.overlay_not_collapsed

    @test label.input_b.dimension == GO.dim_area
    @test label.input_b.on_location == GO.loc_interior
    @test label.input_b.left_location == GO.loc_exterior
    @test label.input_b.right_location == GO.loc_interior
    @test label.input_b.line_state == GO.overlay_boundary_part

    single_line_edge = only(GO.overlay_merge_edges(
        GO.overlay_node_segment_strings(alg, line_a, GI.MultiPoint([(10.0, 10.0)])),
    ))
    missing_b_label = GO.overlay_input_label(single_line_edge, GO.input_b)
    @test missing_b_label.dimension == GO.dim_false
    @test missing_b_label.line_state == GO.overlay_not_part
    @test missing_b_label.on_location == GO.loc_exterior
end

@testset "OverlayNG half-edge graph" begin
    alg = GO.OverlayNG()
    line_a = GI.LineString([(0.0, 0.0), (2.0, 2.0)])
    line_b = GI.LineString([(0.0, 2.0), (2.0, 0.0)])
    edges = GO.overlay_merge_edges(alg, line_a, line_b)
    graph = GO.overlay_graph(edges)

    @test length(edges) == 4
    @test length(graph.half_edges) == 8
    @test all(half_edge -> half_edge.sym.sym === half_edge, graph.half_edges)
    @test all(half_edge -> half_edge.label isa GO.OverlayLabel, graph.half_edges)
    @test all(half_edge -> !half_edge.result_area && !half_edge.result_line, graph.half_edges)
    @test all(half_edge -> !half_edge.visited && half_edge.ring_id == 0, graph.half_edges)

    center_star = GO.overlay_node_star(graph, (1.0, 1.0))
    @test length(center_star) == 4
    @test getproperty.(center_star, :destination) == [
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
    ]
    @test issorted(getproperty.(center_star, :angle))

    corner_star = GO.overlay_node_star(graph, (0.0, 0.0))
    @test length(corner_star) == 1
    @test only(corner_star).destination == (1.0, 1.0)
    @test isempty(GO.overlay_node_star(graph, (10.0, 10.0)))
end

@testset "OverlayNG result edge marking" begin
    alg = GO.OverlayNG()
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
    graph = GO.overlay_graph(GO.overlay_merge_edges(alg, left, right))
    shared_key = GO.OverlayEdgeKey((1.0, 0.0), (1.0, 1.0))
    shared_half_edges = filter(half_edge -> half_edge.edge.key == shared_key, graph.half_edges)

    duplicate_graph = GO.overlay_graph(GO.overlay_merge_edges(alg, left, right))
    duplicate_shared_half_edges = filter(
        half_edge -> half_edge.edge.key == shared_key,
        duplicate_graph.half_edges,
    )
    GO.overlay_mark_result_area_edges!(duplicate_graph, GO.overlay_union)
    @test all(half_edge -> half_edge.result_area, duplicate_shared_half_edges)
    GO.overlay_unmark_duplicate_result_area_edges!(duplicate_graph)
    @test all(half_edge -> !half_edge.result_area, duplicate_shared_half_edges)

    GO.overlay_mark_result_edges!(graph, GO.overlay_union)
    @test count(half_edge -> half_edge.result_area, graph.half_edges) == 6
    @test all(half_edge -> !half_edge.result_area, shared_half_edges)

    GO.overlay_mark_result_edges!(graph, GO.overlay_intersection)
    @test count(half_edge -> half_edge.result_area, graph.half_edges) == 0

    GO.overlay_mark_result_edges!(graph, GO.overlay_difference)
    @test count(half_edge -> half_edge.result_area, graph.half_edges) == 4
    @test count(half_edge -> half_edge.result_area, shared_half_edges) == 1

    line_a = GI.LineString([(0.0, 0.0), (2.0, 0.0)])
    line_b = GI.LineString([(2.0, 0.0), (0.0, 0.0)])
    line_graph = GO.overlay_graph(GO.overlay_merge_edges(alg, line_a, line_b))

    GO.overlay_mark_result_edges!(line_graph, GO.overlay_intersection)
    @test count(half_edge -> half_edge.result_line, line_graph.half_edges) == 2

    GO.overlay_mark_result_edges!(line_graph, GO.overlay_difference)
    @test count(half_edge -> half_edge.result_line, line_graph.half_edges) == 0

    single_line_graph = GO.overlay_graph(GO.overlay_merge_edges(
        GO.overlay_node_segment_strings(alg, line_a, GI.MultiPoint([(10.0, 10.0)])),
    ))

    GO.overlay_mark_result_edges!(single_line_graph, GO.overlay_difference)
    @test count(half_edge -> half_edge.result_line, single_line_graph.half_edges) == 2

    GO.overlay_mark_result_edges!(single_line_graph, GO.overlay_intersection)
    @test count(half_edge -> half_edge.result_line, single_line_graph.half_edges) == 0
end

@testset "OverlayNG point-point dispatch" begin
    alg = GO.OverlayNG()
    points_a = GI.MultiPoint([(0.0, 0.0), (1.0, 1.0), (1.0, 1.0)])
    points_b = GI.MultiPoint([(1.0, 1.0), (2.0, 2.0)])

    @test _overlay_tuples(GO.intersection(alg, points_a, points_b)) == [(1.0, 1.0)]
    @test _overlay_tuples(GO.union(alg, points_a, points_b)) ==
        [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)]
    @test _overlay_tuples(GO.difference(alg, points_a, points_b)) == [(0.0, 0.0)]
    @test _overlay_tuples(GO.symdifference(alg, points_a, points_b)) ==
        [(0.0, 0.0), (2.0, 2.0)]
end

@testset "OverlayNG point-line dispatch" begin
    alg = GO.OverlayNG()
    points = GI.MultiPoint([(0.0, 0.0), (1.0, 0.0), (3.0, 0.0)])
    line = GI.LineString([(0.0, 0.0), (2.0, 0.0)])

    @test _overlay_tuples(GO.intersection(alg, points, line)) ==
        [(0.0, 0.0), (1.0, 0.0)]
    @test _overlay_tuples(GO.difference(alg, points, line)) == [(3.0, 0.0)]

    line_minus_points = GO.difference(alg, line, points)
    @test length(line_minus_points) == 1
    @test line_minus_points[1] === line

    point_union = GO.union(alg, points, line)
    @test length(point_union) == 2
    @test point_union[1] === line
    @test GO.tuples(point_union[2]) == (3.0, 0.0)

    point_target_union = GO.union(alg, points, line; target = GI.PointTrait())
    @test _overlay_tuples(point_target_union) == [(3.0, 0.0)]

    point_symdiff = GO.symdifference(alg, points, line)
    @test length(point_symdiff) == 2
    @test point_symdiff[1] === line
    @test GO.tuples(point_symdiff[2]) == (3.0, 0.0)

    self_crossing_line = GI.LineString([
        (20.0, 20.0),
        (110.0, 110.0),
        (170.0, 50.0),
        (130.0, 10.0),
        (70.0, 70.0),
    ])
    self_crossing_points = GI.MultiPoint([(40.0, 90.0), (20.0, 20.0), (70.0, 70.0)])
    self_crossing_union = GO.union(alg, self_crossing_points, self_crossing_line)
    @test _result_points(self_crossing_union) == [(40.0, 90.0)]
    @test Set(_overlayng_segment_key(first(line), last(line)) for line in _result_lines(self_crossing_union)) == Set([
        _overlayng_segment_key((20.0, 20.0), (70.0, 70.0)),
        _overlayng_segment_key((70.0, 70.0), (110.0, 110.0)),
        _overlayng_segment_key((110.0, 110.0), (170.0, 50.0)),
        _overlayng_segment_key((170.0, 50.0), (130.0, 10.0)),
        _overlayng_segment_key((130.0, 10.0), (70.0, 70.0)),
    ])
end

@testset "OverlayNG point-area dispatch" begin
    alg = GO.OverlayNG()
    polygon = GI.Polygon([[
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.0, 0.0),
    ]])
    points = GI.MultiPoint([(1.0, 1.0), (0.0, 1.0), (3.0, 3.0)])

    @test _overlay_tuples(GO.intersection(alg, points, polygon)) ==
        [(1.0, 1.0), (0.0, 1.0)]
    @test _overlay_tuples(GO.difference(alg, points, polygon)) == [(3.0, 3.0)]

    union_result = GO.union(alg, polygon, points)
    @test length(union_result) == 2
    @test union_result[1] === polygon
    @test GO.tuples(union_result[2]) == (3.0, 3.0)

    area_only_union = GO.union(GO.OverlayNG(; area_result_only = true), polygon, points)
    @test area_only_union == Any[polygon]
end

@testset "OverlayNG empty input dispatch" begin
    alg = GO.OverlayNG()
    empty = GI.FeatureCollection(Any[])
    polygon = GI.Polygon([[
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.0, 0.0),
    ]])
    line = GI.LineString([(0.0, 0.0), (2.0, 0.0)])
    point = GI.Point(1.0, 1.0)

    @test isempty(GO.intersection(alg, polygon, empty))
    @test GO.union(alg, polygon, empty) == Any[polygon]
    @test GO.difference(alg, polygon, empty) == Any[polygon]
    @test isempty(GO.difference(alg, empty, polygon))
    @test GO.symdifference(alg, empty, polygon) == Any[polygon]

    @test GO.union(alg, empty, line) == Any[line]
    @test GO.difference(alg, line, empty) == Any[line]
    @test isempty(GO.intersection(alg, empty, point))
    @test GO.symdifference(alg, empty, point) == Any[point]
    @test isempty(GO.union(alg, empty, empty))
end

@testset "OverlayNG area-area dispatch" begin
    alg = GO.OverlayNG()
    square_a = GI.Polygon([[
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.0, 0.0),
    ]])
    square_b = GI.Polygon([[
        (1.0, 1.0),
        (3.0, 1.0),
        (3.0, 3.0),
        (1.0, 3.0),
        (1.0, 1.0),
    ]])

    intersection_result = GO.intersection(alg, square_a, square_b)
    @test length(intersection_result) == 1
    @test GI.nring(only(intersection_result)) == 1
    @test _overlay_area(intersection_result) ≈ 1.0
    @test Set(_ring_tuples(only(intersection_result))) == Set([
        (1.0, 1.0),
        (2.0, 1.0),
        (2.0, 2.0),
        (1.0, 2.0),
    ])

    adjacent_a = GI.Polygon([[
        (0.0, 0.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (0.0, 0.0),
    ]])
    adjacent_b = GI.Polygon([[
        (1.0, 0.0),
        (2.0, 0.0),
        (2.0, 1.0),
        (1.0, 1.0),
        (1.0, 0.0),
    ]])
    union_result = GO.union(alg, adjacent_a, adjacent_b)
    @test length(union_result) == 1
    @test GI.nring(only(union_result)) == 1
    @test _overlay_area(union_result) ≈ 2.0

    outer = GI.Polygon([[
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 4.0),
        (0.0, 4.0),
        (0.0, 0.0),
    ]])
    inner = GI.Polygon([[
        (1.0, 1.0),
        (3.0, 1.0),
        (3.0, 3.0),
        (1.0, 3.0),
        (1.0, 1.0),
    ]])
    hole_result = GO.difference(alg, outer, inner)
    @test length(hole_result) == 1
    @test GI.nring(only(hole_result)) == 2
    @test _overlay_area(hole_result) ≈ 12.0

    touch_line_point_a = GI.Polygon([[
        (10.0, 10.0),
        (10.0, 30.0),
        (30.0, 30.0),
        (30.0, 10.0),
        (10.0, 10.0),
    ]])
    touch_line_point_b = GI.Polygon([[
        (40.0, 25.0),
        (30.0, 25.0),
        (30.0, 20.0),
        (35.0, 20.0),
        (30.0, 15.0),
        (40.0, 15.0),
        (40.0, 25.0),
    ]])
    touch_intersection = GO.intersection(alg, touch_line_point_a, touch_line_point_b)
    @test Set(_result_lines(touch_intersection)) == Set([[(30.0, 20.0), (30.0, 25.0)]])
    @test _result_points(touch_intersection) == [(30.0, 15.0)]
    @test _overlay_area(touch_intersection) == 0.0

    touch_union = GO.union(alg, touch_line_point_a, touch_line_point_b)
    @test length(touch_union) == 1
    @test GI.nring(only(touch_union)) == 2
    @test _overlay_area(touch_union) ≈ 487.5

    touch_symdifference = GO.symdifference(alg, touch_line_point_a, touch_line_point_b)
    @test length(touch_symdifference) == 1
    @test GI.nring(only(touch_symdifference)) == 2
    @test _overlay_area(touch_symdifference) ≈ 487.5

    overlap_touch_a = touch_line_point_a
    overlap_touch_b = GI.Polygon([[
        (40.0, 25.0),
        (25.0, 25.0),
        (35.0, 15.0),
        (30.0, 15.0),
        (30.0, 10.0),
        (40.0, 10.0),
        (40.0, 25.0),
    ]])
    overlap_touch_symdifference = GO.symdifference(alg, overlap_touch_a, overlap_touch_b)
    @test length(overlap_touch_symdifference) == 1
    @test GI.nring(only(overlap_touch_symdifference)) == 3
    @test _overlay_area(overlap_touch_symdifference) ≈ 525.0
end

@testset "OverlayNG line result extraction" begin
    alg = GO.OverlayNG()
    line_a = GI.LineString([(0.0, 0.0), (2.0, 0.0)])
    line_b = GI.LineString([(1.0, 0.0), (3.0, 0.0)])

    @test Set(_result_lines(GO.intersection(alg, line_a, line_b))) ==
        Set([[(1.0, 0.0), (2.0, 0.0)]])
    @test Set(_result_lines(GO.union(alg, line_a, line_b))) == Set([
        [(0.0, 0.0), (1.0, 0.0)],
        [(1.0, 0.0), (2.0, 0.0)],
        [(2.0, 0.0), (3.0, 0.0)],
    ])
    @test Set(_result_lines(GO.difference(alg, line_a, line_b))) ==
        Set([[(0.0, 0.0), (1.0, 0.0)]])
    @test Set(_result_lines(GO.symdifference(alg, line_a, line_b))) == Set([
        [(0.0, 0.0), (1.0, 0.0)],
        [(2.0, 0.0), (3.0, 0.0)],
    ])

    crossing_a = GI.LineString([(0.0, 0.0), (2.0, 2.0)])
    crossing_b = GI.LineString([(0.0, 2.0), (2.0, 0.0)])
    crossing_intersection = GO.intersection(alg, crossing_a, crossing_b)
    @test _result_points(crossing_intersection) == [(1.0, 1.0)]
    @test isempty(_result_lines(crossing_intersection))

    mixed_a = GI.LineString([(0.0, 0.0), (2.0, 0.0), (3.0, 1.0), (3.0, -1.0)])
    mixed_b = GI.LineString([(1.0, 0.0), (4.0, 0.0)])
    mixed_intersection = GO.intersection(alg, mixed_a, mixed_b)
    @test Set(_result_lines(mixed_intersection)) == Set([[(1.0, 0.0), (2.0, 0.0)]])
    @test _result_points(mixed_intersection) == [(3.0, 0.0)]

    strict_intersection = GO.intersection(GO.OverlayNG(; strict = true), mixed_a, mixed_b)
    @test Set(_result_lines(strict_intersection)) == Set([[(1.0, 0.0), (2.0, 0.0)]])
    @test isempty(_result_points(strict_intersection))
end

@testset "OverlayNG line-area extraction" begin
    alg = GO.OverlayNG()
    square = GI.Polygon([[
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.0, 0.0),
    ]])
    crossing = GI.LineString([(-1.0, 1.0), (3.0, 1.0)])

    @test Set(_result_lines(GO.intersection(alg, crossing, square))) ==
        Set([[(0.0, 1.0), (2.0, 1.0)]])
    @test Set(_result_lines(GO.difference(alg, crossing, square))) == Set([
        [(-1.0, 1.0), (0.0, 1.0)],
        [(2.0, 1.0), (3.0, 1.0)],
    ])

    touch = GI.LineString([(-1.0, -1.0), (0.0, 0.0)])
    touch_intersection = GO.intersection(alg, touch, square)
    @test _result_points(touch_intersection) == [(0.0, 0.0)]
    @test isempty(_result_lines(touch_intersection))

    fixed_alg = GO.OverlayNG(; precision_model = GO.FixedPrecisionModel(1.0))
    line = GI.LineString([(240.0, 190.0), (120.0, 120.0)])
    triangle = GI.Polygon([[
        (110.0, 240.0),
        (50.0, 80.0),
        (240.0, 70.0),
        (110.0, 240.0),
    ]])
    @test Set(_overlayng_segment_key(first(line), last(line)) for line in _result_lines(
        GO.intersection(fixed_alg, line, triangle),
    )) == Set([_overlayng_segment_key((177.0, 153.0), (120.0, 120.0))])
    @test Set(_overlayng_segment_key(first(line), last(line)) for line in _result_lines(
        GO.difference(fixed_alg, line, triangle),
    )) == Set([_overlayng_segment_key((240.0, 190.0), (177.0, 153.0))])

    sliver = GI.Polygon([[
        (95.0, 9.0),
        (81.0, 414.0),
        (87.0, 414.0),
        (95.0, 9.0),
    ]])
    sliver_line = GI.LineString([(93.0, 13.0), (96.0, 13.0)])
    @test Set(_overlayng_segment_key(first(line), last(line)) for line in _result_lines(
        GO.difference(fixed_alg, sliver, sliver_line),
    )) == Set([_overlayng_segment_key((95.0, 9.0), (95.0, 13.0))])
end

@testset "OverlayNG JTS XML smoke fixtures" begin
    alg = GO.OverlayNG()
    fixtures = (
        ("TestNGOverlayP.xml", 1:12, nothing),
        ("TestNGOverlayL.xml", 1:11, nothing),
        ("TestNGOverlayA.xml", 1:20, nothing),
        ("TestNGOverlayEmpty.xml", 1:16, nothing),
        ("TestNGOverlayGC.xml", 1:4, nothing),
        ("TestOverlayPP.xml", 1:8, nothing),
        ("TestOverlayPL.xml", 1:5, nothing),
        ("TestOverlayPA.xml", 1:3, nothing),
        ("TestOverlayLL.xml", 1:7, nothing),
        ("TestOverlayLA.xml", 1:4, nothing),
        ("TestOverlayAA.xml", 1:13, nothing),
        ("TestOverlayEmpty.xml", 1:144, nothing),
        ("TestOverlayLLPrec.xml", 1:2, nothing),
        ("TestOverlayLAPrec.xml", 1:4, nothing),
        ("TestOverlayAAPrec.xml", 1:1, nothing),
        ("TestOverlayAAPrec.xml", 2:2, ("intersection",)),
        ("TestOverlayAAPrec.xml", 3:5, ("union", "difference", "symdifference")),
        ("TestOverlayAAPrec.xml", 6:6, nothing),
        ("TestOverlayAAPrec.xml", 7:7, ("intersection", "union")),
        ("TestOverlayAAPrec.xml", 11:11, nothing),
        ("TestOverlayAAPrec.xml", 12:12, ("intersection",)),
        ("TestOverlayAAPrec.xml", 13:14, nothing),
        ("TestOverlayAAPrec.xml", 15:16, ("intersection",)),
        ("TestOverlayAAPrec.xml", 17:17, nothing),
        ("TestOverlayAAPrec.xml", 18:18, ("intersection", "symdifference")),
    )

    matched_operations = 0
    for (filename, case_indices, operations) in fixtures
        test_set = load_test_set(joinpath(_JTS_OVERLAY_FIXTURE_DIR, filename))
        alg = _overlayng_fixture_algorithm(test_set)
        for case_index in case_indices
            case = test_set.cases[case_index]
            for op in case.operations
                _overlayng_fixture_operation_allowed(operations, op) || continue
                matched_operations += 1
                @test _overlayng_fixture_summary(_overlayng_fixture_value(alg, op)) ==
                    _overlayng_fixture_summary(op.expected)
            end
        end
    end
    @test matched_operations == 798
end
