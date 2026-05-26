using Test
import GeometryOps as GO
import GeoInterface as GI

_overlay_tuples(geoms) = map(geom -> GO.tuples(geom), geoms)

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
    @test count(half_edge -> half_edge.result_line, line_graph.half_edges) == 1

    GO.overlay_mark_result_edges!(line_graph, GO.overlay_difference)
    @test count(half_edge -> half_edge.result_line, line_graph.half_edges) == 0

    single_line_graph = GO.overlay_graph(GO.overlay_merge_edges(
        GO.overlay_node_segment_strings(alg, line_a, GI.MultiPoint([(10.0, 10.0)])),
    ))

    GO.overlay_mark_result_edges!(single_line_graph, GO.overlay_difference)
    @test count(half_edge -> half_edge.result_line, single_line_graph.half_edges) == 1

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

@testset "OverlayNG unsupported edge overlay" begin
    alg = GO.OverlayNG()
    line_a = GI.LineString([(0.0, 0.0), (1.0, 1.0)])
    line_b = GI.LineString([(0.0, 1.0), (1.0, 0.0)])

    @test_throws ArgumentError GO.intersection(alg, line_a, line_b)
end
