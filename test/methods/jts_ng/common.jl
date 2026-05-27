using Test
import GeometryOps as GO
import GeoInterface as GI
import GeoInterface.Extents: Extents

struct TestGridPrecision end
function GO.apply_ng_precision(::TestGridPrecision, point, ::Type{T} = Float64) where {T}
    return (round(T, GI.x(point)), round(T, GI.y(point)))
end

@testset "Algorithm markers" begin
    relate = GO.RelateNG()
    overlay = GO.OverlayNG()
    fixed_overlay = GO.OverlayNG(; precision_model = GO.FixedPrecisionModel(10))

    @test GO.manifold(relate) isa GO.Planar
    @test GO.manifold(overlay) isa GO.Planar
    @test GO.rebuild(relate, GO.Planar()) === relate
    @test GO.rebuild(overlay, GO.Planar()) === overlay
    @test relate.boundary_node_rule isa GO.Mod2BoundaryNodeRule
    @test !relate.prepared
    @test !overlay.strict
    @test !overlay.area_result_only
    @test overlay.optimized
    @test overlay.precision_model isa GO.NoPrecisionModel
    @test fixed_overlay.precision_model isa GO.FixedPrecisionModel
    @test fixed_overlay.precision_model.scale == 10.0
end

@testset "Topology vocabulary" begin
    @test GO.location_index(GO.loc_interior) == 1
    @test GO.location_index(GO.loc_boundary) == 2
    @test GO.location_index(GO.loc_exterior) == 3

    @test GO.dimension_char(GO.dim_false) == 'F'
    @test GO.dimension_char(GO.dim_point) == '0'
    @test GO.dimension_char(GO.dim_line) == '1'
    @test GO.dimension_char(GO.dim_area) == '2'
    @test GO.dimension_from_char('f') == GO.dim_false
    @test GO.max_dimension(GO.dim_point, GO.dim_area) == GO.dim_area

    @test !GO.is_in_boundary(GO.Mod2BoundaryNodeRule(), 2)
    @test GO.is_in_boundary(GO.Mod2BoundaryNodeRule(), 3)
    @test GO.is_in_boundary(GO.EndpointBoundaryNodeRule(), 1)
    @test GO.is_in_boundary(GO.MultivalentEndpointBoundaryNodeRule(), 2)
    @test GO.is_in_boundary(GO.MonovalentEndpointBoundaryNodeRule(), 1)
end

@testset "IntersectionMatrix" begin
    matrix = GO.IntersectionMatrix()
    @test GO.de9im_string(matrix) == "FFFFFFFFF"
    @test matrix[GO.loc_interior, GO.loc_interior] == GO.dim_false

    matrix[GO.loc_interior, GO.loc_interior] = GO.dim_point
    @test GO.de9im_string(matrix) == "0FFFFFFFF"
    @test GO.matches(matrix, "T********")
    @test !GO.matches(matrix, "F********")

    GO.set_at_least!(matrix, GO.loc_interior, GO.loc_interior, GO.dim_line)
    @test matrix[GO.loc_interior, GO.loc_interior] == GO.dim_line
    GO.set_at_least!(matrix, GO.loc_interior, GO.loc_interior, GO.dim_point)
    @test matrix[GO.loc_interior, GO.loc_interior] == GO.dim_line

    other = GO.IntersectionMatrix("F1FFFFFFF")
    GO.set_at_least!(matrix, other)
    @test matrix[GO.loc_interior, GO.loc_boundary] == GO.dim_line
    @test sprint(show, matrix) == GO.de9im_string(matrix)
end

@testset "NG extraction" begin
    line = GI.LineString([(0.0, 0.0), (1.0, 1.0), (1.0, 1.0), (2.0, 2.0)])
    segments = GO.extract_ng_segment_strings(line, Float32; input_side = GO.input_b)
    @test length(segments) == 1

    segment = only(segments)
    @test segment.points isa Vector{Tuple{Float32,Float32}}
    @test segment.points == [(0.0f0, 0.0f0), (1.0f0, 1.0f0), (2.0f0, 2.0f0)]
    @test segment.had_repeated_coordinates
    @test !segment.is_zero_length
    @test segment.source.input_side == GO.input_b
    @test segment.source.source_dimension == GO.dim_line
    @test segment.source.element_id == 1
    @test segment.source.ring_id == GO.NG_NO_RING_ID
    @test segment.source.ring_role == GO.ring_none
    @test segment.source.source_orientation == GO.ring_orientation_none
    @test segment.source.depth_delta == 0
    @test segment.source.parent_polygonal === nothing

    zero_length = GI.LineString([(5.0, 5.0), (5.0, 5.0), (5.0, 5.0)])
    zero_segments = GO.extract_ng_segment_strings(zero_length)
    @test length(zero_segments) == 1
    @test only(zero_segments).points == [(5.0, 5.0)]
    @test only(zero_segments).is_zero_length

    extent = Extents.Extent(X = (9.0, 11.0), Y = (9.0, 11.0))
    multiline = GI.MultiLineString([
        [(0.0, 0.0), (1.0, 1.0)],
        [(10.0, 10.0), (11.0, 11.0)],
    ])
    clipped_segments = GO.extract_ng_segment_strings(multiline; extent)
    @test length(clipped_segments) == 1
    @test only(clipped_segments).points == [(10.0, 10.0), (11.0, 11.0)]
    @test only(clipped_segments).source.element_id == 1
end

@testset "NG polygon ring provenance" begin
    polygon = GI.Polygon([
        [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)],
        [(1.0, 1.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0), (1.0, 1.0)],
    ])

    segments = GO.extract_ng_segment_strings(polygon; orient_rings = :relateng)
    @test length(segments) == 2
    shell, hole = segments

    @test shell.source.source_dimension == GO.dim_area
    @test shell.source.element_id == 1
    @test hole.source.element_id == 1
    @test shell.source.ring_id == 0
    @test hole.source.ring_id == 1
    @test shell.source.ring_role == GO.ring_shell
    @test hole.source.ring_role == GO.ring_hole
    @test shell.source.source_orientation == GO.ring_counterclockwise
    @test hole.source.source_orientation == GO.ring_clockwise
    @test shell.source.geometry === polygon
    @test hole.source.geometry === polygon
    @test shell.source.parent_polygonal === polygon
    @test hole.source.parent_polygonal === polygon
    @test shell.source.coordinates_reversed
    @test hole.source.coordinates_reversed
    @test GO.ng_is_clockwise(shell.points)
    @test !GO.ng_is_clockwise(hole.points)
    @test shell.source.depth_delta == -1
    @test hole.source.depth_delta == -1

    multipolygon = GI.MultiPolygon([polygon])
    mp_segments = GO.extract_ng_segment_strings(multipolygon)
    @test length(mp_segments) == 2
    @test all(segment -> segment.source.parent_polygonal === multipolygon, mp_segments)
end

@testset "NG point extraction and dimensions" begin
    line = GI.LineString([(0.0, 0.0), (1.0, 1.0)])
    polygon = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 0.0)]])
    collection = GI.GeometryCollection([
        GI.Point(1.0, 2.0),
        GI.MultiPoint([(3.0, 4.0), (5.0, 6.0)]),
        line,
    ])

    points = GO.extract_ng_points(collection; input_side = GO.input_b)
    @test getproperty.(points, :point) == [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)]
    @test all(point -> point.source.input_side == GO.input_b, points)
    @test all(point -> point.source.source_dimension == GO.dim_point, points)
    @test getproperty.(getproperty.(points, :source), :element_id) == [1, 2, 3]

    @test GO.ng_source_dimension(GI.Point(1.0, 1.0)) == GO.dim_point
    @test GO.ng_source_dimension(line) == GO.dim_line
    @test GO.ng_source_dimension(polygon) == GO.dim_area
    @test GO.ng_source_dimension(collection) == GO.dim_line
end

@testset "NG segment primitives" begin
    @test GO.ng_orientation((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)) > 0
    @test GO.ng_cross((1.0, 0.0), (0.0, 1.0)) > 0

    cross = GO.ng_segment_intersection(
        ((0.0, 0.0), (2.0, 2.0)),
        ((0.0, 2.0), (2.0, 0.0)),
    )
    @test cross.orientation == GO.line_cross
    @test GO.ng_has_intersection(cross)
    @test all(cross.point1 .≈ (1.0, 1.0))
    @test all(cross.fraction1 .≈ (0.5, 0.5))
    @test GO.ng_intersection_points(cross) == [cross.point1]

    hinge = GO.ng_segment_intersection(
        ((0.0, 0.0), (1.0, 1.0)),
        ((1.0, 1.0), (2.0, 0.0)),
    )
    @test hinge.orientation == GO.line_hinge
    @test hinge.point1 == (1.0, 1.0)

    overlap = GO.ng_segment_intersection(
        ((2.0, 0.0), (0.0, 0.0)),
        ((0.5, 0.0), (1.5, 0.0)),
    )
    @test overlap.orientation == GO.line_over
    @test overlap.fraction1[1] <= overlap.fraction2[1]
    @test GO.ng_intersection_points(overlap) == [(1.5, 0.0), (0.5, 0.0)]

    degenerate = GO.ng_segment_intersection(
        ((0.0, 0.0), (0.0, 0.0)),
        ((0.0, 0.0), (1.0, 1.0)),
    )
    @test degenerate.orientation == GO.line_out
    @test !GO.ng_has_intersection(degenerate)
    @test degenerate.is_degenerate_a
    @test !degenerate.is_degenerate_b

    @test GO.ng_segments_maybe_intersect(
        ((0.0, 0.0), (1.0, 1.0)),
        ((0.5, 0.5), (2.0, 2.0)),
    )
    @test !GO.ng_segments_maybe_intersect(
        ((0.0, 0.0), (1.0, 1.0)),
        ((2.0, 2.0), (3.0, 3.0)),
    )
end

@testset "NG precision hooks and half-closed ownership" begin
    snapped = GO.ng_segment_intersection(
        ((0.1, 0.1), (1.9, 1.9)),
        ((0.1, 1.9), (1.9, 0.1));
        precision_model = TestGridPrecision(),
    )
    @test snapped.orientation == GO.line_cross
    @test snapped.point1 == (1.0, 1.0)

    fixed = GO.ng_segment_intersection(
        ((0.0, 11.0), (620.0, 10.0)),
        ((400.0, 60.0), (400.0, 10.0));
        precision_model = GO.FixedPrecisionModel(1.0),
    )
    @test fixed.orientation == GO.line_cross
    @test fixed.point1 == (400.0, 10.0)

    open_string = only(GO.extract_ng_segment_strings(
        GI.LineString([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]),
    ))
    @test GO.ng_is_containing_segment(open_string, 1, (0.0, 0.0))
    @test !GO.ng_is_containing_segment(open_string, 1, (1.0, 0.0))
    @test GO.ng_is_containing_segment(open_string, 2, (1.0, 0.0))
    @test GO.ng_is_containing_segment(open_string, 2, (2.0, 0.0))
    @test GO.ng_is_containing_segment(open_string, 1, (0.5, 0.0))

    closed_points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]
    @test GO.ng_is_closed_segment_string(closed_points)
    @test !GO.ng_is_containing_segment(closed_points, 1, (1.0, 0.0))
    @test GO.ng_is_containing_segment(closed_points, 2, (1.0, 0.0))
    @test !GO.ng_is_containing_segment(closed_points, 3, (0.0, 0.0))
    @test_throws BoundsError GO.ng_is_containing_segment(closed_points, 4, (0.0, 0.0))
end
