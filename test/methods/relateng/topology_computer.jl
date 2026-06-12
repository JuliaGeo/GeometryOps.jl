# Tests for the RelateNG TopologyComputer (topology_computer.jl), the port
# of JTS TopologyComputer.java. JTS has no dedicated JUnit file for it; per
# the implementation plan (Task 18) the computer is tested through its public
# entry points with a `RelateMatrixPredicate` attached, asserting the
# resulting (partial) IM strings:
#
# - `init_exterior_dims!` a-priori exterior entries for all dim pairs
#   (hand-derived from TopologyComputer.java:44-102),
# - empty-geometry initialization,
# - the addX entry points for every target dimension,
# - `add_intersection!` + `evaluate_nodes!` (symbolic node-section grouping),
# - `updateAreaAreaCross` via `rk_is_crossing`,
# - the D3 coincidence-merge pass (including the multi-segment-pair
#   crossing wheel, which exercises the exact rational apex fallback in
#   `rk_compare_edge_dir`),
# - short-circuiting once the predicate value is known.

using Test
import GeometryOps as GO
import GeometryOps: Planar, True, False
import GeoInterface as GI
import LibGEOS as LG  # only for EMPTY geometries — GI wrappers cannot be empty

const M = Planar()
const EX = True()

rgeom(g) = GO.RelateGeometry(M, g; exact = EX)

# A TopologyComputer with a full-matrix predicate attached; returns both.
function im_computer(rga::GO.RelateGeometry, rgb::GO.RelateGeometry)
    pred = GO.RelateMatrixPredicate()
    return GO.TopologyComputer(pred, rga, rgb), pred
end
im_computer(ga, gb) = im_computer(rgeom(ga), rgeom(gb))

imstr(pred) = string(GO.result_im(pred))

# The partial IM accumulated by `init_exterior_dims!` alone.
init_im(ga, gb) = imstr(im_computer(ga, gb)[2])

# 1-based index of the segment (p, q) (either orientation) in a segment string.
function find_seg(ss, p, q)
    pts = ss.pts
    for i in 1:(length(pts) - 1)
        (pts[i] == p && pts[i + 1] == q) && return i
        (pts[i] == q && pts[i + 1] == p) && return i
    end
    error("segment not found")
end

# Fixtures
const PT_A = GI.Point(1.0, 1.0)
const PT_B = GI.Point(5.0, 5.0)
const LINE_A = GI.LineString([(0.0, 0.0), (4.0, 0.0)])
const LINE_B = GI.LineString([(0.0, 5.0), (4.0, 5.0)])
const RING_LINE = GI.LineString([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 0.0)])
const POLY_A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
const POLY_B = GI.Polygon([[(5.0, 5.0), (7.0, 5.0), (7.0, 7.0), (5.0, 7.0), (5.0, 5.0)]])
const ZERO_LEN_LINE = GI.LineString([(3.0, 3.0), (3.0, 3.0)])

@testset "init_exterior_dims!: a-priori exterior entries" begin
    # Derived by hand from TopologyComputer.java:44-83. The matrix starts
    # all-F except E/E = 2 (IMPredicate constructor).
    # P/P: no a-priori entries
    @test init_im(PT_A, PT_B) == "FFFFFFFF2"
    # P/L: P exterior intersects L interior
    @test init_im(PT_A, LINE_B) == "FFFFFF1F2"
    @test init_im(LINE_A, PT_B) == "FF1FFFFF2"
    # P/A: the area Int and Bdy intersect the point exterior
    @test init_im(PT_A, POLY_B) == "FFFFFF212"
    @test init_im(POLY_A, PT_B) == "FF2FF1FF2"
    # L/A: the area interior intersects the line exterior
    @test init_im(LINE_A, POLY_B) == "FFFFFF2F2"
    @test init_im(POLY_A, LINE_B) == "FF2FFFFF2"
    # L/L and A/A: no a-priori entries
    @test init_im(LINE_A, LINE_B) == "FFFFFFFF2"
    @test init_im(POLY_A, POLY_B) == "FFFFFFFF2"
    # zero-length lines have real dimension P (getDimensionReal)
    @test init_im(ZERO_LEN_LINE, PT_B) == "FFFFFFFF2"
end

@testset "init_exterior_dims!: empty geometries" begin
    # Derived by hand from TopologyComputer.java:74-102 (initExteriorEmpty).
    empty_a = LG.readgeom("POLYGON EMPTY")
    # empty vs P: non-empty interior intersects the exterior (dim P)
    @test init_im(empty_a, PT_B) == "FFFFFF0F2"
    @test init_im(PT_A, empty_a) == "FF0FFFFF2"
    # empty vs L (open line has a boundary)
    @test init_im(empty_a, LINE_B) == "FFFFFF102"
    @test init_im(LINE_A, empty_a) == "FF1FF0FF2"
    # empty vs closed line (ring): no boundary entry
    @test init_im(empty_a, RING_LINE) == "FFFFFF1F2"
    # empty vs A
    @test init_im(empty_a, POLY_B) == "FFFFFF212"
    @test init_im(POLY_A, empty_a) == "FF2FF1FF2"
    # empty vs empty: nothing
    @test init_im(empty_a, LG.readgeom("LINESTRING EMPTY")) == "FFFFFFFF2"
end

@testset "accessors and flags" begin
    tc, _ = im_computer(POLY_A, POLY_B)
    @test GO.get_geometry(tc, GO.GEOM_A) === tc.geom_a
    @test GO.get_geometry(tc, GO.GEOM_B) === tc.geom_b
    @test GO.get_dimension(tc, true) == GO.DIM_A
    @test GO.get_dimension(tc, false) == GO.DIM_A
    @test GO.is_area_area(tc)
    @test !GO.is_area_area(im_computer(LINE_A, POLY_B)[1])
    # exterior-check requirement forwards to the predicate
    @test GO.is_exterior_check_required(tc, true)
    @test GO.is_exterior_check_required(tc, false)

    # isSelfNodingRequired (TopologyComputer.java:142-154)
    # - predicate does not require self-noding
    pred = GO.pred_intersects()
    @test !GO.is_self_noding_required(GO.TopologyComputer(pred, rgeom(LINE_A), rgeom(LINE_B)))
    # - A requires self-noding (a line may self-cross)
    @test GO.is_self_noding_required(im_computer(LINE_A, POLY_B)[1])
    # - polygonal A and B never require self-noding
    @test !GO.is_self_noding_required(im_computer(POLY_A, POLY_B)[1])
    # - B a mixed A/L GC requires full noding
    mixed_b = GI.GeometryCollection([POLY_B, LINE_B])
    @test GO.is_self_noding_required(im_computer(POLY_A, mixed_b)[1])

    # the manifold/exact settings of both inputs must agree
    @test_throws ArgumentError GO.TopologyComputer(GO.RelateMatrixPredicate(),
        GO.RelateGeometry(M, PT_A; exact = True()),
        GO.RelateGeometry(M, PT_B; exact = False()))
end

@testset "add_point_on_point_interior!/_exterior!" begin
    tc, pred = im_computer(PT_A, PT_B)
    GO.add_point_on_point_interior!(tc, (1.0, 1.0))
    @test imstr(pred) == "0FFFFFFF2"

    tc, pred = im_computer(PT_A, PT_B)
    GO.add_point_on_point_exterior!(tc, GO.GEOM_A, (1.0, 1.0))
    @test imstr(pred) == "FF0FFFFF2"

    tc, pred = im_computer(PT_A, PT_B)
    GO.add_point_on_point_exterior!(tc, GO.GEOM_B, (5.0, 5.0))
    @test imstr(pred) == "FFFFFF0F2"
end

@testset "add_point_on_geometry!" begin
    # target dim P: only the point-interior entry
    tc, pred = im_computer(PT_A, PT_B)
    GO.add_point_on_geometry!(tc, true, GO.LOC_INTERIOR, GO.DIM_P, (5.0, 5.0))
    @test imstr(pred) == "0FFFFFFF2"

    # target dim L: no extra inference (zero-length-line caveat in Java)
    tc, pred = im_computer(PT_A, LINE_B)
    GO.add_point_on_geometry!(tc, true, GO.LOC_INTERIOR, GO.DIM_L, (1.0, 5.0))
    @test imstr(pred) == "0FFFFF1F2"

    # target dim A: area interior and boundary extend beyond the point
    tc, pred = im_computer(PT_A, POLY_B)
    GO.add_point_on_geometry!(tc, true, GO.LOC_INTERIOR, GO.DIM_A, (6.0, 6.0))
    @test imstr(pred) == "0FFFFF212"

    # B point in A area: entries transposed
    tc, pred = im_computer(POLY_A, PT_B)
    GO.add_point_on_geometry!(tc, false, GO.LOC_INTERIOR, GO.DIM_A, (1.0, 1.0))
    @test imstr(pred) == "0F2FF1FF2"

    # empty target: no entries to infer, the dim switch is never reached
    tc, pred = im_computer(PT_A, LG.readgeom("POLYGON EMPTY"))
    GO.add_point_on_geometry!(tc, true, GO.LOC_EXTERIOR, GO.DIM_FALSE, (1.0, 1.0))
    @test imstr(pred) == "FF0FFFFF2"

    # unknown target dimension throws (Java IllegalStateException)
    tc, _ = im_computer(PT_A, POLY_B)
    @test_throws ArgumentError GO.add_point_on_geometry!(tc, true, GO.LOC_INTERIOR, 5, (6.0, 6.0))
end

@testset "add_line_end_on_geometry!" begin
    # target dim P: only the line-end entry
    tc, pred = im_computer(LINE_A, PT_B)
    GO.add_line_end_on_geometry!(tc, true, GO.LOC_BOUNDARY, GO.LOC_INTERIOR, GO.DIM_P, (0.0, 0.0))
    @test imstr(pred) == "FF10FFFF2"

    # target dim L, line end in line EXTERIOR: source line interior extends
    # into the target exterior
    tc, pred = im_computer(LINE_A, LINE_B)
    GO.add_line_end_on_geometry!(tc, true, GO.LOC_BOUNDARY, GO.LOC_EXTERIOR, GO.DIM_L, (0.0, 0.0))
    @test imstr(pred) == "FF1FF0FF2"

    # target dim L, line end on line INTERIOR: only the end-point entry
    tc, pred = im_computer(LINE_A, LINE_B)
    GO.add_line_end_on_geometry!(tc, true, GO.LOC_BOUNDARY, GO.LOC_INTERIOR, GO.DIM_L, (0.0, 0.0))
    @test imstr(pred) == "FFF0FFFF2"

    # target dim A, line end in area INTERIOR
    tc, pred = im_computer(LINE_A, POLY_B)
    GO.add_line_end_on_geometry!(tc, true, GO.LOC_BOUNDARY, GO.LOC_INTERIOR, GO.DIM_A, (6.0, 6.0))
    @test imstr(pred) == "1FF0FF2F2"

    # target dim A, line end on area BOUNDARY: no further inference
    tc, pred = im_computer(LINE_A, POLY_B)
    GO.add_line_end_on_geometry!(tc, true, GO.LOC_BOUNDARY, GO.LOC_BOUNDARY, GO.DIM_A, (5.0, 5.0))
    @test imstr(pred) == "FFFF0F2F2"

    # target dim A, line end in area EXTERIOR
    tc, pred = im_computer(LINE_A, POLY_B)
    GO.add_line_end_on_geometry!(tc, true, GO.LOC_INTERIOR, GO.LOC_EXTERIOR, GO.DIM_A, (1.0, 0.0))
    @test imstr(pred) == "FF1FFF2F2"

    # empty target: only the end-point entry, dim switch never reached
    tc, pred = im_computer(LINE_A, LG.readgeom("POLYGON EMPTY"))
    GO.add_line_end_on_geometry!(tc, true, GO.LOC_BOUNDARY, GO.LOC_EXTERIOR, GO.DIM_FALSE, (0.0, 0.0))
    @test imstr(pred) == "FF1FF0FF2"

    # unknown target dimension throws
    tc, _ = im_computer(LINE_A, POLY_B)
    @test_throws ArgumentError GO.add_line_end_on_geometry!(tc, true, GO.LOC_BOUNDARY, GO.LOC_INTERIOR, 5, (6.0, 6.0))
end

@testset "add_area_vertex!" begin
    # locTarget EXTERIOR, boundary vertex: full neighbourhood inference
    tc, pred = im_computer(POLY_A, LINE_B)
    GO.add_area_vertex!(tc, true, GO.LOC_BOUNDARY, GO.LOC_EXTERIOR, GO.DIM_L, (0.0, 0.0))
    @test imstr(pred) == "FF2FF1FF2"

    # locTarget EXTERIOR, interior vertex (overlapping-GC case): only Int/Ext
    tc, pred = im_computer(POLY_A, LINE_B)
    GO.add_area_vertex!(tc, true, GO.LOC_INTERIOR, GO.LOC_EXTERIOR, GO.DIM_L, (1.0, 1.0))
    @test imstr(pred) == "FF2FFFFF2"

    # on point (addAreaVertexOnPoint): boundary vertex
    tc, pred = im_computer(POLY_A, PT_B)
    GO.add_area_vertex!(tc, true, GO.LOC_BOUNDARY, GO.LOC_INTERIOR, GO.DIM_P, (0.0, 0.0))
    @test imstr(pred) == "FF20F1FF2"

    # on line (addAreaVertexOnLine): boundary vertex — only the point entry
    tc, pred = im_computer(POLY_A, LINE_B)
    GO.add_area_vertex!(tc, true, GO.LOC_BOUNDARY, GO.LOC_INTERIOR, GO.DIM_L, (0.0, 0.0))
    @test imstr(pred) == "FF20FFFF2"

    # on line, interior vertex: area interior beyond the line
    tc, pred = im_computer(POLY_A, LINE_B)
    GO.add_area_vertex!(tc, true, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_L, (1.0, 1.0))
    @test imstr(pred) == "0F2FFFFF2"

    # on area (addAreaVertexOnArea):
    # B/B: deferred to node analysis (point entry only)
    tc, pred = im_computer(POLY_A, POLY_B)
    GO.add_area_vertex!(tc, true, GO.LOC_BOUNDARY, GO.LOC_BOUNDARY, GO.DIM_A, (0.0, 0.0))
    @test imstr(pred) == "FFFF0FFF2"
    # I/B
    tc, pred = im_computer(POLY_A, POLY_B)
    GO.add_area_vertex!(tc, true, GO.LOC_INTERIOR, GO.LOC_BOUNDARY, GO.DIM_A, (1.0, 1.0))
    @test imstr(pred) == "212FFFFF2"
    # B/I
    tc, pred = im_computer(POLY_A, POLY_B)
    GO.add_area_vertex!(tc, true, GO.LOC_BOUNDARY, GO.LOC_INTERIOR, GO.DIM_A, (0.0, 0.0))
    @test imstr(pred) == "2FF1FF2F2"
    # I/I
    tc, pred = im_computer(POLY_A, POLY_B)
    GO.add_area_vertex!(tc, true, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_A, (1.0, 1.0))
    @test imstr(pred) == "2FFFFFFF2"

    # unknown target dimension throws
    tc, _ = im_computer(POLY_A, POLY_B)
    @test_throws ArgumentError GO.add_area_vertex!(tc, true, GO.LOC_BOUNDARY, GO.LOC_INTERIOR, 5, (0.0, 0.0))
end

@testset "add_intersection! + evaluate_nodes!: L/L proper crossing" begin
    line_a = GI.LineString([(-1.0, 0.0), (1.0, 0.0)])
    line_b = GI.LineString([(0.0, -1.0), (0.0, 1.0)])
    rga, rgb = rgeom(line_a), rgeom(line_b)
    tc, pred = im_computer(rga, rgb)
    ssa = GO.extract_segment_strings(rga, true, nothing)[1]
    ssb = GO.extract_segment_strings(rgb, false, nothing)[1]
    key = GO.crossing_node((-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0))
    GO.add_intersection!(tc, GO.create_node_section(ssa, 1, key),
        GO.create_node_section(ssb, 1, key))
    # updateNodeLocation: the crossing lies in both line interiors
    @test imstr(pred) == "0FFFFFFF2"
    @test length(tc.node_sections) == 1
    GO.evaluate_nodes!(tc)
    # node wheel: each line's edges lie in the other line's exterior
    @test imstr(pred) == "0F1FFF1F2"
end

@testset "add_intersection! + evaluate_nodes!: A/A proper crossing" begin
    # Overlapping unit-offset squares; one boundary crossing node carries the
    # entire standard overlap matrix.
    sq_a = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    sq_b = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    rga, rgb = rgeom(sq_a), rgeom(sq_b)
    tc, pred = im_computer(rga, rgb)
    ssa = GO.extract_segment_strings(rga, true, nothing)[1]
    ssb = GO.extract_segment_strings(rgb, false, nothing)[1]
    # A's top edge crosses B's left edge properly at (1, 2)
    ia = find_seg(ssa, (0.0, 2.0), (2.0, 2.0))
    ib = find_seg(ssb, (1.0, 1.0), (1.0, 3.0))
    key = GO.crossing_node((0.0, 2.0), (2.0, 2.0), (1.0, 1.0), (1.0, 3.0))
    GO.add_intersection!(tc, GO.create_node_section(ssa, ia, key),
        GO.create_node_section(ssb, ib, key))
    # proper crossing of two area boundaries: interiors intersect (dim 2),
    # and the node lies on both boundaries (dim 0)
    @test imstr(pred) == "2FFF0FFF2"
    GO.evaluate_nodes!(tc)
    @test imstr(pred) == "212101212"
end

@testset "updateAreaAreaCross via rk_is_crossing (vertex node)" begin
    # Two polygons whose boundaries cross at a shared vertex (origin).
    # A's corner directions: (-2,-1) and (1,2); B's: (1,-2) and (-1,2) —
    # interleaved around the node, so the boundaries cross.
    poly_a = GI.Polygon([[(-2.0, -1.0), (0.0, 0.0), (1.0, 2.0), (3.0, 2.0), (3.0, -1.0), (-2.0, -1.0)]])
    poly_b = GI.Polygon([[(1.0, -2.0), (0.0, 0.0), (-1.0, 2.0), (-3.0, 2.0), (-3.0, -2.0), (1.0, -2.0)]])
    vkey = GO.vertex_node((0.0, 0.0))
    nsa = GO.NodeSection(true, Int8(2), Int32(1), Int32(0), nothing, true, (-2.0, -1.0), vkey, (1.0, 2.0))
    nsb = GO.NodeSection(false, Int8(2), Int32(1), Int32(0), nothing, true, (1.0, -2.0), vkey, (-1.0, 2.0))
    tc, pred = im_computer(poly_a, poly_b)
    GO.add_intersection!(tc, nsa, nsb)
    # crossing: interiors intersect; node is on both boundaries
    @test imstr(pred) == "2FFF0FFF2"

    # touching variant: B's corner directions both on the same side of A's
    # corner — boundaries do not cross, interiors entry is not inferred
    nsb_touch = GO.NodeSection(false, Int8(2), Int32(1), Int32(0), nothing, true, (-1.0, 2.0), vkey, (-2.0, 1.0))
    tc, pred = im_computer(poly_a, poly_b)
    GO.add_intersection!(tc, nsa, nsb_touch)
    @test imstr(pred) == "FFFF0FFF2"
end

@testset "rk_compare_edge_dir: exact rational apex fallback" begin
    # Direction points which are NOT endpoints of the crossing node's
    # defining segments (they arise from D3 coincidence-merged nodes).
    key = GO.crossing_node((-1.0, -1.0), (1.0, 1.0), (-1.0, 0.0), (1.0, 0.0))  # apex (0,0)
    @test GO.rk_compare_edge_dir(M, key, (-1.0, 1.0), (1.0, -1.0); exact = EX) < 0   # 135° vs 315°
    @test GO.rk_compare_edge_dir(M, key, (1.0, -1.0), (-1.0, 1.0); exact = EX) > 0
    # mixed defining/foreign directions
    @test GO.rk_compare_edge_dir(M, key, (1.0, 1.0), (-1.0, 1.0); exact = EX) < 0    # 45° vs 135°
    @test GO.rk_compare_edge_dir(M, key, (-1.0, 1.0), (1.0, 1.0); exact = EX) > 0
    # equal angles: a foreign direction on the same ray as a defining endpoint
    @test GO.rk_compare_edge_dir(M, key, (2.0, 2.0), (1.0, 1.0); exact = EX) == 0

    # apex not representable in Float64: (2/3, 2/3)
    key2 = GO.crossing_node((0.0, 0.0), (3.0, 3.0), (0.0, 1.0), (2.0, 0.0))
    @test GO.rk_compare_edge_dir(M, key2, (1.0, 0.5), (0.0, 0.5); exact = EX) > 0    # SE vs SW quadrant
    @test GO.rk_compare_edge_dir(M, key2, (0.0, 0.5), (1.0, 0.5); exact = EX) < 0
    # exact collinearity through the irrational apex
    @test GO.rk_compare_edge_dir(M, key2, (1.0, 1.0), (3.0, 3.0); exact = EX) == 0
end

@testset "D3 coincidence merge: crossings from different segment pairs" begin
    # A self-crosses at the origin; B passes through the same point with a
    # third direction. Three distinct crossing keys denote one node.
    line1 = GI.LineString([(-1.0, -1.0), (1.0, 1.0)])
    line2 = GI.LineString([(-1.0, 1.0), (1.0, -1.0)])
    ml_a = GI.MultiLineString([line1, line2])
    line_b = GI.LineString([(-1.0, 0.0), (1.0, 0.0)])
    rga, rgb = rgeom(ml_a), rgeom(line_b)
    tc, pred = im_computer(rga, rgb)
    @test GO.is_self_noding_required(tc)

    ssa1, ssa2 = GO.extract_segment_strings(rga, true, nothing)
    ssb = GO.extract_segment_strings(rgb, false, nothing)[1]
    k1 = GO.crossing_node((-1.0, -1.0), (1.0, 1.0), (-1.0, 0.0), (1.0, 0.0))
    k2 = GO.crossing_node((-1.0, 1.0), (1.0, -1.0), (-1.0, 0.0), (1.0, 0.0))
    k3 = GO.crossing_node((-1.0, -1.0), (1.0, 1.0), (-1.0, 1.0), (1.0, -1.0))
    @test length(Set([k1, k2, k3])) == 3

    GO.add_intersection!(tc, GO.create_node_section(ssa1, 1, k1), GO.create_node_section(ssb, 1, k1))
    GO.add_intersection!(tc, GO.create_node_section(ssa2, 1, k2), GO.create_node_section(ssb, 1, k2))
    # A self-intersection: contributes sections but no direct IM update
    GO.add_intersection!(tc, GO.create_node_section(ssa1, 1, k3), GO.create_node_section(ssa2, 1, k3))
    @test imstr(pred) == "0FFFFFFF2"
    @test length(tc.node_sections) == 3

    # evaluate_nodes! merges the coinciding keys into one node, then builds
    # a single 6-direction wheel. Directions from the non-canonical segment
    # pairs exercise the rational-apex comparison in rk_compare_edge_dir.
    GO.evaluate_nodes!(tc)
    @test length(tc.node_sections) == 1
    @test imstr(pred) == "0F1FFF1F2"
end

@testset "D3 coincidence merge: vertex node absorbs crossing node" begin
    # B has a vertex exactly at A's self-crossing point. The vertex key is
    # preferred as the canonical merged node.
    line1 = GI.LineString([(-1.0, -1.0), (1.0, 1.0)])
    line2 = GI.LineString([(-1.0, 1.0), (1.0, -1.0)])
    ml_a = GI.MultiLineString([line1, line2])
    line_b = GI.LineString([(-1.0, 0.0), (0.0, 0.0), (1.0, 0.0)])
    rga, rgb = rgeom(ml_a), rgeom(line_b)
    tc, pred = im_computer(rga, rgb)

    ssa1, ssa2 = GO.extract_segment_strings(rga, true, nothing)
    ssb = GO.extract_segment_strings(rgb, false, nothing)[1]
    vkey = GO.vertex_node((0.0, 0.0))
    k3 = GO.crossing_node((-1.0, -1.0), (1.0, 1.0), (-1.0, 1.0), (1.0, -1.0))

    GO.add_intersection!(tc, GO.create_node_section(ssa1, 1, vkey), GO.create_node_section(ssb, 1, vkey))
    GO.add_intersection!(tc, GO.create_node_section(ssa1, 1, k3), GO.create_node_section(ssa2, 1, k3))
    @test length(tc.node_sections) == 2

    GO.evaluate_nodes!(tc)
    @test length(tc.node_sections) == 1
    @test haskey(tc.node_sections, vkey)
    @test imstr(pred) == "0F1FFF1F2"
end

@testset "no merge of distinct crossing points" begin
    # The coincidence merge runs whenever crossing keys exist (Task 21:
    # JTS merges via its coordinate-keyed node map in every mode), but it
    # only merges keys denoting the SAME exact point — the two genuine
    # square-boundary crossings here are distinct points and stay distinct.
    sq_a = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    sq_b = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    rga, rgb = rgeom(sq_a), rgeom(sq_b)
    tc, pred = im_computer(rga, rgb)
    @test !GO.is_self_noding_required(tc)
    ssa = GO.extract_segment_strings(rga, true, nothing)[1]
    ssb = GO.extract_segment_strings(rgb, false, nothing)[1]
    # the two genuine boundary crossings of the squares: (1,2) and (2,1)
    ia1 = find_seg(ssa, (0.0, 2.0), (2.0, 2.0))
    ib1 = find_seg(ssb, (1.0, 1.0), (1.0, 3.0))
    key1 = GO.crossing_node((0.0, 2.0), (2.0, 2.0), (1.0, 1.0), (1.0, 3.0))
    ia2 = find_seg(ssa, (2.0, 0.0), (2.0, 2.0))
    ib2 = find_seg(ssb, (1.0, 1.0), (3.0, 1.0))
    key2 = GO.crossing_node((2.0, 0.0), (2.0, 2.0), (1.0, 1.0), (3.0, 1.0))
    GO.add_intersection!(tc, GO.create_node_section(ssa, ia1, key1), GO.create_node_section(ssb, ib1, key1))
    GO.add_intersection!(tc, GO.create_node_section(ssa, ia2, key2), GO.create_node_section(ssb, ib2, key2))
    @test length(tc.node_sections) == 2
    GO.evaluate_nodes!(tc)
    @test length(tc.node_sections) == 2
    @test imstr(pred) == "212101212"
end

@testset "short-circuit: known result freezes the predicate" begin
    pred = GO.pred_intersects()
    tc = GO.TopologyComputer(pred, rgeom(LINE_A), rgeom(LINE_B))
    @test !GO.is_result_known(tc)
    GO.add_point_on_point_interior!(tc, (0.0, 0.0))
    @test GO.is_result_known(tc)
    @test GO.get_result(tc)
    # further updates are no-ops on the determined value
    GO.add_point_on_point_exterior!(tc, GO.GEOM_A, (4.0, 0.0))
    GO.add_area_vertex!(tc, true, GO.LOC_BOUNDARY, GO.LOC_EXTERIOR, GO.DIM_A, (0.0, 0.0))
    GO.evaluate_nodes!(tc)
    @test GO.is_result_known(tc) && GO.get_result(tc)

    # finish! finalizes an undetermined predicate
    pred = GO.pred_intersects()
    tc = GO.TopologyComputer(pred, rgeom(LINE_A), rgeom(LINE_B))
    GO.finish!(tc)
    @test GO.is_result_known(tc)
    @test !GO.get_result(tc)
end

@testset "unknown target dimension throws ArgumentError" begin
    a = GI.Point(0.0, 0.0)
    b = GI.Point(1.0, 1.0)
    rga = GO.RelateGeometry(GO.Planar(), a; exact = GO.True())
    rgb = GO.RelateGeometry(GO.Planar(), b; exact = GO.True())
    tc = GO.TopologyComputer(GO.pred_intersects(), rga, rgb)
    @test_throws ArgumentError GO.add_point_on_geometry!(
        tc, true, GO.LOC_INTERIOR, Int8(9), (0.0, 0.0))
end
