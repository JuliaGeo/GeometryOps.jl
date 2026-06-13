# Tests for the RelateNG point-location machinery (point_locator.jl).
# LinearBoundary section ports JTS LinearBoundaryTest.java; the
# AdjacentEdgeLocator section ports JTS AdjacentEdgeLocatorTest.java; the
# RelatePointLocator section ports JTS RelatePointLocatorTest.java.

using Test
import GeometryOps as GO
import GeometryOps: Planar, True
import GeoInterface as GI

# Port of LinearBoundaryTest.checkLinearBoundary. The Java tests parse WKT;
# here each WKT literal is translated by hand into GI.LineStrings and a set
# of expected boundary coordinate tuples (nothing ⇔ Java's null = no boundary).
function check_linear_boundary(lines, rule, expected_boundary)
    lb = GO.LinearBoundary(lines, rule)
    has_boundary_expected = expected_boundary === nothing ? false : true
    @test GO.has_boundary(lb) == has_boundary_expected
    check_boundary_points(lb, lines, expected_boundary)
end

# Port of LinearBoundaryTest.checkBoundaryPoints (+ extractPoints).
function check_boundary_points(lb, lines, expected_boundary)
    bdy_set = expected_boundary === nothing ? Set{Tuple{Float64, Float64}}() :
        Set{Tuple{Float64, Float64}}(expected_boundary)
    for p in bdy_set
        @test GO.is_boundary(lb, p)
    end
    for line in lines, p in GI.getpoint(line)
        pt = (GI.x(p), GI.y(p))
        if !(pt in bdy_set)
            @test !GO.is_boundary(lb, pt)
        end
    end
end

@testset "LinearBoundary" begin
    # testLineMod2: LINESTRING (0 0, 9 9), boundary MULTIPOINT((0 0), (9 9))
    @testset "line Mod2" begin
        lines = [GI.LineString([(0.0, 0.0), (9.0, 9.0)])]
        check_linear_boundary(lines, GO.Mod2Boundary(), [(0.0, 0.0), (9.0, 9.0)])
    end

    # testLines2Mod2: MULTILINESTRING ((0 0, 9 9), (9 9, 5 1)),
    # boundary MULTIPOINT((0 0), (5 1))
    @testset "lines2 Mod2" begin
        lines = [
            GI.LineString([(0.0, 0.0), (9.0, 9.0)]),
            GI.LineString([(9.0, 9.0), (5.0, 1.0)]),
        ]
        check_linear_boundary(lines, GO.Mod2Boundary(), [(0.0, 0.0), (5.0, 1.0)])
    end

    # testLines3Mod2: MULTILINESTRING ((0 0, 9 9), (9 9, 5 1), (9 9, 1 5)),
    # boundary MULTIPOINT((0 0), (5 1), (1 5), (9 9))
    @testset "lines3 Mod2" begin
        lines = [
            GI.LineString([(0.0, 0.0), (9.0, 9.0)]),
            GI.LineString([(9.0, 9.0), (5.0, 1.0)]),
            GI.LineString([(9.0, 9.0), (1.0, 5.0)]),
        ]
        check_linear_boundary(lines, GO.Mod2Boundary(),
            [(0.0, 0.0), (5.0, 1.0), (1.0, 5.0), (9.0, 9.0)])
    end

    # testLines3Monvalent: same lines, MONOVALENT_ENDPOINT_BOUNDARY_RULE,
    # boundary MULTIPOINT((0 0), (5 1), (1 5)) — (9 9) has degree 3, not 1
    @testset "lines3 Monovalent" begin
        lines = [
            GI.LineString([(0.0, 0.0), (9.0, 9.0)]),
            GI.LineString([(9.0, 9.0), (5.0, 1.0)]),
            GI.LineString([(9.0, 9.0), (1.0, 5.0)]),
        ]
        check_linear_boundary(lines, GO.MonovalentEndpointBoundary(),
            [(0.0, 0.0), (5.0, 1.0), (1.0, 5.0)])
    end
end

# Port of AdjacentEdgeLocatorTest.checkLocation. The Java tests parse WKT;
# here each WKT literal is hand-translated into GI geometries.
function check_location(geom, x, y, expected_loc)
    ael = GO.AdjacentEdgeLocator(Planar(), geom; exact = True())
    loc = GO.locate(ael, (Float64(x), Float64(y)))
    @test loc == expected_loc
end

@testset "AdjacentEdgeLocator" begin
    # testAdjacent2:
    # GEOMETRYCOLLECTION (POLYGON ((1 9, 5 9, 5 1, 1 1, 1 9)),
    #                     POLYGON ((9 9, 9 1, 5 1, 5 9, 9 9)))
    @testset "adjacent 2" begin
        gc = GI.GeometryCollection([
            GI.Polygon([[(1.0, 9.0), (5.0, 9.0), (5.0, 1.0), (1.0, 1.0), (1.0, 9.0)]]),
            GI.Polygon([[(9.0, 9.0), (9.0, 1.0), (5.0, 1.0), (5.0, 9.0), (9.0, 9.0)]]),
        ])
        check_location(gc, 5, 5, GO.LOC_INTERIOR)
    end

    # testNonAdjacent:
    # GEOMETRYCOLLECTION (POLYGON ((1 9, 4 9, 5 1, 1 1, 1 9)),
    #                     POLYGON ((9 9, 9 1, 5 1, 5 9, 9 9)))
    @testset "non-adjacent" begin
        gc = GI.GeometryCollection([
            GI.Polygon([[(1.0, 9.0), (4.0, 9.0), (5.0, 1.0), (1.0, 1.0), (1.0, 9.0)]]),
            GI.Polygon([[(9.0, 9.0), (9.0, 1.0), (5.0, 1.0), (5.0, 9.0), (9.0, 9.0)]]),
        ])
        check_location(gc, 5, 5, GO.LOC_BOUNDARY)
    end

    # testAdjacent6WithFilledHoles:
    # GEOMETRYCOLLECTION (
    #   POLYGON ((1 9, 5 9, 6 6, 1 5, 1 9), (2 6, 4 8, 6 6, 2 6)),
    #   POLYGON ((2 6, 4 8, 6 6, 2 6)),
    #   POLYGON ((9 9, 9 5, 6 6, 5 9, 9 9)),
    #   POLYGON ((9 1, 5 1, 6 6, 9 5, 9 1), (7 2, 6 6, 8 3, 7 2)),
    #   POLYGON ((7 2, 6 6, 8 3, 7 2)),
    #   POLYGON ((1 1, 1 5, 6 6, 5 1, 1 1)))
    @testset "adjacent 6 with filled holes" begin
        gc = GI.GeometryCollection([
            GI.Polygon([
                [(1.0, 9.0), (5.0, 9.0), (6.0, 6.0), (1.0, 5.0), (1.0, 9.0)],
                [(2.0, 6.0), (4.0, 8.0), (6.0, 6.0), (2.0, 6.0)],
            ]),
            GI.Polygon([[(2.0, 6.0), (4.0, 8.0), (6.0, 6.0), (2.0, 6.0)]]),
            GI.Polygon([[(9.0, 9.0), (9.0, 5.0), (6.0, 6.0), (5.0, 9.0), (9.0, 9.0)]]),
            GI.Polygon([
                [(9.0, 1.0), (5.0, 1.0), (6.0, 6.0), (9.0, 5.0), (9.0, 1.0)],
                [(7.0, 2.0), (6.0, 6.0), (8.0, 3.0), (7.0, 2.0)],
            ]),
            GI.Polygon([[(7.0, 2.0), (6.0, 6.0), (8.0, 3.0), (7.0, 2.0)]]),
            GI.Polygon([[(1.0, 1.0), (1.0, 5.0), (6.0, 6.0), (5.0, 1.0), (1.0, 1.0)]]),
        ])
        check_location(gc, 6, 6, GO.LOC_INTERIOR)
    end

    # testAdjacent5WithEmptyHole: as above, but the (7 2, 6 6, 8 3) hole is
    # not filled by a matching polygon, so the node is on the boundary.
    @testset "adjacent 5 with empty hole" begin
        gc = GI.GeometryCollection([
            GI.Polygon([
                [(1.0, 9.0), (5.0, 9.0), (6.0, 6.0), (1.0, 5.0), (1.0, 9.0)],
                [(2.0, 6.0), (4.0, 8.0), (6.0, 6.0), (2.0, 6.0)],
            ]),
            GI.Polygon([[(2.0, 6.0), (4.0, 8.0), (6.0, 6.0), (2.0, 6.0)]]),
            GI.Polygon([[(9.0, 9.0), (9.0, 5.0), (6.0, 6.0), (5.0, 9.0), (9.0, 9.0)]]),
            GI.Polygon([
                [(9.0, 1.0), (5.0, 1.0), (6.0, 6.0), (9.0, 5.0), (9.0, 1.0)],
                [(7.0, 2.0), (6.0, 6.0), (8.0, 3.0), (7.0, 2.0)],
            ]),
            GI.Polygon([[(1.0, 1.0), (1.0, 5.0), (6.0, 6.0), (5.0, 1.0), (1.0, 1.0)]]),
        ])
        check_location(gc, 6, 6, GO.LOC_BOUNDARY)
    end

    # testContainedAndAdjacent:
    # GEOMETRYCOLLECTION (POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9)),
    #                     POLYGON ((9 2, 2 2, 2 8, 9 8, 9 2)))
    @testset "contained and adjacent" begin
        gc = GI.GeometryCollection([
            GI.Polygon([[(1.0, 9.0), (9.0, 9.0), (9.0, 1.0), (1.0, 1.0), (1.0, 9.0)]]),
            GI.Polygon([[(9.0, 2.0), (2.0, 2.0), (2.0, 8.0), (9.0, 8.0), (9.0, 2.0)]]),
        ])
        check_location(gc, 9, 5, GO.LOC_BOUNDARY)
        check_location(gc, 9, 8, GO.LOC_BOUNDARY)
    end

    # testDisjointCollinear (bug caused by incorrect point-on-segment logic):
    # GEOMETRYCOLLECTION (MULTIPOLYGON (((1 4, 4 4, 4 1, 1 1, 1 4)),
    #                                   ((5 4, 8 4, 8 1, 5 1, 5 4))))
    @testset "disjoint collinear" begin
        gc = GI.GeometryCollection([
            GI.MultiPolygon([
                GI.Polygon([[(1.0, 4.0), (4.0, 4.0), (4.0, 1.0), (1.0, 1.0), (1.0, 4.0)]]),
                GI.Polygon([[(5.0, 4.0), (8.0, 4.0), (8.0, 1.0), (5.0, 1.0), (5.0, 4.0)]]),
            ]),
        ])
        check_location(gc, 2, 4, GO.LOC_BOUNDARY)
    end
end

# Port of RelatePointLocatorTest.java. The Java fixture WKT
# `gcPLA` is hand-translated into GI constructors:
# GEOMETRYCOLLECTION (POINT (1 1), POINT (2 1), LINESTRING (3 1, 3 9),
#   LINESTRING (4 1, 5 4, 7 1, 4 1), LINESTRING (12 12, 14 14),
#   POLYGON ((6 5, 6 9, 9 9, 9 5, 6 5)),
#   POLYGON ((10 10, 10 16, 16 16, 16 10, 10 10)),
#   POLYGON ((11 11, 11 17, 17 17, 17 11, 11 11)),
#   POLYGON ((12 12, 12 16, 16 16, 16 12, 12 12)))
const gc_PLA = GI.GeometryCollection([
    GI.Point(1.0, 1.0),
    GI.Point(2.0, 1.0),
    GI.LineString([(3.0, 1.0), (3.0, 9.0)]),
    GI.LineString([(4.0, 1.0), (5.0, 4.0), (7.0, 1.0), (4.0, 1.0)]),
    GI.LineString([(12.0, 12.0), (14.0, 14.0)]),
    GI.Polygon([[(6.0, 5.0), (6.0, 9.0), (9.0, 9.0), (9.0, 5.0), (6.0, 5.0)]]),
    GI.Polygon([[(10.0, 10.0), (10.0, 16.0), (16.0, 16.0), (16.0, 10.0), (10.0, 10.0)]]),
    GI.Polygon([[(11.0, 11.0), (11.0, 17.0), (17.0, 17.0), (17.0, 11.0), (11.0, 11.0)]]),
    GI.Polygon([[(12.0, 12.0), (12.0, 16.0), (16.0, 16.0), (16.0, 12.0), (12.0, 12.0)]]),
])

# Port of RelatePointLocatorTest.checkDimLocation.
function check_dim_location(geom, x, y, expected_dim_loc)
    locator = GO.RelatePointLocator(Planar(), geom; exact = True())
    actual = GO.locate_with_dim(locator, (Float64(x), Float64(y)))
    @test actual == expected_dim_loc
end

# Port of RelatePointLocatorTest.checkLineEndDimLocation.
function check_line_end_dim_location(geom, x, y, expected_dim_loc)
    locator = GO.RelatePointLocator(Planar(), geom; exact = True())
    actual = GO.locate_line_end_with_dim(locator, (Float64(x), Float64(y)))
    @test actual == expected_dim_loc
end

# Port of RelatePointLocatorTest.checkNodeLocation.
function check_node_location(geom, x, y, expected_loc)
    locator = GO.RelatePointLocator(Planar(), geom; exact = True())
    actual = GO.locate_node(locator, (Float64(x), Float64(y)), nothing)
    @test actual == expected_loc
end

@testset "RelatePointLocator" begin
    # testPoint
    @testset "point" begin
        check_dim_location(gc_PLA, 1, 1, GO.DL_POINT_INTERIOR)
        check_dim_location(gc_PLA, 0, 1, GO.DL_EXTERIOR)
    end

    # testPointInLine
    @testset "point in line" begin
        check_dim_location(gc_PLA, 3, 8, GO.DL_LINE_INTERIOR)
    end

    # testPointInArea
    @testset "point in area" begin
        check_dim_location(gc_PLA, 8, 8, GO.DL_AREA_INTERIOR)
    end

    # testLine
    @testset "line" begin
        check_dim_location(gc_PLA, 3, 3, GO.DL_LINE_INTERIOR)
        check_dim_location(gc_PLA, 3, 1, GO.DL_LINE_BOUNDARY)
    end

    # testLineInArea
    @testset "line in area" begin
        check_dim_location(gc_PLA, 11, 11, GO.DL_AREA_INTERIOR)
        check_dim_location(gc_PLA, 14, 14, GO.DL_AREA_INTERIOR)
    end

    # testArea
    @testset "area" begin
        check_dim_location(gc_PLA, 8, 8, GO.DL_AREA_INTERIOR)
        check_dim_location(gc_PLA, 9, 9, GO.DL_AREA_BOUNDARY)
    end

    # testAreaInArea
    @testset "area in area" begin
        check_dim_location(gc_PLA, 11, 11, GO.DL_AREA_INTERIOR)
        check_dim_location(gc_PLA, 12, 12, GO.DL_AREA_INTERIOR)
        check_dim_location(gc_PLA, 10, 10, GO.DL_AREA_BOUNDARY)
        check_dim_location(gc_PLA, 16, 16, GO.DL_AREA_INTERIOR)
    end

    # testLineNode
    @testset "line node" begin
        # checkNodeLocation(gcPLA, 12.1, 12.2, Location.INTERIOR) — commented
        # out in the Java test too.
        check_node_location(gc_PLA, 3, 1, GO.LOC_BOUNDARY)
    end

    # testLineEndInGCLA:
    # GEOMETRYCOLLECTION (POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0)),
    #   LINESTRING (12 2, 0 2, 0 5, 5 5), LINESTRING (12 10, 12 2))
    @testset "line end in GC LA" begin
        gc = GI.GeometryCollection([
            GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]]),
            GI.LineString([(12.0, 2.0), (0.0, 2.0), (0.0, 5.0), (5.0, 5.0)]),
            GI.LineString([(12.0, 10.0), (12.0, 2.0)]),
        ])
        check_line_end_dim_location(gc, 5, 5, GO.DL_AREA_INTERIOR)
        check_line_end_dim_location(gc, 12, 2, GO.DL_LINE_INTERIOR)
        check_line_end_dim_location(gc, 12, 10, GO.DL_LINE_BOUNDARY)
    end
end
