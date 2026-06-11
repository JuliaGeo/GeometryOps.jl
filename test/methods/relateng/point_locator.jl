# Tests for the RelateNG point-location machinery (point_locator.jl).
# LinearBoundary section ports JTS LinearBoundaryTest.java; later tasks
# append AdjacentEdgeLocator and RelatePointLocator tests to this file.

using Test
import GeometryOps as GO
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
