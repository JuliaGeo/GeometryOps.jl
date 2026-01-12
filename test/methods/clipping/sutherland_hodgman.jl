using Test
import GeometryOps as GO
import GeoInterface as GI

@testset "ConvexConvexSutherlandHodgman" begin
    @testset "Basic intersection" begin
        # Two overlapping squares - intersection is 1x1 square
        square1 = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
        square2 = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square1, square2)
        @test result isa GI.Polygon
        @test GO.area(result) ≈ 1.0 atol=1e-10
    end

    @testset "No intersection" begin
        # Disjoint squares
        square1 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        square2 = GI.Polygon([[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 6.0), (5.0, 5.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square1, square2)
        @test result isa GI.Polygon
        @test GO.area(result) ≈ 0.0 atol=1e-10
    end

    @testset "One contains other" begin
        # Large square contains small square
        large = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
        small = GI.Polygon([[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0), (2.0, 2.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), large, small)
        @test result isa GI.Polygon
        @test GO.area(result) ≈ 4.0 atol=1e-10

        # Reverse order should give same result
        result2 = GO.intersection(GO.ConvexConvexSutherlandHodgman(), small, large)
        @test GO.area(result2) ≈ 4.0 atol=1e-10
    end

    @testset "Triangles" begin
        # Two overlapping triangles (both CCW winding)
        tri1 = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (2.0, 4.0), (0.0, 0.0)]])
        tri2 = GI.Polygon([[(0.0, 2.0), (2.0, -2.0), (4.0, 2.0), (0.0, 2.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), tri1, tri2)
        @test result isa GI.Polygon
        @test GO.area(result) > 0
    end

    @testset "Identical polygons" begin
        # Same polygon should return itself
        square = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square, square)
        @test result isa GI.Polygon
        @test GO.area(result) ≈ 4.0 atol=1e-10
    end

    @testset "Shared edge" begin
        # Two squares sharing an edge
        square1 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        square2 = GI.Polygon([[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 0.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square1, square2)
        @test result isa GI.Polygon
        # Shared edge only - area should be 0 or near 0
        @test GO.area(result) ≈ 0.0 atol=1e-10
    end

    @testset "Unsupported geometry types" begin
        square = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        point = GI.Point(0.5, 0.5)

        @test_throws ArgumentError GO.intersection(GO.ConvexConvexSutherlandHodgman(), square, point)
        @test_throws ArgumentError GO.intersection(GO.ConvexConvexSutherlandHodgman(), point, square)
    end
end
