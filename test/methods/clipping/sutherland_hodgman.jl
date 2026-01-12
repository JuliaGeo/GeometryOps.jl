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
end
