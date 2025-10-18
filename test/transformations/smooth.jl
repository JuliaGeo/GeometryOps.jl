using Test
import GeometryOps as GO
import GeoInterface as GI

@testset "Chaikin" begin
    @testset "LineString" begin
        line = GI.LineString([(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)])
        smoothed = GO.smooth(line)
        @test GI.npoint(smoothed) == 6
        @test GI.getpoint(smoothed, 1) == (0.0, 0.0)
        @test GI.getpoint(smoothed, 2) == (0.25, 0.25)
        @test GI.getpoint(smoothed, 3) == (0.75, 0.75)
        @test GI.getpoint(smoothed, 4) == (1.25, 0.75)
        @test GI.getpoint(smoothed, 5) == (1.75, 0.25)
        @test GI.getpoint(smoothed, 6) == (2.0, 0.0)

        smoothed = GO.smooth(line; iterations=2)
        @test GI.npoint(smoothed) == 10
        @test GI.getpoint(smoothed, 1) == (0.0, 0.0)
        @test GI.getpoint(smoothed, 2) == (0.0625, 0.0625)
        @test GI.getpoint(smoothed, 3) == (0.1875, 0.1875)
        @test GI.getpoint(smoothed, 4) == (0.5, 0.5)
        @test GI.getpoint(smoothed, 5) == (1.0, 1.0)
        @test GI.getpoint(smoothed, 6) == (1.5, 0.5)
        @test GI.getpoint(smoothed, 7) == (1.8125, 0.1875)
        @test GI.getpoint(smoothed, 8) == (1.9375, 0.0625)
        @test GI.getpoint(smoothed, 9) == (2.0, 0.0)
    end

    @testset "Polygon" begin
        poly = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        smoothed = GO.smooth(poly)
        exterior = GI.getexterior(smoothed)
        @test GI.npoint(exterior) == 9
        @test GI.getpoint(exterior, 1) == (0.25, 0.0)
        @test GI.getpoint(exterior, 2) == (0.75, 0.0)
        @test GI.getpoint(exterior, 3) == (1.0, 0.25)
        @test GI.getpoint(exterior, 4) == (1.0, 0.75)
        @test GI.getpoint(exterior, 5) == (0.75, 1.0)
        @test GI.getpoint(exterior, 6) == (0.25, 1.0)
        @test GI.getpoint(exterior, 7) == (0.0, 0.75)
        @test GI.getpoint(exterior, 8) == (0.0, 0.25)
        @test GI.getpoint(exterior, 9) == (0.25, 0.0)
    end

    @testset "MultiPolygon" begin
        poly1 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        poly2 = GI.Polygon([[(2.0, 2.0, (3.0, 2.0), (3.0, 3.0), (2.0, 3.0), (2.0, 2.0)]])
        mpoly = GI.MultiPolygon([poly1, poly2])
        smoothed = GO.smooth(mpoly)
        @test GI.ngeom(smoothed) == 2
        @test GI.npoint(smoothed) == 18
    end
end
