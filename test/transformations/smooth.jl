using Test
import GeometryOps as GO
import GeoInterface as GI

@testset "Chaikin" begin
    @testset "LineString" begin
        line = GI.LineString([(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)])
        expected = [(0.0, 0.0), (0.25, 0.25), (0.75, 0.75), (1.25, 0.75), (1.75, 0.25), (2.0, 0.0)]
        smoothed = GO.smooth(line)
        @test GI.npoint(smoothed) == length(expected)
        for i in 1:length(expected)
            @testset let i = i
                @test GI.getpoint(smoothed, i) == expected[i]
            end
        end

        smoothed = GO.smooth(line; iterations=2)
        @test GI.npoint(smoothed) == 12
        expected = [(0.0, 0.0), (0.0625, 0.0625), (0.1875, 0.1875), (0.375, 0.375), (0.625, 0.625), (0.875, 0.75), (1.125, 0.75), (1.375, 0.625), (1.625, 0.375), (1.8125, 0.1875), (1.9375, 0.0625), (2.0, 0.0)]
        for i in 1:length(expected)
            @testset let i = i
                @test GI.getpoint(smoothed, i) == expected[i]
            end
        end
    end

    @testset "Polygon" begin
        poly = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        smoothed = GO.smooth(poly)
        exterior = GI.getexterior(smoothed)
        expected = [(0.25, 0.0), (0.75, 0.0), (1.0, 0.25), (1.0, 0.75), (0.75, 1.0), (0.25, 1.0), (0.0, 0.75), (0.0, 0.25), (0.25, 0.0)]
        @test GI.npoint(exterior) == length(expected)
        for i in 1:length(expected)
            @testset let i = i
                @test GI.getpoint(exterior, i) == expected[i]
            end
        end
    end

    @testset "MultiPolygon" begin
        poly1 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        poly2 = GI.Polygon([[(2.0, 2.0), (3.0, 2.0), (3.0, 3.0), (2.0, 3.0), (2.0, 2.0)]])
        mpoly = GI.MultiPolygon([poly1, poly2])
        smoothed = GO.smooth(mpoly)
        @test GI.ngeom(smoothed) == 2
        @test GI.npoint(smoothed) == 18
    end

    @testset "Spherical" begin
        line = GO.transform(GO.UnitSphereFromGeographic(), GI.LineString([(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)]))
        expected = [line.geom[1], GO.slerp(line.geom[1], line.geom[2], 0.25), GO.slerp(line.geom[1], line.geom[2], 0.75), GO.slerp(line.geom[2], line.geom[3], 0.25), GO.slerp(line.geom[2], line.geom[3], 0.75), line.geom[end]]
        smoothed = GO.smooth(GO.Spherical(), line)
        @test GI.npoint(smoothed) == length(expected)
        for i in 1:length(expected)
            @testset let i = i
                @test GI.getpoint(smoothed, i) == expected[i]
            end
        end
    end
end
