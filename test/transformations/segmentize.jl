using Test

import GeometryOps as GO
import GeoInterface as GI

@testset "Segmentation on multiple geometry levels" begin
    ls = GI.LineString([(0, 0), (1, 1), (2, 2), (3, 3)])
    lr = GI.LinearRing([(0, 0), (1, 1), (1, 0), (0, 1), (0, 0)])
    p = GI.Polygon([lr])
    mp = GI.MultiPolygon([p, p, p])
    mls = GI.MultiLineString([ls, ls, ls])

    @test_all_implementations (ls, lr, p, mp, mls) begin

        @testset "LinearSegments" begin
            @test GO.segmentize(ls; max_distance = 0.5) isa GI.LineString
            if GI.trait(lr) isa GI.LinearRingTrait
                @test GO.segmentize(lr; max_distance = 0.5) isa GI.LinearRing
            end
            # Test that linear rings are closed after segmentization
            segmentized_ring = GO.segmentize(lr; max_distance = 0.5)
            @test GI.getpoint(segmentized_ring, 1) == GI.getpoint(segmentized_ring, GI.npoint(segmentized_ring))

            @test GO.segmentize(p; max_distance = 0.5) isa GI.Polygon
            @test GO.segmentize(mp; max_distance = 0.5) isa GI.MultiPolygon
            @test GI.ngeom(GO.segmentize(mp; max_distance = 0.5)) == 3

            # Now test multilinestrings
            @test GO.segmentize(mls; max_distance = 0.5) isa GI.MultiLineString
            @test GI.ngeom(GO.segmentize(mls; max_distance = 0.5)) == 3
        end

        @testset "GeodesicSegments" begin

            @test GO.segmentize(GO.GeodesicSegments(; max_distance = 0.5*900), ls) isa GI.LineString
            if GI.trait(lr) isa GI.LinearRingTrait
                @test GO.segmentize(GO.GeodesicSegments(; max_distance = 0.5*900), lr) isa GI.LinearRing
            end
            # Test that linear rings are closed after segmentization
            segmentized_ring = GO.segmentize(GO.GeodesicSegments(; max_distance = 0.5*900), lr)
            @test GI.getpoint(segmentized_ring, 1) == GI.getpoint(segmentized_ring, GI.npoint(segmentized_ring))

            p = GI.Polygon([lr])
            mp = GI.MultiPolygon([p, p, p])
            @test GO.segmentize(GO.GeodesicSegments(; max_distance = 0.5*900), p) isa GI.Polygon
            @test GO.segmentize(GO.GeodesicSegments(; max_distance = 0.5*900), mp) isa GI.MultiPolygon
            @test GI.ngeom(GO.segmentize(GO.GeodesicSegments(; max_distance = 0.5*900), mp)) == 3

            # Now test multilinestrings
            mls = GI.MultiLineString([ls, ls, ls])
            @test GO.segmentize(GO.GeodesicSegments(; max_distance = 0.5*900), mls) isa GI.MultiLineString
            @test GI.ngeom(GO.segmentize(GO.GeodesicSegments(; max_distance = 0.5*900), mls)) == 3
        end
    end
end

@testset "LinearSegments" begin
    lr = GI.LinearRing([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
    @test_all_implementations lr begin
        ct = GO.centroid(lr)
        ar = GO.area(lr)
        for max_distance in exp10.(LinRange(log10(0.01), log10(1), 10))
            segmentized = GO.segmentize(GO.LinearSegments(; max_distance), lr)
            @test all(GO.centroid(segmentized) .≈ ct)
            @test GO.area(segmentized) ≈ ar
        end
    end
end

@testset "GeodesicSegments" begin
    lr = GI.LinearRing([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
    @test_all_implementations lr begin
        for max_distance in exp10.(LinRange(log10(0.01), log10(1), 10)) .* 900
            @test_nowarn segmentized = GO.segmentize(GO.GeodesicSegments(; max_distance), lr)
        end
    end
end
