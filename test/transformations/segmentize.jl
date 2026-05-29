using Test
using Proj
import GeometryOps as GO
import GeoInterface as GI
using ..TestHelpers

@testset "Segmentation on multiple geometry levels" begin
    ls = GI.LineString([(0, 0), (1, 1), (2, 2), (3, 3)])
    lr = GI.LinearRing([(0, 0), (1, 1), (1, 0), (0, 1), (0, 0)])
    p = GI.Polygon([lr])
    mp = GI.MultiPolygon([p, p, p])
    mls = GI.MultiLineString([ls, ls, ls])

    @testset_implementations "LinearSegments max_distance" begin
        @test GO.segmentize($ls; max_distance = 0.5) isa GI.LineString
        if GI.trait($lr) isa GI.LinearRingTrait
            @test GO.segmentize($lr; max_distance = 0.5) isa GI.LinearRing
        end
        # Test that linear rings are closed after segmentization
        segmentized_ring = GO.segmentize($lr; max_distance = 0.5)
        @test GI.getpoint(segmentized_ring, 1) == 
            GI.getpoint(segmentized_ring, GI.npoint(segmentized_ring))

        @test GO.segmentize($p; max_distance = 0.5) isa GI.Polygon
        @test GO.segmentize($mp; max_distance = 0.5) isa GI.MultiPolygon
        @test GI.ngeom(GO.segmentize($mp; max_distance = 0.5)) == 3

        # Now test multilinestrings
        @test GO.segmentize($mls; max_distance = 0.5) isa GI.MultiLineString
        @test GI.ngeom(GO.segmentize($mls; max_distance = 0.5)) == 3
    end

    @testset_implementations "LinearSegments min_number" begin
        @test GO.segmentize($ls; min_number = 7) isa GI.LineString
        if GI.trait($lr) isa GI.LinearRingTrait
            @test GO.segmentize($lr; min_number = 7) isa GI.LinearRing
        end
        # Test that linear rings are closed after segmentization
        segmentized_ring = GO.segmentize($lr; min_number = 7)
        @test GI.getpoint(segmentized_ring, 1) ==
            GI.getpoint(segmentized_ring, GI.npoint(segmentized_ring))

        @test GO.segmentize($p; min_number = 7) isa GI.Polygon
        @test GO.segmentize($mp; min_number = 7) isa GI.MultiPolygon
        @test GI.ngeom(GO.segmentize($mp; min_number = 7)) == 3

        # Now test multilinestrings
        @test GO.segmentize($mls; min_number = 7) isa GI.MultiLineString
        @test GI.ngeom(GO.segmentize($mls; min_number = 7)) == 3
    end

    @testset_implementations "GeodesicSegments" begin
        @test GO.segmentize(GO.Geodesic(), $ls; max_distance = 0.5*900) isa GI.LineString
        if GI.trait($lr) isa GI.LinearRingTrait
            @test GO.segmentize(GO.Geodesic(), $lr; max_distance = 0.5*900) isa GI.LinearRing
        end
        # Test that linear rings are closed after segmentization
        segmentized_ring = GO.segmentize(GO.Geodesic(), $lr; max_distance = 0.5*900)
        @test GI.getpoint(segmentized_ring, 1) == GI.getpoint(segmentized_ring, GI.npoint(segmentized_ring))
        @test GO.segmentize(GO.Geodesic(), $p; max_distance = 0.5*900) isa GI.Polygon
        @test GO.segmentize(GO.Geodesic(), $mp; max_distance = 0.5*900) isa GI.MultiPolygon
        @test GI.ngeom(GO.segmentize(GO.Geodesic(), $mp; max_distance = 0.5*900)) == 3

        # Now test multilinestrings
        @test GO.segmentize(GO.Geodesic(), $mls; max_distance = 0.5*900) isa GI.MultiLineString
        @test GI.ngeom(GO.segmentize(GO.Geodesic(), $mls; max_distance = 0.5*900)) == 3
    end

end

lr = GI.LinearRing([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
@testset_implementations "Planar" begin
    ct = GO.centroid($lr)
    ar = GO.area($lr)
    for max_distance in exp10.(LinRange(log10(0.01), log10(1), 10))
        segmentized = GO.segmentize(GO.Planar(), $lr; max_distance)
        @test all(GO.centroid(segmentized) .≈ ct)
        @test GO.area(segmentized) ≈ ar
    end
end

@testset_implementations "Planar min_number" begin
    ct = GO.centroid($lr)
    ar = GO.area($lr)
    for min_number in 5:15
        segmentized = GO.segmentize(GO.Planar(), $lr; min_number)
        @test all(GO.centroid(segmentized) .≈ ct)
        @test GO.area(segmentized) ≈ ar
    end
end

lr = GI.LinearRing([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
@testset_implementations "Geodesic" begin
    for max_distance in exp10.(LinRange(log10(0.01), log10(1), 10)) .* 900
        @test_nowarn segmentized = GO.segmentize(GO.Geodesic(), $lr; max_distance)
    end
end

@testset "Segmentation with min_number" begin
    @testset_implementations "Simple line" begin
        # If there is no line to be the second longest we risk going
        # into an infinate loop.
        ls = GI.LineString([(0, 0), (10, 2)])
        segmentized = GO.segmentize(ls; min_number = 3)
        @test GI.ngeom(segmentized) == 3
        @test GI.getgeom(segmentized, 2) == (5.0, 1.0)
    end

    @testset_implementations "Geom is larger then min_number" begin
        ls = GI.LineString([(0, 0), (5, 5), (10, 10)])
        segmentized = GO.segmentize(ls; min_number = 2)
        @test GI.ngeom(segmentized) == 3
        @test GI.getgeom(segmentized, 2) == (5, 5)
    end

    @testset_implementations "Cuts long segments in n-parts when second longest segment allows" begin
        ls = GI.LineString([(0.0, 0.0), (7.5, 7.5), (10.0, 10.0)])
        segmentized = GO.segmentize(ls; min_number = 5)
        @test GI.ngeom(segmentized) == 5
        @test GI.getgeom(segmentized, 2) == (2.5, 2.5)
        @test GI.getgeom(segmentized, 3) == (5.0, 5.0)
    end

    @testset_implementations "Points as vectors" begin
        # If there is no line to be the second longest we risk going
        # into an infinate loop.
        ls = GI.LineString([[0, 0], [10, 10]])
        segmentized = GO.segmentize(ls; min_number = 3)
        @test GI.ngeom(segmentized) == 3
        @test GI.getgeom(segmentized, 2) == [5.0, 5.0]
    end

    @testset_implementations "Points as Float32" begin
        # If there is no line to be the second longest we risk going
        # into an infinate loop.
        ls = GI.LineString(Vector{Tuple{Float32, Float32}}([(0.0, 0.0), (10.0, 10.0)]))
        segmentized = GO.segmentize(ls; min_number = 3)
        @test GI.ngeom(segmentized) == 3
        @test GI.getgeom(segmentized, 2) isa Tuple{Float32, Float32}
    end
end

@testset_implementations "Wrong number of keyword arguments" begin
    ls = GI.LineString([(0, 0), (1, 1), (2, 2), (3, 3)])
    @test_throws ArgumentError GO.segmentize(ls)
    @test_throws ArgumentError GO.segmentize(ls; max_distance = 0.5, min_number = 5)
end
