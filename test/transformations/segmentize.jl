using Test
using Proj
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG
import ArchGDAL as AG
import GeometryBasics as GB
using GeometryOpsTestHelpers

@testset "Segmentation on multiple geometry levels" begin
    ls = GI.LineString([(0, 0), (1, 1), (2, 2), (3, 3)])
    lr = GI.LinearRing([(0, 0), (1, 1), (1, 0), (0, 1), (0, 0)])
    p = GI.Polygon([lr])
    mp = GI.MultiPolygon([p, p, p])
    mls = GI.MultiLineString([ls, ls, ls])

    @testset_implementations "LinearSegments" begin
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

    @testset_implementations "SphericalSegments" begin
        @test GO.segmentize(GO.Spherical(), $ls; max_distance = 0.5*900) isa GI.LineString
        if GI.trait($lr) isa GI.LinearRingTrait
            @test GO.segmentize(GO.Spherical(), $lr; max_distance = 0.5*900) isa GI.LinearRing
        end
        # Test that linear rings are closed after segmentization
        segmentized_ring = GO.segmentize(GO.Spherical(), $lr; max_distance = 0.5*900)
        @test GI.getpoint(segmentized_ring, 1) == GI.getpoint(segmentized_ring, GI.npoint(segmentized_ring))
        @test GO.segmentize(GO.Spherical(), $p; max_distance = 0.5*900) isa GI.Polygon
        @test GO.segmentize(GO.Spherical(), $mp; max_distance = 0.5*900) isa GI.MultiPolygon
        @test GI.ngeom(GO.segmentize(GO.Spherical(), $mp; max_distance = 0.5*900)) == 3

        # Now test multilinestrings
        @test GO.segmentize(GO.Spherical(), $mls; max_distance = 0.5*900) isa GI.MultiLineString
        @test GI.ngeom(GO.segmentize(GO.Spherical(), $mls; max_distance = 0.5*900)) == 3
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

lr = GI.LinearRing([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
@testset_implementations "Geodesic" begin
    for max_distance in exp10.(LinRange(log10(0.01), log10(1), 10)) .* 900
        @test_nowarn segmentized = GO.segmentize(GO.Geodesic(), $lr; max_distance)
    end
end

@testset "Spherical" begin
    # A segment along the 52°N parallel.  The great circle joining its endpoints bows
    # poleward of the parallel itself, which is the whole difference from `Planar`.
    parallel = GI.LineString([(-50.0, 52.0), (50.0, 52.0)])
    # The endpoints are ~0.982 radians apart, so on the unit sphere a `max_distance` of
    # 0.6 splits the arc into exactly two segments, i.e. one new point at the midpoint.
    unit_sphere = GO.Spherical(; radius = 1)

    @testset "Great circle midpoint" begin
        segmentized = GO.segmentize(unit_sphere, parallel; max_distance = 0.6)
        @test GI.npoint(segmentized) == 3
        midpoint = GI.getpoint(segmentized, 2)
        # The great-circle midpoint of (±λ, φ) sits at atand(tand(φ) / cosd(λ)).
        @test GI.x(midpoint) ≈ 0.0 atol = 1e-9
        @test GI.y(midpoint) ≈ atand(tand(52) / cosd(50)) # ≈ 63.33°N
        # ...which is well poleward of the parallel a `Planar` segmentization follows.
        @test GI.y(midpoint) - 52 > 11
        planar = GO.segmentize(GO.Planar(), parallel; max_distance = 60)
        @test all(==(52.0), GI.y.(GI.getpoint(planar)))
    end

    @testset "Endpoints are preserved exactly" begin
        segmentized = GO.segmentize(unit_sphere, parallel; max_distance = 0.01)
        @test GI.getpoint(segmentized, 1) == (-50.0, 52.0)
        @test GI.getpoint(segmentized, GI.npoint(segmentized)) == (50.0, 52.0)
    end

    @testset "`max_distance` is in units of the radius" begin
        Ω = GO.UnitSpherical.spherical_distance(
            GO.UnitSpherical.UnitSphereFromGeographic()((-50.0, 52.0)),
            GO.UnitSpherical.UnitSphereFromGeographic()((50.0, 52.0)),
        )
        for max_distance in exp10.(LinRange(log10(0.01), log10(1), 10))
            segmentized = GO.segmentize(unit_sphere, parallel; max_distance)
            @test GI.npoint(segmentized) == ceil(Int, Ω / max_distance) + 1
            # The same arc, on a sphere of radius R, needs `max_distance * R` to split
            # into the same number of segments.
            radius = 6371000
            scaled = GO.segmentize(GO.Spherical(; radius), parallel; max_distance = max_distance * radius)
            @test GI.npoint(scaled) == GI.npoint(segmentized)
            @test collect(GI.getpoint(scaled)) == collect(GI.getpoint(segmentized))
        end
    end

    @testset "No segment exceeds `max_distance`" begin
        max_distance = 0.1
        segmentized = GO.segmentize(unit_sphere, parallel; max_distance)
        points = GO.UnitSpherical.UnitSphereFromGeographic().(GI.getpoint(segmentized))
        for (a, b) in zip(points, Iterators.drop(points, 1))
            @test GO.UnitSpherical.spherical_distance(a, b) <= max_distance + 1e-12
        end
    end

    @testset "The three-argument form keeps the manifold" begin
        @test collect(GI.getpoint(GO.segmentize(unit_sphere, parallel, 0.6))) ==
            collect(GI.getpoint(GO.segmentize(unit_sphere, parallel; max_distance = 0.6)))
    end

    @testset "Short segments are left alone" begin
        # Nothing is added when the arc is already shorter than `max_distance`.
        short = GI.LineString([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])
        segmentized = GO.segmentize(unit_sphere, short; max_distance = 1)
        @test collect(GI.getpoint(segmentized)) == [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]
        # Repeated points are also fine (zero-length arcs).
        degenerate = GI.LineString([(1.0, 2.0), (1.0, 2.0)])
        @test collect(GI.getpoint(GO.segmentize(unit_sphere, degenerate; max_distance = 0.001))) ==
            [(1.0, 2.0), (1.0, 2.0)]
    end
end
