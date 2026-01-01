using Test
import GeometryOps as GO
import GeoInterface as GI

@testset "minimum_bounding_circle" begin
    @testset "PlanarCircle basics" begin
        c = GO.PlanarCircle((1.0, 2.0), 4.0)
        @test c.center == (1.0, 2.0)
        @test c.radius_squared == 4.0
        @test GO.radius(c) == 2.0
    end

    @testset "GeoInterface for PlanarCircle" begin
        c = GO.PlanarCircle((0.0, 0.0), 1.0)
        @test GI.geomtrait(c) == GI.PolygonTrait()
        @test GI.nring(c) == 1

        ring = GI.getexterior(c)
        @test GI.geomtrait(ring) == GI.LinearRingTrait()
        @test GI.npoint(ring) == 101

        # First point should be at angle 0 (1, 0)
        p1 = GI.getpoint(ring, 1)
        @test p1[1] ≈ 1.0 atol=1e-10
        @test p1[2] ≈ 0.0 atol=1e-10

        # Point at 90 degrees (index 26 for 100 segments)
        p26 = GI.getpoint(ring, 26)
        @test p26[1] ≈ 0.0 atol=1e-10
        @test p26[2] ≈ 1.0 atol=1e-10

        # Closing point should equal first point
        p101 = GI.getpoint(ring, 101)
        @test p101[1] ≈ p1[1] atol=1e-10
        @test p101[2] ≈ p1[2] atol=1e-10
    end

    @testset "Welzl algorithm" begin
        @testset "Single point" begin
            c = GO.minimum_bounding_circle([(1.0, 2.0)])
            @test c.center == (1.0, 2.0)
            @test c.radius_squared == 0.0
        end

        @testset "Two points" begin
            c = GO.minimum_bounding_circle([(0.0, 0.0), (2.0, 0.0)])
            @test c.center == (1.0, 0.0)
            @test GO.radius(c) ≈ 1.0
        end

        @testset "Three points (equilateral triangle)" begin
            # Equilateral triangle centered at origin
            pts = [(1.0, 0.0), (-0.5, sqrt(3)/2), (-0.5, -sqrt(3)/2)]
            c = GO.minimum_bounding_circle(pts)
            @test c.center[1] ≈ 0.0 atol=1e-10
            @test c.center[2] ≈ 0.0 atol=1e-10
            @test GO.radius(c) ≈ 1.0 atol=1e-10
        end

        @testset "Four points (square)" begin
            pts = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
            c = GO.minimum_bounding_circle(pts)
            @test c.center[1] ≈ 0.5 atol=1e-10
            @test c.center[2] ≈ 0.5 atol=1e-10
            @test GO.radius(c) ≈ sqrt(2)/2 atol=1e-10
        end

        @testset "Points with interior point" begin
            # Square with center point - center point shouldn't affect result
            pts = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.5, 0.5)]
            c = GO.minimum_bounding_circle(pts)
            @test c.center[1] ≈ 0.5 atol=1e-10
            @test c.center[2] ≈ 0.5 atol=1e-10
            @test GO.radius(c) ≈ sqrt(2)/2 atol=1e-10
        end

        @testset "Collinear points" begin
            pts = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]
            c = GO.minimum_bounding_circle(pts)
            @test c.center[1] ≈ 1.5 atol=1e-10
            @test c.center[2] ≈ 0.0 atol=1e-10
            @test GO.radius(c) ≈ 1.5 atol=1e-10
        end

        @testset "From polygon" begin
            poly = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
            c = GO.minimum_bounding_circle(poly)
            @test c.center[1] ≈ 0.5 atol=1e-10
            @test c.center[2] ≈ 0.5 atol=1e-10
            @test GO.radius(c) ≈ sqrt(2)/2 atol=1e-10
        end

        @testset "Empty input" begin
            c = GO.minimum_bounding_circle(Tuple{Float64, Float64}[])
            @test isnan(c.center[1])
            @test isnan(c.center[2])
            @test isnan(c.radius_squared)
        end
    end

    @testset "All points inside circle" begin
        # Generate random points and verify all are inside
        pts = [(rand(), rand()) for _ in 1:50]
        c = GO.minimum_bounding_circle(pts)
        for p in pts
            dx = p[1] - c.center[1]
            dy = p[2] - c.center[2]
            dist_sq = dx^2 + dy^2
            @test dist_sq <= c.radius_squared * (1 + 1e-10)
        end
    end
end

@testset "Spherical minimum_bounding_circle" begin
    @testset "Single point" begin
        # Geographic point (lon, lat)
        c = GO.minimum_bounding_circle(GO.Welzl(; manifold=GO.Spherical()), [(0.0, 0.0)])
        @test c.radius == 0.0
        @test c.radiuslike == 1.0  # cos(0) = 1
    end

    @testset "Two points" begin
        # Two points on the equator, 90 degrees apart
        pts = [(0.0, 0.0), (90.0, 0.0)]  # (lon, lat)
        c = GO.minimum_bounding_circle(GO.Welzl(; manifold=GO.Spherical()), pts)
        # Center should be at (45, 0), radius should be π/4 radians (45 degrees)
        @test c.radius ≈ π/4 atol=1e-10
    end

    @testset "Three points (equilateral spherical triangle)" begin
        # Three points forming an equilateral triangle on the sphere
        # Points at 120 degree intervals on the equator
        pts = [(0.0, 0.0), (120.0, 0.0), (240.0, 0.0)]
        c = GO.minimum_bounding_circle(GO.Welzl(; manifold=GO.Spherical()), pts)
        # All three points should be equidistant from center
        p1 = GO.UnitSpherical.UnitSphereFromGeographic()(pts[1])
        p2 = GO.UnitSpherical.UnitSphereFromGeographic()(pts[2])
        p3 = GO.UnitSpherical.UnitSphereFromGeographic()(pts[3])
        d1 = GO.UnitSpherical.spherical_distance(p1, c.point)
        d2 = GO.UnitSpherical.spherical_distance(p2, c.point)
        d3 = GO.UnitSpherical.spherical_distance(p3, c.point)
        @test d1 ≈ d2 atol=1e-10
        @test d2 ≈ d3 atol=1e-10
        @test d1 ≈ c.radius atol=1e-10
    end

    @testset "Four points with interior point" begin
        # Three boundary points plus one interior point
        pts = [(0.0, 0.0), (10.0, 0.0), (5.0, 8.66), (5.0, 3.0)]  # last point is interior
        c = GO.minimum_bounding_circle(GO.Welzl(; manifold=GO.Spherical()), pts)
        # All points should be inside or on the circle
        for pt in pts
            p = GO.UnitSpherical.UnitSphereFromGeographic()(pt)
            d = GO.UnitSpherical.spherical_distance(p, c.point)
            @test d <= c.radius * (1 + 1e-10)
        end
    end

    @testset "Empty input" begin
        c = GO.minimum_bounding_circle(GO.Welzl(; manifold=GO.Spherical()), Tuple{Float64, Float64}[])
        @test isnan(c.radius)
        @test isnan(c.radiuslike)
    end

    @testset "All points inside circle" begin
        # Generate random geographic points and verify all are inside
        pts = [(rand() * 20 - 10, rand() * 20 - 10) for _ in 1:20]  # small region
        c = GO.minimum_bounding_circle(GO.Welzl(; manifold=GO.Spherical()), pts)
        for pt in pts
            p = GO.UnitSpherical.UnitSphereFromGeographic()(pt)
            d = GO.UnitSpherical.spherical_distance(p, c.point)
            @test d <= c.radius * (1 + 1e-10)
        end
    end

    @testset "From polygon" begin
        poly = GI.Polygon([[(-5.0, -5.0), (5.0, -5.0), (5.0, 5.0), (-5.0, 5.0), (-5.0, -5.0)]])
        c = GO.minimum_bounding_circle(GO.Welzl(; manifold=GO.Spherical()), poly)
        # Verify corners are on or inside the circle
        corners = [(-5.0, -5.0), (5.0, -5.0), (5.0, 5.0), (-5.0, 5.0)]
        for pt in corners
            p = GO.UnitSpherical.UnitSphereFromGeographic()(pt)
            d = GO.UnitSpherical.spherical_distance(p, c.point)
            @test d <= c.radius * (1 + 1e-10)
        end
    end
end
