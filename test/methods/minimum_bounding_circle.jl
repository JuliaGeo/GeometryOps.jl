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
