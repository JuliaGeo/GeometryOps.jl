using Test
using Proj
import GeometryOps as GO
import GeoInterface as GI

@testset "Arclength Functionality with Geodesic" begin
    
    # Test with Geodesic manifold (requires Proj)
    @testset "Geodesic manifold" begin
        # Use geographic coordinates (lat/lon)
        line = GI.LineString([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])
        
        # These should now work without errors
        @test_nowarn GO.arclength_to_point(GO.Geodesic(), line, (0.5, 0.0))
        @test_nowarn GO.point_at_arclength(GO.Geodesic(), line, 50000.0)  # 50km
        
        # Test specific values - geodesic distances will be different from planar
        distance_geo = GO.arclength_to_point(GO.Geodesic(), line, (1.0, 0.0))
        @test distance_geo > 0
        @test distance_geo != 1.0  # Should be different from planar distance
        
        point_geo = GO.point_at_arclength(GO.Geodesic(), line, distance_geo)
        @test abs(point_geo[1] - 1.0) < 0.01  # Should be close to (1.0, 0.0)
        @test abs(point_geo[2] - 0.0) < 0.01
    end
    
    # Test comparison between Planar and Geodesic for small distances
    @testset "Planar vs Geodesic comparison" begin
        # Small line segment where differences should be minimal
        line = GI.LineString([(0.0, 0.0), (0.01, 0.0)])  # About 1.1km at equator
        
        distance_planar = GO.arclength_to_point(GO.Planar(), line, (0.005, 0.0))
        distance_geodesic = GO.arclength_to_point(GO.Geodesic(), line, (0.005, 0.0))
        
        # For small distances, they should be fairly similar but not identical
        @test distance_planar ≈ 0.005
        @test distance_geodesic > 0
        @test abs(distance_planar - distance_geodesic) < distance_planar  # Geodesic should be somewhat different
    end
    
    # Test LinearRing with Geodesic
    @testset "Geodesic LinearRing" begin
        ring = GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)])
        
        @test_nowarn GO.arclength_to_point(GO.Geodesic(), ring, (0.5, 0.0))
        @test_nowarn GO.point_at_arclength(GO.Geodesic(), ring, 100000.0)  # 100km
        
        distance = GO.arclength_to_point(GO.Geodesic(), ring, (1.0, 0.0))
        @test distance > 0
        
        point = GO.point_at_arclength(GO.Geodesic(), ring, distance)
        @test abs(point[1] - 1.0) < 0.01
        @test abs(point[2] - 0.0) < 0.01
    end
    
end