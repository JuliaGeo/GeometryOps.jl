using Test
import GeometryOps as GO
import GeoInterface as GI

@testset "Arclength Functionality" begin
    
    # Test simple horizontal line
    @testset "Simple horizontal line" begin
        line = GI.LineString([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)])
        
        # Test arclength_to_point
        @test GO.arclength_to_point(line, (0.0, 0.0)) ≈ 0.0
        @test GO.arclength_to_point(line, (1.0, 0.0)) ≈ 1.0
        @test GO.arclength_to_point(line, (2.0, 0.0)) ≈ 2.0
        @test GO.arclength_to_point(line, (0.5, 0.0)) ≈ 0.5
        @test GO.arclength_to_point(line, (1.5, 0.0)) ≈ 1.5
        
        # Test point_at_arclength
        point = GO.point_at_arclength(line, 0.0)
        @test point[1] ≈ 0.0 && point[2] ≈ 0.0
        point = GO.point_at_arclength(line, 1.0)
        @test point[1] ≈ 1.0 && point[2] ≈ 0.0
        point = GO.point_at_arclength(line, 2.0)
        @test point[1] ≈ 2.0 && point[2] ≈ 0.0
        point = GO.point_at_arclength(line, 0.5)
        @test point[1] ≈ 0.5 && point[2] ≈ 0.0
        point = GO.point_at_arclength(line, 1.5)
        @test point[1] ≈ 1.5 && point[2] ≈ 0.0
    end
    
    # Test L-shaped line
    @testset "L-shaped line" begin
        line = GI.LineString([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])
        
        # Test arclength_to_point
        @test GO.arclength_to_point(line, (0.0, 0.0)) ≈ 0.0
        @test GO.arclength_to_point(line, (1.0, 0.0)) ≈ 1.0
        @test GO.arclength_to_point(line, (1.0, 1.0)) ≈ 2.0
        @test GO.arclength_to_point(line, (0.5, 0.0)) ≈ 0.5
        @test GO.arclength_to_point(line, (1.0, 0.5)) ≈ 1.5
        
        # Test point_at_arclength
        point = GO.point_at_arclength(line, 0.0)
        @test point[1] ≈ 0.0 && point[2] ≈ 0.0
        point = GO.point_at_arclength(line, 1.0)
        @test point[1] ≈ 1.0 && point[2] ≈ 0.0
        point = GO.point_at_arclength(line, 2.0)
        @test point[1] ≈ 1.0 && point[2] ≈ 1.0
        point = GO.point_at_arclength(line, 0.5)
        @test point[1] ≈ 0.5 && point[2] ≈ 0.0
        point = GO.point_at_arclength(line, 1.5)
        @test point[1] ≈ 1.0 && point[2] ≈ 0.5
    end
    
    # Test with explicit Planar manifold
    @testset "Explicit Planar manifold" begin
        line = GI.LineString([(0.0, 0.0), (3.0, 4.0)])  # 3-4-5 triangle
        
        @test GO.arclength_to_point(GO.Planar(), line, (0.0, 0.0)) ≈ 0.0
        @test GO.arclength_to_point(GO.Planar(), line, (3.0, 4.0)) ≈ 5.0
        @test GO.arclength_to_point(GO.Planar(), line, (1.5, 2.0)) ≈ 2.5
        
        point = GO.point_at_arclength(GO.Planar(), line, 0.0)
        @test point[1] ≈ 0.0 && point[2] ≈ 0.0
        point = GO.point_at_arclength(GO.Planar(), line, 5.0)
        @test point[1] ≈ 3.0 && point[2] ≈ 4.0
        point = GO.point_at_arclength(GO.Planar(), line, 2.5)
        @test point[1] ≈ 1.5 && point[2] ≈ 2.0
    end
    
    # Test edge cases
    @testset "Edge cases" begin
        # Two point line (minimum valid linestring)
        two_point_line = GI.LineString([(0.0, 0.0), (1.0, 0.0)])
        @test GO.arclength_to_point(two_point_line, (0.0, 0.0)) ≈ 0.0
        @test GO.arclength_to_point(two_point_line, (1.0, 0.0)) ≈ 1.0
        
        # Distance beyond line length
        line = GI.LineString([(0.0, 0.0), (1.0, 0.0)])
        point = GO.point_at_arclength(line, 10.0)
        @test point[1] ≈ 1.0 && point[2] ≈ 0.0
        
        # Negative distance
        point = GO.point_at_arclength(line, -1.0)
        @test point[1] ≈ 0.0 && point[2] ≈ 0.0
    end
    
    # Test with LinearRing
    @testset "LinearRing" begin
        ring = GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)])
        
        # Test some points on the ring
        @test GO.arclength_to_point(ring, (0.0, 0.0)) ≈ 0.0
        @test GO.arclength_to_point(ring, (1.0, 0.0)) ≈ 1.0
        @test GO.arclength_to_point(ring, (1.0, 1.0)) ≈ 2.0
        @test GO.arclength_to_point(ring, (0.0, 1.0)) ≈ 3.0
        
        # Test point_at_arclength
        point = GO.point_at_arclength(ring, 0.0)
        @test point[1] ≈ 0.0 && point[2] ≈ 0.0
        point = GO.point_at_arclength(ring, 1.0)
        @test point[1] ≈ 1.0 && point[2] ≈ 0.0
        point = GO.point_at_arclength(ring, 2.5)
        @test point[1] ≈ 0.5 && point[2] ≈ 1.0
    end
    
end