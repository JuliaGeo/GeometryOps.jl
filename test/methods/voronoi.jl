using Test
using DelaunayTriangulation

import GeometryOps as GO
import GeometryInterface as GI

@testset "Voronoi" begin
    @testset "Basic voronoi tessellation" begin
        # Test with simple points
        points = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]
        polygons = GO.voronoi(points)
        
        # Should return 3 polygons, one for each point
        @test length(polygons) == 3
        
        # Each should be a valid polygon
        for poly in polygons
            @test GI.isgeometry(poly)
            @test GI.geomtrait(poly) isa GI.PolygonTrait
        end
    end
    
    @testset "Voronoi with various input types" begin
        # Test with GeoInterface points
        points = [GI.Point(0.0, 0.0), GI.Point(1.0, 0.0), GI.Point(0.5, 1.0)]
        polygons = GO.voronoi(points)
        @test length(polygons) == 3
        
        # Test with mixed geometry collection
        geoms = [
            GI.Point(0.0, 0.0),
            GI.LineString([(1.0, 0.0), (1.5, 0.5)]),  # Will extract endpoints
            GI.Point(0.5, 1.0)
        ]
        polygons = GO.voronoi(geoms)
        @test length(polygons) == 4  # 1 + 2 + 1 points
    end
    
    # @testset "Voronoi with clipping" begin
    #     # Test clipped vs unclipped
    #     points = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0), (0.5, 0.5)]
        
    #     # Clipped (default)
    #     clipped_polygons = GO.voronoi(points; clip = true)
    #     @test length(clipped_polygons) == 4
        
    #     # Unclipped  
    #     unclipped_polygons = GO.voronoi(points; clip = false)
    #     @test length(unclipped_polygons) == 4
        
    #     # The polygons should be different (unclipped extends to infinity conceptually)
    #     # But both should have valid polygons
    #     for poly in clipped_polygons
    #         @test GI.isgeometry(poly)
    #     end
    #     for poly in unclipped_polygons
    #         @test GI.isgeometry(poly)
    #     end
    # end
    
    # @testset "Voronoi with custom boundary" begin
    #     # Create points inside a square
    #     points = [(0.25, 0.25), (0.75, 0.25), (0.75, 0.75), (0.25, 0.75), (0.5, 0.5)]
        
    #     # Create a square boundary
    #     boundary = GI.Polygon([
    #         GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)])
    #     ])
        
    #     polygons = GO.voronoi(points, boundary)
    #     @test length(polygons) == 5
        
    #     # All polygons should be valid
    #     for poly in polygons
    #         @test GI.isgeometry(poly)
    #         @test GI.geomtrait(poly) isa GI.PolygonTrait
    #     end
    # end
    
    @testset "Grid of points" begin
        # Create a regular grid
        xs = 0.0:0.5:2.0
        ys = 0.0:0.5:2.0
        points = [(x, y) for x in xs for y in ys]
        
        polygons = GO.voronoi(points)
        @test length(polygons) == length(points)
        
        # Each polygon should be valid
        for poly in polygons
            @test GI.isgeometry(poly)
            
            # Get the exterior ring to check it's closed
            ring = GI.getexterior(poly)
            coords = GI.coordinates(ring)
            @test first(coords) ≈ last(coords)  # Ring should be closed
        end
    end
    
    @testset "Error handling" begin
        # Too few points
        @test_throws ArgumentError GO.voronoi([(0.0, 0.0)])
        @test_throws ArgumentError GO.voronoi([(0.0, 0.0), (1.0, 0.0)])
        
        # Empty input
        @test_throws ArgumentError GO.voronoi([])
    end
    
    @testset "Random points stress test" begin
        # Test with more points
        n = 50
        points = [(rand(), rand()) for _ in 1:n]
        
        polygons = GO.voronoi(points)
        @test length(polygons) == n
        
        # All should be valid polygons
        for poly in polygons
            @test GI.isgeometry(poly)
            @test GI.geomtrait(poly) isa GI.PolygonTrait
        end
    end
end