using Test
using DelaunayTriangulation

import GeometryOps as GO
import GeoInterface as GI
using Random

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
    
    @testset "Voronoi with custom boundary" begin
        rng = Xoshiro(0)
    
        points = [(0.25, 0.25), (0.75, 0.25), (0.75, 0.75), (0.25, 0.75), (0.5, 0.5)]
        boundary1 = (((0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0)), (1, 2, 3, 4, 1))
        boundary2 = ([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)], [1, 2, 3, 4, 1])
        boundary3 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])

        vorn0 = GO.voronoi(points) # clipped to convex hull
        vorn1 = GO.voronoi(points, clip_polygon = boundary1, rng = Xoshiro(0))
        vorn2 = GO.voronoi(points, clip_polygon = boundary2, rng = Xoshiro(0))
        vorn3 = GO.voronoi(points, clip_polygon = boundary3, rng = Xoshiro(0))

        for vorn in (vorn0, vorn1, vorn2, vorn3)
            @test length(vorn) == 5
            # All polygons should be valid
            for poly in vorn
                @test GI.isgeometry(poly)
                @test GI.geomtrait(poly) isa GI.PolygonTrait
            end
        end
    end
    
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

    @testset "Cells follow flattened input order" begin
        rng = Xoshiro(21)
        points = [(rand(rng), rand(rng)) for _ in 1:20]
        boundary = GI.Polygon([[(-1.0, -1.0), (2.0, -1.0), (2.0, 2.0), (-1.0, 2.0), (-1.0, -1.0)]])
        # Include a permutation and mixed geometries to exercise flattening order.
        for ordered_points in (points, reverse(points))
            geoms = [GI.Point(first(ordered_points)), GI.MultiPoint(ordered_points[2:end])]
            for clip_polygon in (nothing, boundary)
                polygons = GO.voronoi(geoms; clip_polygon, rng = Xoshiro(0))
                @test length(polygons) == length(ordered_points)
                # A cell's centroid must be closest to its own generator.
                for (i, poly) in enumerate(polygons)
                    center = GO.centroid(poly)
                    distances = [sum(abs2, center .- p) for p in ordered_points]
                    @test argmin(distances) == i
                end
            end
        end
    end

    @testset "Clipping can omit cells without reordering survivors" begin
        points = [(0.0, 0.0), (10.0, 0.0), (0.0, 10.0), (10.0, 10.0)]
        boundary = GI.Polygon([[(-1.0, 6.0), (11.0, 6.0), (11.0, 11.0), (-1.0, 11.0), (-1.0, 6.0)]])
        polygons = GO.voronoi(points; clip_polygon = boundary, rng = Xoshiro(0))
        @test length(polygons) == 2
        for (poly, i) in zip(polygons, (3, 4))
            center = GO.centroid(poly)
            @test argmin([sum(abs2, center .- p) for p in points]) == i
        end
    end

    @testset "Explicit RNG controls triangulation and clipping" begin
        points = [(0.1, 0.2), (0.8, 0.1), (0.9, 0.9), (0.2, 0.8), (0.4, 0.5)]
        boundary = GI.Polygon([[(-1.0, -1.0), (2.0, -1.0), (2.0, 2.0), (-1.0, 2.0), (-1.0, -1.0)]])
        for clip_polygon in (nothing, boundary)
            default_rng_before = copy(Random.default_rng())
            first_result = GO.voronoi(points; clip_polygon, rng = Xoshiro(17))
            @test rand(copy(Random.default_rng()), UInt) == rand(default_rng_before, UInt)
            second_result = GO.voronoi(points; clip_polygon, rng = Xoshiro(17))
            @test GI.coordinates.(first_result) == GI.coordinates.(second_result)
        end
    end

    @testset "Clean clipping input" begin
        points = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0))
        order = (1, 2, 3, 4, 1)
        new_points, new_order = GO._clean_voronoi_clip_point_inputs((points, order))
        @test all(points .== new_points)
        @test all(order .== new_order)

        reverse_points = reverse(points)
        new_points, new_order = GO._clean_voronoi_clip_point_inputs((reverse_points, order))
        @test all(points .== new_points)
        @test all(order .== new_order)

        short_points = points[1:end-1]
        short_order = order[1:end-1]
        new_points, new_order = GO._clean_voronoi_clip_point_inputs((short_points, short_order))
        @test all(points .== new_points)
        @test all(order .== new_order)

        shuffled_combos = shuffle(Xoshiro(0), collect(zip(points, order)))
        shuffled_points, shuffled_order = first.(shuffled_combos), last.(shuffled_combos)
        new_points, new_order = GO._clean_voronoi_clip_point_inputs((shuffled_points, shuffled_order))
        @test all(points .== new_points)
        @test all(order .== new_order)
    end
end