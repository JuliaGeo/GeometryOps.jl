using Test
using CoordinateTransformations
import GeoInterface as GI
import GeometryOps as GO
using ..TestHelpers

geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]), 
                   GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])])

@testset_implementations "transform" begin
    translated = GI.Polygon([GI.LinearRing([[4.5, 3.5], [6.5, 5.5], [8.5, 7.5], [4.5, 3.5]]), 
                             GI.LinearRing([[6.5, 5.5], [8.5, 7.5], [9.5, 8.5], [6.5, 5.5]])])
    f = CoordinateTransformations.Translation(3.5, 1.5)
    @test GO.transform(f, $geom) == translated
end

@testset_implementations "transform 2D to 3D" begin
    flat_points_raw = collect(GO.flatten(GI.PointTrait, $geom))
    flat_points_transformed = map(flat_points_raw) do p
        (GI.x(p), GI.y(p), hypot(GI.x(p), GI.y(p)))
    end

    geom_transformed = GO.transform($geom) do p
        (GI.x(p), GI.y(p), hypot(GI.x(p), GI.y(p)))
    end
    @test collect(GO.flatten(GI.PointTrait, geom_transformed)) == flat_points_transformed
    @test GI.is3d(geom_transformed)
    @test !GI.ismeasured(geom_transformed)
end

@testset_implementations "rotation transformations" begin
    using StaticArrays
    
    # Test simple rotation around origin - rotate by 90 degrees
    rotation_90 = @SMatrix [0.0 -1.0; 1.0 0.0]
    square = GI.Polygon([GI.LinearRing([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])])
    rotated_square = GO.transform(p -> rotation_90 * p, square)
    
    # After 90-degree rotation: (0,0) -> (0,0), (1,0) -> (0,1), (1,1) -> (-1,1), (0,1) -> (-1,0)
    expected_points = [(0.0, 0.0), (0.0, 1.0), (-1.0, 1.0), (-1.0, 0.0), (0.0, 0.0)]
    rotated_points = collect(GI.getpoint(GI.getexterior(rotated_square)))
    
    for (actual, expected) in zip(rotated_points, expected_points)
        @test GI.x(actual) ≈ expected[1] atol=1e-10
        @test GI.y(actual) ≈ expected[2] atol=1e-10
    end
    
    # Test rotation around centroid
    center = GO.centroid(square)  # Should be (0.5, 0.5)
    rotated_around_center = GO.transform(square) do p
        rotated = rotation_90 * (p .- center)
        return rotated .+ center
    end
    
    # Verify that centroid remains the same after rotation around center
    new_center = GO.centroid(rotated_around_center)
    @test new_center[1] ≈ center[1] atol=1e-10
    @test new_center[2] ≈ center[2] atol=1e-10
end

@testset_implementations "rotation with CoordinateTransformations" begin
    using StaticArrays
    
    # Test using CoordinateTransformations for rotation around centroid
    square = GI.Polygon([GI.LinearRing([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])])
    center = GO.centroid(square)
    rotation_matrix = @SMatrix [0.0 -1.0; 1.0 0.0]  # 90 degree rotation
    
    # Compose transformations
    rotation_transform = CoordinateTransformations.Translation(center) ∘ 
                        CoordinateTransformations.LinearMap(rotation_matrix) ∘ 
                        CoordinateTransformations.Translation(-center[1], -center[2])
    
    rotated_square = GO.transform(rotation_transform, square)
    
    # Verify centroid preservation
    new_center = GO.centroid(rotated_square)
    @test new_center[1] ≈ center[1] atol=1e-10
    @test new_center[2] ≈ center[2] atol=1e-10
end

@testset_implementations "rotate convenience function" begin
    using StaticArrays
    
    square = GI.Polygon([GI.LinearRing([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])])
    
    # Test rotation around centroid (default behavior)
    center = GO.centroid(square)
    rotated_default = GO.rotate(square, π/2)
    new_center = GO.centroid(rotated_default) 
    
    # Centroid should be preserved
    @test new_center[1] ≈ center[1] atol=1e-10
    @test new_center[2] ≈ center[2] atol=1e-10
    
    # Test rotation around origin
    rotated_origin = GO.rotate(square, π/2; origin = (0, 0))
    origin_points = collect(GI.getpoint(GI.getexterior(rotated_origin)))
    
    # Check specific point transformations: (1,0) -> (0,1), (1,1) -> (-1,1), etc.
    expected_points = [(0.0, 0.0), (0.0, 1.0), (-1.0, 1.0), (-1.0, 0.0), (0.0, 0.0)]
    
    for (actual, expected) in zip(origin_points, expected_points)
        @test GI.x(actual) ≈ expected[1] atol=1e-10
        @test GI.y(actual) ≈ expected[2] atol=1e-10
    end
    
    # Test rotation around custom point
    custom_origin = (0.5, 0.5)
    rotated_custom = GO.rotate(square, π/2; origin = custom_origin)
    custom_center = GO.centroid(rotated_custom)
    
    # When rotating around (0.5, 0.5), the centroid should stay at (0.5, 0.5)
    @test custom_center[1] ≈ custom_origin[1] atol=1e-10
    @test custom_center[2] ≈ custom_origin[2] atol=1e-10
end
