using Test
using CoordinateTransformations
import GeoInterface as GI
import GeometryOps as GO
using ..TestHelpers

geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]), 
                   GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])])
rotation_geom = GI.Polygon([GI.LinearRing([(1.0, 1.0), (3.0, 1.0), (3.0, 2.0), (1.0, 2.0), (1.0, 1.0)])])

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

@testset_implementations "transform polygon rotation around the origin" begin
    rotation = LinearMap([0.0 -1.0; 1.0 0.0])
    rotated = GO.transform(rotation, $rotation_geom)
    expected_points = [
        (-1.0, 1.0), (-1.0, 3.0), (-2.0, 3.0), (-2.0, 1.0), (-1.0, 1.0),
    ]
    rotated_points = map(collect(GO.flatten(GI.PointTrait, rotated))) do p
        (GI.x(p), GI.y(p))
    end
    @test rotated_points == expected_points
end

@testset_implementations "transform polygon rotation around centroid preserves centroid and area" begin
    rotation = LinearMap([0.0 -1.0; 1.0 0.0])
    center = GO.centroid($rotation_geom)
    rotated = GO.transform(
        Translation(center...) ∘ rotation ∘ Translation((-).(center)...),
        $rotation_geom,
    )
    @test all(GO.centroid(rotated) .≈ center)
    @test GO.area(rotated) ≈ GO.area($rotation_geom)
end
