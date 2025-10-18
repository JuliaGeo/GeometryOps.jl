using Test
import GeometryOps as GO
import LibGEOS as LG
import Proj

# Basic geometries for testing
point = LG.Point([0.0, 0.0])
line = LG.LineString([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]])
square = LG.Polygon([[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]]])

@testset "Basic perimeter tests" begin
    # Points have zero perimeter
    @test GO.perimeter(point) == 0
    
    # Lines have perimeter equal to their length
    @test GO.perimeter(line) == 2.0  # 1.0 + 1.0
    
    # Square has perimeter of 4 (each side is 1)
    @test GO.perimeter(square) == 4.0
end

@testset "Spherical and geodesic" begin
    highlat_poly = LG.Polygon([[[70., 70.], [70., 80.], [80., 80.], [80., 70.], [70., 70.]]])
    @test GO.perimeter(GO.Planar(), highlat_poly) == 40
    @test GO.perimeter(GO.Planar(), highlat_poly) < GO.perimeter(GO.Spherical(), highlat_poly)
    @test GO.perimeter(GO.Spherical(), highlat_poly) < GO.perimeter(GO.Geodesic(), highlat_poly)
end