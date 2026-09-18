using Test
import GeometryOps as GO
import LibGEOS as LG

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
@testset "CRS-aware automatic perimeter" begin
    import GeoFormatTypes, Proj
    import GeoInterface as GI

    ring = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]
    geographic_square = GI.Polygon([ring]; crs = GeoFormatTypes.EPSG(4326))
    @test GO.perimeter(geographic_square) == GO.perimeter(GO.Geodesic(), geographic_square)
    @test GO.perimeter(geographic_square) ≈ 443770.91724830196 rtol = 1e-10

    @test GO.perimeter(GI.Polygon([ring])) == 4.0
    @test GO.perimeter(GI.Polygon([ring]; crs = GeoFormatTypes.EPSG(2263))) == 4.0
    @test GO.perimeter([geographic_square]) == 4.0

    # A geodesic line is not closed.
    geographic_line = GI.LineString(ring[1:2]; crs = GeoFormatTypes.EPSG(4326))
    @test GO.perimeter(geographic_line) == GO.perimeter(GO.Geodesic(), geographic_line)

    # EPSG:4807 has coordinates in grads on the Clarke 1880 (IGN) ellipsoid.
    clarke = GO.Geodesic(semimajor_axis = 6378249.2, inv_flattening = 293.466021293627)
    grads_square = GI.Polygon([ring]; crs = GeoFormatTypes.EPSG(4807))
    degree_square = GI.Polygon([[(0.9x, 0.9y) for (x, y) in ring]])
    @test GO.perimeter(grads_square) ≈ GO.perimeter(clarke, degree_square) rtol = 1e-10

    geographic_multipolygon = GI.MultiPolygon(fill(GI.Polygon([ring]), 256); crs = GeoFormatTypes.EPSG(4326))
    @test GO.perimeter(geographic_multipolygon; threaded = true) ≈ 256 * GO.perimeter(geographic_square)
end
