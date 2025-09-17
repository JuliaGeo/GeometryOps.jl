using Test
import GeoInterface as GI
import GeometryOps as GO 
using GeometryOps: perimeter, Planar, Spherical, Geodesic
using Proj
import LibGEOS as LG
using ..TestHelpers

pt = LG.Point([0.0, 0.0])
empty_pt = LG.readgeom("POINT EMPTY")
mpt = LG.MultiPoint([[0.0, 0.0], [1.0, 0.0]])
empty_mpt = LG.readgeom("MULTIPOINT EMPTY")
l1 = LG.LineString([[0.0, 0.0], [0.5, 0.5], [1.0, 0.5]])
empty_l = LG.readgeom("LINESTRING EMPTY")
ml1 = LG.MultiLineString([[[0.0, 0.0], [0.5, 0.5], [1.0, 0.5]], [[0.0, 0.0], [0.0, 0.1]]])
empty_ml = LG.readgeom("MULTILINESTRING EMPTY")
r1 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [1.0, 2.0], [0.0, 0.0]])
empty_r = LG.readgeom("LINEARRING EMPTY")
p1 = LG.Polygon([
    [[10.0, 0.0], [30.0, 0.0], [30.0, 20.0], [10.0, 20.0], [10.0, 0.0]],
])
p2 = LG.Polygon([
    [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
    [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]],
])
p3 = LG.Polygon([
    [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
    [[15.0, 1.0], [25.0, 1.0], [25.0, 11.0], [15.0, 11.0], [15.0, 1.0]],
])
p4 = LG.Polygon([
    [
        [0.0, 5.0], [2.0, 2.0], [5.0, 2.0], [2.0, -2.0], [5.0, -5.0],
        [0.0, -2.0], [-5.0, -5.0], [-2.0, -2.0], [-5.0, 2.0], [-2.0, 2.0],
        [0.0, 5.0],
    ],
])
empty_p = LG.readgeom("POLYGON EMPTY")
mp1 = LG.MultiPolygon([p2, p4])
empty_mp = LG.readgeom("MULTIPOLYGON EMPTY")
c = LG.GeometryCollection([p1, p2, r1])
c_with_epty_l = LG.GeometryCollection([p1, p2, r1, empty_l])
empty_c = LG.readgeom("GEOMETRYCOLLECTION EMPTY")

# Simple square for testing
square = LG.Polygon([[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]]])

@testset_implementations "That handle empty geoms" begin 
    @test GO.perimeter($empty_pt) == 0
    @test GO.perimeter($empty_mpt) == 0
    @test GO.perimeter($empty_l) == 0
    @test GO.perimeter($empty_ml) == 0
    @test GO.perimeter($empty_r) == 0
    @test GO.perimeter($empty_p) == 0
    @test GO.perimeter($empty_mp) == 0
    @test GO.perimeter(c_with_epty_l) ≈ GO.perimeter(c)
    @test GO.perimeter(empty_c) == 0
end

@testset "With GeometryCollection" begin 
    # Geometry collection summed perimeter
    @test GO.perimeter(c) ≈ GO.perimeter(p1) + GO.perimeter(p2) + GO.perimeter(r1)
    @test GO.perimeter(c, Float32) isa Float32
end  

@testset_implementations "all" begin 
    # Points have zero perimeter
    @test GO.perimeter($pt) == 0
    @test GO.perimeter($pt) isa Float64
    @test GO.perimeter($pt, Float32) isa Float32
    @test GO.perimeter($mpt) == 0
    
    # Lines have non-zero perimeter (length)
    @test GO.perimeter($l1) > 0
    @test GO.perimeter($l1) ≈ sqrt(0.5^2 + 0.5^2) + sqrt(0.5^2 + 0.0^2)  # Distance between points
    @test GO.perimeter($ml1) > GO.perimeter($l1)
    
    # Rings have non-zero perimeter
    @test GO.perimeter($r1) > 0
    @test GO.perimeter($r1) ≈ 1.0 + 2.0 + sqrt(1.0^2 + 2.0^2)  # 3 sides of triangle
    
    # Simple square perimeter test
    @test GO.perimeter($square) ≈ 4.0  # Square with side length 1
    @test GO.perimeter($square, Float32) isa Float32
    
    # Polygons have non-zero perimeter
    @test GO.perimeter($p1) > 0
    @test GO.perimeter($p1) ≈ 2 * 20.0 + 2 * 20.0  # Rectangle perimeter: 2*(width + height)
    
    # Polygon with hole
    outer_perimeter = 2 * 20.0 + 2 * 20.0  # Same as p1
    hole_perimeter = 2 * 10.0 + 2 * 10.0   # Hole perimeter
    @test GO.perimeter($p2) ≈ outer_perimeter + hole_perimeter
    
    # Multipolygon calculations work
    @test GO.perimeter($mp1) ≈ GO.perimeter($p2) + GO.perimeter($p4)
    @test GO.perimeter($mp1, Float32) isa Float32
end

@testset "Geodesic perimeter tests" begin
    # Test geodesic perimeter on simple geometries
    geodesic = GO.Geodesic()
    
    # Simple square in geographic coordinates (degrees)
    geo_square = LG.Polygon([[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]]])
    
    # Geodesic perimeter should be different from planar
    planar_perim = GO.perimeter(GO.Planar(), geo_square)
    geodesic_perim = GO.perimeter(geodesic, geo_square)
    @test geodesic_perim != planar_perim
    @test geodesic_perim > 0
    
    # Test type conversion
    @test GO.perimeter(geodesic, geo_square, Float32) isa Float32
    
    # Test with a larger polygon that would show more difference
    large_square = LG.Polygon([[[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0], [0.0, 0.0]]])
    large_planar_perim = GO.perimeter(GO.Planar(), large_square)
    large_geodesic_perim = GO.perimeter(geodesic, large_square)
    @test large_geodesic_perim != large_planar_perim
    @test large_geodesic_perim > 0
    
    # The geodesic perimeter should be larger than planar for this case
    # (due to Earth's curvature over longer distances)
    @test large_geodesic_perim > large_planar_perim
end

@testset "Spherical perimeter tests" begin
    # Test spherical perimeter calculations
    spherical = GO.Spherical(radius = 6371000.0)  # Earth radius in meters
    
    # Simple square in geographic coordinates
    geo_square = LG.Polygon([[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]]])
    
    # Spherical perimeter should be different from planar
    planar_perim = GO.perimeter(GO.Planar(), geo_square)
    spherical_perim = GO.perimeter(spherical, geo_square)
    @test spherical_perim != planar_perim
    @test spherical_perim > 0
    
    # Test type conversion
    @test GO.perimeter(spherical, geo_square, Float32) isa Float32
end
