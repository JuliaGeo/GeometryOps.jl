using Test
import GeoInterface as GI
import GeometryOps as GO 
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
empty_l = LG.readgeom("LINESTRING EMPTY")
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

@testset_implementations "That handle empty geoms" begin 
    @test GO.area($empty_pt) == LG.area($empty_pt) == 0
    @test GO.area($empty_mpt) == LG.area($empty_mpt) == 0
    @test GO.area($empty_l) == LG.area($empty_l) == 0
    @test GO.area($empty_ml) == LG.area($empty_ml) == 0
    @test GO.area($empty_r) == LG.area($empty_r) == 0
    # Empty polygon
    @test GO.signed_area($empty_p) == 0
    @test GO.area($empty_p) == LG.area($empty_p) == 0
    # Empty multipolygon
    @test GO.area($empty_mp) == LG.area($empty_mp) == 0
    # Empty collection
    @test GO.area(c_with_epty_l) == LG.area(c_with_epty_l)
    @test GO.area(c_with_epty_l, Float32) isa Float32
    @test GO.area(empty_c) == LG.area(empty_c) == 0
end

@testset "With GeometryCollection" begin 
    # Geometry collection summed area
    @test GO.area(c) == LG.area(c)
    @test GO.area(c, Float32) isa Float32
end  

@testset_implementations "all" begin 
    # Points, lines, and rings have zero area
    @test GO.area($pt) == GO.signed_area($pt) == LG.area($pt) == 0
    @test GO.area($pt) isa Float64
    @test GO.signed_area($pt, Float32) isa Float32
    @test GO.signed_area($pt) isa Float64
    @test GO.area($pt, Float32) isa Float32
    @test GO.area($mpt) == GO.signed_area($mpt) == LG.area($mpt) == 0
    @test GO.area($l1) == GO.signed_area($l1) == LG.area($l1) == 0
    @test GO.area($ml1) == GO.signed_area($ml1) == LG.area($ml1) == 0
    @test GO.area($r1) == GO.signed_area($r1) == LG.area($r1) == 0

    # Polygons have non-zero area
    # CCW polygons have positive signed area
    @test GO.area($p1) == GO.signed_area($p1) == LG.area($p1)
    @test GO.signed_area($p1) > 0
    # Float32 calculations
    @test GO.area($p1) isa Float64
    @test GO.area($p1, Float32) isa Float32
    # CW polygons have negative signed area
    a2 = LG.area($p2)
    @test GO.area($p2) == a2
    @test GO.signed_area($p2) == -a2
    # Winding order of holes doesn't affect sign of signed area
    @test GO.signed_area($p3) == -a2
    # Concave polygon correctly calculates area
    a4 = LG.area($p4)
    @test GO.area($p4) == a4
    @test GO.signed_area($p4) == -a4

    # Multipolygon calculations work
    @test GO.area($mp1) == a2 + a4
    @test GO.area($mp1, Float32) isa Float32
end


highlat_poly = LG.Polygon([[[70., 70.], [70., 80.], [80., 80.], [80., 70.], [70., 70.]]])

@testset_implementations "Spherical/geodesic" begin
    @test GO.area(GO.Planar(), $highlat_poly) == 100
    @test GO.area(GO.Planar(), $highlat_poly) < GO.area(GO.Geodesic(), $highlat_poly)
end

@testset "GirardSphericalArea algorithm type" begin
    # Test that the algorithm type exists and is a SingleManifoldAlgorithm for Spherical
    @test GO.GirardSphericalArea <: GO.GeometryOpsCore.SingleManifoldAlgorithm{<:GO.Spherical}

    # Test construction with default manifold
    alg = GO.GirardSphericalArea()
    @test alg isa GO.GirardSphericalArea

    # Test manifold accessor
    @test GO.GeometryOpsCore.manifold(alg) isa GO.Spherical
end

@testset "Spherical triangle area" begin
    using GeometryOps.UnitSpherical: UnitSphericalPoint, UnitSphereFromGeographic

    # Test a known spherical triangle: 1/8 of the unit sphere
    # Vertices at (0,0), (90,0), (0,90) in lon,lat
    # This is an octant of the sphere with area = 4πR²/8 = πR²/2
    p1 = UnitSphereFromGeographic()((0.0, 0.0))   # lon=0, lat=0
    p2 = UnitSphereFromGeographic()((90.0, 0.0))  # lon=90, lat=0
    p3 = UnitSphereFromGeographic()((0.0, 90.0))  # lon=0, lat=90 (north pole)

    # On unit sphere, area should be π/2
    area = GO._girard_spherical_triangle_area(p1, p2, p3)
    @test area ≈ π/2 atol=1e-10

    # Test with reversed winding (should give same absolute area)
    area_rev = GO._girard_spherical_triangle_area(p1, p3, p2)
    @test abs(area_rev) ≈ π/2 atol=1e-10
end

@testset "Spherical ring area" begin
    # Test the octant triangle as a ring
    ring = GI.LinearRing([(0.0, 0.0), (90.0, 0.0), (0.0, 90.0), (0.0, 0.0)])

    # On unit sphere, area should be π/2
    area = GO._girard_spherical_ring_area(Float64, ring)
    @test abs(area) ≈ π/2 atol=1e-10

    # Test a small polygon where spherical ≈ planar (for sanity check)
    small_ring = GI.LinearRing([
        (0.0, 0.0), (0.001, 0.0), (0.001, 0.001), (0.0, 0.001), (0.0, 0.0)
    ])
    small_area = GO._girard_spherical_ring_area(Float64, small_ring)
    # For very small regions, spherical area should be positive for CCW ring
    @test small_area > 0
end