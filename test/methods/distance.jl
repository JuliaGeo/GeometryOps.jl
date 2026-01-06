using Test
import ArchGDAL as AG
import GeoInterface as GI
import GeometryOps as GO 
import LibGEOS as LG
using ..TestHelpers

pt1 = LG.Point([0.0, 0.0])
pt2 = LG.Point([0.0, 1.0])
pt3 = LG.Point([2.5, 2.5])
pt4 = LG.Point([3.0, 3.0])
pt5 = LG.Point([5.1, 5.0])
pt6 = LG.Point([3.0, 1.0])
pt7 = LG.Point([0.1, 4.9])
pt8 = LG.Point([2.0, 1.1])
pt9 = LG.Point([3.5, 3.1])
pt10 = LG.Point([10.0, 10.0])
pt11 = LG.Point([2.5, 7.0])

mpt1 = LG.MultiPoint([pt1, pt2, pt3])

l1 = LG.LineString([[0.0, 0.0], [0.0, 5.0], [5.0, 5.0]])

r1 = LG.LinearRing([[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [0.0, 0.0]])
r2 = LG.LinearRing([[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [5.0, 0.0], [0.0, 0.0]])
r3 = LG.LinearRing([[1.0, 1.0], [3.0, 2.0], [4.0, 1.0], [1.0, 1.0]])
r4 = LG.LinearRing([[4.0, 3.0], [3.0, 3.0], [4.0, 4.0], [4.0, 3.0]])
r5 = LG.LinearRing([[0.0, 6.0], [2.5, 8.0], [5.0, 6.0], [0.0, 6.0]])

p1 = LG.Polygon(r2, [r3, r4])
p2 = LG.Polygon(r5)

mp1 = LG.MultiPolygon([p1, p2])

c1 = LG.GeometryCollection([pt1, r1, p1])

@testset_implementations "Where LinearRing exists" [LG, GI] begin
    # Point on linear ring
    @test GO.distance($pt1, $r1) == LG.distance($pt1, $r1)
    @test GO.distance($pt3, $r1) == LG.distance($pt3, $r1)
    # Point outside of linear ring
    @test GO.distance($pt5, $r1) ≈ LG.distance($pt5, $r1)
    # Point inside of hole created by linear ring
    @test GO.distance($pt3, $r1) ≈ LG.distance($pt3, $r1)
    @test GO.distance($pt4, $r1) ≈ LG.distance($pt4, $r1)
end

@testset_implementations "Where GeometryCollection exists" [LG, AG, GI] begin
    @test GO.distance($pt1, c1) == LG.distance($pt1, c1)
end

@testset_implementations "Point and Point" begin
    # Distance from point to same point
    @test GO.distance($pt1, $pt1) == LG.distance($pt1, $pt1)
    # Distance from point to different point
    @test GO.distance($pt1, $pt2) ≈ GO.distance($pt2, $pt1) ≈ LG.distance($pt1, $pt2)
    # Return types
    @test GO.distance($pt1, $pt1) isa Float64
    @test GO.distance($pt1, $pt1, Float32) isa Float32
end

@testset_implementations "Point and Line" begin
    #Point on line vertex
    @test GO.distance($pt1, $l1) == GO.distance($l1, $pt1) == LG.distance($pt1, $l1)
    # Point on line edge
    @test GO.distance($pt2, $l1) == GO.distance($l1, $pt2) == LG.distance($pt2, $l1)
    # Point equidistant from both segments
    @test GO.distance($pt3, $l1) ≈  GO.distance($l1, $pt3) ≈  LG.distance($pt3, $l1)
    # Point closer to one segment than another
    @test GO.distance($pt4, $l1) ≈  GO.distance($l1, $pt4) ≈  LG.distance($pt4, $l1)
    # Return types
    @test GO.distance($pt1, $l1) isa Float64
    @test GO.distance($pt1, $l1, Float32) isa Float32
end

@testset_implementations "Point and Polygon" begin
    # Point on polygon exterior edge
    @test GO.distance($pt1, $p1) == LG.distance($pt1, $p1)
    @test GO.signed_distance($pt1, $p1) == 0
    @test GO.distance($pt2, $p1) == LG.distance($pt2, $p1)
    # Point on polygon hole edge
    @test GO.distance($pt4, $p1) == LG.distance($pt4, $p1)
    @test GO.signed_distance($pt4, $p1) == 0
    @test GO.distance($pt6, $p1) == LG.distance($pt6, $p1)
    # Point inside of polygon
    @test GO.distance($pt3, $p1) == LG.distance($pt3, $p1)
    @test GO.signed_distance($pt3, $p1) ≈
        -(min(LG.distance($pt3, r2), LG.distance($pt3, r3), LG.distance($pt3, r4)))
    @test GO.distance($pt7, $p1) == LG.distance($pt7, $p1)
    @test GO.signed_distance($pt7, $p1) ≈
        -(min(LG.distance($pt7, r2), LG.distance($pt7, r3), LG.distance($pt7, r4)))
    # Point outside of polygon exterior
    @test GO.distance($pt5, $p1) ≈ LG.distance($pt5, $p1)
    @test GO.signed_distance($pt5, $p1) ≈ LG.distance($pt5, $p1)
    # Point inside of polygon hole
    @test GO.distance($pt8, $p1) ≈ LG.distance($pt8, $p1)
    @test GO.signed_distance($pt8, $p1) ≈ LG.distance($pt8, $p1)
    @test GO.distance($pt9, $p1) ≈ LG.distance($pt9, $p1)
    # Return types
    @test GO.distance($pt1, $p1) isa Float64
    @test GO.distance($pt1, $p1, Float32) isa Float32
end

@testset_implementations "Point and MultiPoint" begin
    @test GO.distance($pt4, $mpt1) == LG.distance(pt4, mpt1)
    @test GO.distance($pt4, $mpt1) isa Float64
    @test GO.distance($pt4, $mpt1, Float32) isa Float32
end

@testset_implementations "Point and MultiPolygon" begin
    # Point outside of either polygon
    @test GO.distance($pt5, mp1) ≈ LG.distance($pt5, mp1)
    @test GO.distance($pt10, mp1) ≈ LG.distance($pt10, mp1)
    # Point within one polygon
    @test GO.distance($pt3, mp1) == LG.distance($pt3, mp1)
    @test GO.signed_distance($pt3, mp1) ≈
        -(min(LG.distance($pt3, r2), LG.distance($pt3, r3), LG.distance($pt3, r4), LG.distance($pt3, r5)))
    @test GO.distance($pt11, mp1) == LG.distance($pt11, mp1)
    @test GO.signed_distance($pt11, mp1) ≈
        -(min(LG.distance($pt11, r2), LG.distance($pt11, r3), LG.distance($pt11, r4), LG.distance($pt11, r5)))
end

# Tests for spherical distance using Spherical() manifold
@testset "Spherical distance" begin
    # Earth's mean radius in meters (WGS84)
    EARTH_RADIUS = 6371008.8

    # Helper function to compute expected great-circle distance using Haversine formula
    function haversine_distance(lon1, lat1, lon2, lat2)
        # Convert to radians
        λ1, φ1 = deg2rad(lon1), deg2rad(lat1)
        λ2, φ2 = deg2rad(lon2), deg2rad(lat2)

        Δλ = λ2 - λ1
        Δφ = φ2 - φ1

        a = sin(Δφ/2)^2 + cos(φ1) * cos(φ2) * sin(Δλ/2)^2
        c = 2 * atan(sqrt(a), sqrt(1-a))

        return EARTH_RADIUS * c
    end

    @testset "Point to Point - Spherical" begin
        # Same point should have zero distance
        pt_origin = GI.Point(0.0, 0.0)  # (lon, lat)
        @test GO.distance(GO.Spherical(), pt_origin, pt_origin) ≈ 0.0 atol=1e-10

        # Two points on the equator separated by 90° longitude
        # Distance should be 1/4 of Earth's circumference ≈ 10,007.5 km
        pt_eq1 = GI.Point(0.0, 0.0)   # Prime meridian, equator
        pt_eq2 = GI.Point(90.0, 0.0)  # 90°E, equator
        expected_quarter_circumference = EARTH_RADIUS * π / 2
        @test GO.distance(GO.Spherical(), pt_eq1, pt_eq2) ≈ expected_quarter_circumference rtol=1e-6

        # North Pole to South Pole - half of Earth's circumference
        pt_north = GI.Point(0.0, 90.0)   # North Pole
        pt_south = GI.Point(0.0, -90.0)  # South Pole
        expected_half_circumference = EARTH_RADIUS * π
        @test GO.distance(GO.Spherical(), pt_north, pt_south) ≈ expected_half_circumference rtol=1e-6

        # Distance from equator to North Pole (quarter circumference)
        @test GO.distance(GO.Spherical(), pt_eq1, pt_north) ≈ expected_quarter_circumference rtol=1e-6

        # Antipodal points (opposite sides of Earth)
        pt_antipode = GI.Point(180.0, 0.0)
        expected_half_circumference = EARTH_RADIUS * π
        @test GO.distance(GO.Spherical(), pt_eq1, pt_antipode) ≈ expected_half_circumference rtol=1e-6

        # Real-world known distance: London to New York
        # London: approximately (51.5074°N, -0.1278°W)
        # New York: approximately (40.7128°N, -74.0060°W)
        pt_london = GI.Point(-0.1278, 51.5074)
        pt_nyc = GI.Point(-74.0060, 40.7128)
        expected_distance = haversine_distance(-0.1278, 51.5074, -74.0060, 40.7128)
        @test GO.distance(GO.Spherical(), pt_london, pt_nyc) ≈ expected_distance rtol=1e-4

        # Short distance: London to Paris
        # Paris: approximately (48.8566°N, 2.3522°E)
        pt_paris = GI.Point(2.3522, 48.8566)
        expected_distance_paris = haversine_distance(-0.1278, 51.5074, 2.3522, 48.8566)
        @test GO.distance(GO.Spherical(), pt_london, pt_paris) ≈ expected_distance_paris rtol=1e-4

        # Test commutativity: distance(a, b) == distance(b, a)
        @test GO.distance(GO.Spherical(), pt_london, pt_nyc) ≈ GO.distance(GO.Spherical(), pt_nyc, pt_london)

        # Test type parameter
        @test GO.distance(GO.Spherical(), pt_eq1, pt_eq2, Float64) isa Float64
        @test GO.distance(GO.Spherical(), pt_eq1, pt_eq2, Float32) isa Float32
    end

    @testset "Points along meridians and parallels" begin
        # Points along the same meridian (longitude constant)
        # 1 degree of latitude ≈ 111.32 km at any location
        pt_lat0 = GI.Point(0.0, 0.0)
        pt_lat1 = GI.Point(0.0, 1.0)
        expected_1deg_lat = haversine_distance(0.0, 0.0, 0.0, 1.0)
        @test GO.distance(GO.Spherical(), pt_lat0, pt_lat1) ≈ expected_1deg_lat rtol=1e-4

        # Points along the equator (latitude = 0)
        # 1 degree of longitude at equator ≈ 111.32 km
        pt_lon0 = GI.Point(0.0, 0.0)
        pt_lon1 = GI.Point(1.0, 0.0)
        expected_1deg_lon_eq = haversine_distance(0.0, 0.0, 1.0, 0.0)
        @test GO.distance(GO.Spherical(), pt_lon0, pt_lon1) ≈ expected_1deg_lon_eq rtol=1e-4

        # Distance at higher latitude (longitude circles are smaller)
        # At 60°N, 1 degree of longitude ≈ 55.66 km (half of equatorial)
        pt_60n_lon0 = GI.Point(0.0, 60.0)
        pt_60n_lon1 = GI.Point(1.0, 60.0)
        expected_1deg_lon_60 = haversine_distance(0.0, 60.0, 1.0, 60.0)
        @test GO.distance(GO.Spherical(), pt_60n_lon0, pt_60n_lon1) ≈ expected_1deg_lon_60 rtol=1e-4
    end
end
