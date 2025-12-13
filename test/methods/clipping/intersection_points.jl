using Test
import GeoInterface as GI
import GeometryOps as GO
import LibGEOS as LG
using ..TestHelpers

l1 = GI.LineString([(90000.0, 1000.0), (90000.0, 22500.0), (95000.0, 22500.0), (95000.0, 1000.0), (90000.0, 1000.0)])
l2 = GI.LineString([(90000.0, 7500.0), (107500.0, 27500.0), (112500.0, 27500.0), (95000.0, 7500.0), (90000.0, 7500.0)])
l3 = GI.LineString([(90000.0, 90000.0), (90000.0, 105000.0), (105000.0, 105000.0), (105000.0, 90000.0), (90000.0, 90000.0)])
l4 = GI.LineString([(-98000.0, 90000.0), (-98000.0, 105000.0), (98000.0, 105000.0), (98000.0, 90000.0), (-98000.0, 90000.0)])
l5 = GI.LineString([(19999.999, 25000.0), (19999.999, 29000.0), (39999.998999999996, 29000.0), (39999.998999999996, 25000.0), (19999.999, 25000.0)])
l6 = GI.LineString([(0.0, 25000.0), (0.0, 29000.0), (20000.0, 29000.0), (20000.0, 25000.0), (0.0, 25000.0)])

p1, p2 = GI.Polygon([l1]), GI.Polygon([l2])

@testset_implementations begin
    # Three intersection points
    LG_l1_l2_mp = GI.MultiPoint(collect(GI.getpoint(LG.intersection($l1, $l2))))
    @test GO.equals(GI.MultiPoint(GO.intersection_points($l1, $l2)), LG_l1_l2_mp)

    # Four intersection points with large intersection
    LG_l3_l4_mp = GI.MultiPoint(collect(GI.getpoint(LG.intersection($l3, $l4))))
    @test GO.equals(GI.MultiPoint(GO.intersection_points($l3, $l4)), LG_l3_l4_mp)

    # Four intersection points with very small intersection
    LG_l5_l6_mp = GI.MultiPoint(collect(GI.getpoint(LG.intersection($l5, $l6))))
    @test GO.equals(GI.MultiPoint(GO.intersection_points($l5, $l6)), LG_l5_l6_mp)

    # Test that intersection points between lines and polygons is equivalent
    @test GO.equals(GI.MultiPoint(GO.intersection_points($p1, $p2)), GI.MultiPoint(GO.intersection_points($l1, $l2)))

    # No intersection points between polygon and line
    @test isempty(GO.intersection_points($p1, $l6))
end

@testset "Spherical intersection points" begin
    # Two lines crossing on the sphere: equator meets prime meridian
    # Line along equator (from -45° to 45° longitude at 0° latitude)
    line1 = GI.LineString([(-45.0, 0.0), (45.0, 0.0)])
    # Line along prime meridian (from -45° to 45° latitude at 0° longitude)
    line2 = GI.LineString([(0.0, -45.0), (0.0, 45.0)])

    pts = GO.intersection_points(GO.Spherical(), line1, line2)
    @test length(pts) == 1
    @test isapprox(pts[1][1], 0.0, atol=1e-6)  # longitude
    @test isapprox(pts[1][2], 0.0, atol=1e-6)  # latitude

    # Two disjoint lines on the sphere (parallel lines on different latitudes)
    line3 = GI.LineString([(-45.0, 10.0), (45.0, 10.0)])
    line4 = GI.LineString([(-45.0, 20.0), (45.0, 20.0)])

    pts = GO.intersection_points(GO.Spherical(), line3, line4)
    @test isempty(pts)

    # Lines that share an endpoint (hinge)
    line5 = GI.LineString([(0.0, 0.0), (30.0, 0.0)])
    line6 = GI.LineString([(0.0, 0.0), (0.0, 30.0)])

    pts = GO.intersection_points(GO.Spherical(), line5, line6)
    @test length(pts) == 1
    @test isapprox(pts[1][1], 0.0, atol=1e-6)
    @test isapprox(pts[1][2], 0.0, atol=1e-6)

    # Test crossing at a non-origin point
    # Two great circle arcs that cross at approximately (45°, 22.5°)
    line7 = GI.LineString([(0.0, 0.0), (90.0, 45.0)])
    line8 = GI.LineString([(90.0, 0.0), (0.0, 45.0)])

    pts = GO.intersection_points(GO.Spherical(), line7, line8)
    @test length(pts) == 1
    # The intersection should be somewhere between the endpoints
    @test 0.0 < pts[1][1] < 90.0  # longitude
    @test 0.0 < pts[1][2] < 45.0  # latitude
end
