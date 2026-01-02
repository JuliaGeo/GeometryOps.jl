using Test
using LinearAlgebra

using GeometryOps.UnitSpherical

import GeoInterface as GI

@testset "Coordinate transforms" begin
    @testset "UnitSphereFromGeographic" begin
        # Test with GeoInterface Point
        point = GI.Point(45, 45)
        result = UnitSphereFromGeographic()(point)
        @test result isa UnitSphericalPoint{Float64}
        @test length(result) == 3
        @test isapprox(result[1], 0.5, atol=1e-10)
        @test isapprox(result[2], 0.5, atol=1e-10)
        @test isapprox(result[3], 1/√2, atol=1e-10)

        # Test with tuple
        result = UnitSphereFromGeographic()((45, 45))
        @test result isa UnitSphericalPoint{Float64}
        @test length(result) == 3
        @test isapprox(result[1], 0.5, atol=1e-10)
        @test isapprox(result[2], 0.5, atol=1e-10)
        @test isapprox(result[3], 1/√2, atol=1e-10)

        # Test edge cases
        # North pole
        result = UnitSphereFromGeographic()((0, 90))
        @test isapprox(result[1], 0.0, atol=1e-10)
        @test isapprox(result[2], 0.0, atol=1e-10)
        @test isapprox(result[3], 1.0, atol=1e-10)

        # South pole
        result = UnitSphereFromGeographic()((0, -90))
        @test isapprox(result[1], 0.0, atol=1e-10)
        @test isapprox(result[2], 0.0, atol=1e-10)
        @test isapprox(result[3], -1.0, atol=1e-10)

        # Equator
        result = UnitSphereFromGeographic()((0, 0))
        @test isapprox(result[1], 1.0, atol=1e-10)
        @test isapprox(result[2], 0.0, atol=1e-10)
        @test isapprox(result[3], 0.0, atol=1e-10)
    end

    @testset "GeographicFromUnitSphere" begin
        # Test basic conversion
        point = UnitSphericalPoint(0.5, 0.5, 1/√2)
        result = GeographicFromUnitSphere()(point)
        @test result isa Tuple{Float64,Float64}
        @test isapprox(result[1], 45.0, atol=1e-10)  # longitude
        @test isapprox(result[2], 45.0, atol=1e-10)  # latitude

        # Test edge cases
        # North pole
        result = GeographicFromUnitSphere()(UnitSphericalPoint(0.0, 0.0, 1.0))
        @test isapprox(result[1], 0.0, atol=1e-10)  # longitude (undefined at poles, convention is 0)
        @test isapprox(result[2], 90.0, atol=1e-10)  # latitude

        # South pole
        result = GeographicFromUnitSphere()(UnitSphericalPoint(0.0, 0.0, -1.0))
        @test isapprox(result[1], 0.0, atol=1e-10)  # longitude (undefined at poles, convention is 0)
        @test isapprox(result[2], -90.0, atol=1e-10)  # latitude

        # Equator
        result = GeographicFromUnitSphere()(UnitSphericalPoint(1.0, 0.0, 0.0))
        @test isapprox(result[1], 0.0, atol=1e-10)  # longitude
        @test isapprox(result[2], 0.0, atol=1e-10)  # latitude

        # Test with regular vector
        result = GeographicFromUnitSphere()([0.5, 0.5, 1/√2])
        @test result isa Tuple{Float64,Float64}
        @test isapprox(result[1], 45.0, atol=1e-10)
        @test isapprox(result[2], 45.0, atol=1e-10)

        # Test error handling for non-3D vectors
        @test_throws AssertionError GeographicFromUnitSphere()([1.0, 0.0])
        @test_throws AssertionError GeographicFromUnitSphere()([1.0, 0.0, 0.0, 0.0])
    end
end

@testset "Spherical caps" begin
    @testset "Construction" begin
        # Test construction from UnitSphericalPoint
        point = UnitSphericalPoint(1.0, 0.0, 0.0)
        cap = SphericalCap(point, π/4)
        @test cap.point == point
        @test cap.radius == π/4

        # Test construction from geographic point
        geo_point = GI.Point(45, 45)
        cap = SphericalCap(geo_point, π/4)
        @test cap.point isa UnitSphericalPoint{Float64}
        @test cap.radius == π/4

        # Test construction from tuple
        cap = SphericalCap((45, 45), π/4)
        @test cap.point isa UnitSphericalPoint{Float64}
        @test cap.radius == π/4
    end

    @testset "Intersection and containment" begin
        # Create two caps that intersect
        cap1 = SphericalCap(UnitSphericalPoint(1.0, 0.0, 0.0), π/4)
        cap2 = SphericalCap(UnitSphericalPoint(1/√2, 1/√2, 0.0), π/4)
        @test UnitSpherical._intersects(cap1, cap2)
        @test UnitSpherical._intersects(cap2, cap1)
        @test !UnitSpherical._disjoint(cap1, cap2)

        # Create two caps that don't intersect
        cap3 = SphericalCap(UnitSphericalPoint(1.0, 0.0, 0.0), π/8)
        cap4 = SphericalCap(UnitSphericalPoint(0.0, 0.0, 1.0), π/8)
        @test !UnitSpherical._intersects(cap3, cap4)
        @test UnitSpherical._disjoint(cap3, cap4)

        # Test containment
        big_cap = SphericalCap(UnitSphericalPoint(1.0, 0.0, 0.0), π/2)
        small_cap = SphericalCap(UnitSphericalPoint(1/√2, 1/√2, 0.0), π/4)
        @test_broken UnitSpherical._contains(big_cap, small_cap)
        @test !UnitSpherical._contains(small_cap, big_cap)
    end

    @testset "Circumcenter and circumradius" begin
        # Test with an equilateral triangle on the equator
        a = UnitSphericalPoint(1.0, 0.0, 0.0)
        b = UnitSphericalPoint(-0.5, √3/2, 0.0)
        c = UnitSphericalPoint(-0.5, -√3/2, 0.0)
        cap = SphericalCap(a, b, c)
        
        # The circumcenter should be at the north pole
        @test isapprox(cap.point[1], 0.0, atol=1e-10)
        @test isapprox(cap.point[2], 0.0, atol=1e-10)
        @test isapprox(cap.point[3], 1.0, atol=1e-10)
        # The radius should be π/2 (90 degrees)
        @test isapprox(cap.radius, π/2, atol=1e-10)

        # Test with a triangle in the northern hemisphere
        a = UnitSphericalPoint(1.0, 0.0, 0.0)
        b = UnitSphericalPoint(0.0, 1.0, 0.0)
        c = UnitSphericalPoint(0.0, 0.0, 1.0)
        cap = SphericalCap(a, b, c)
        
        # The circumcenter should be at (1/√3, 1/√3, 1/√3)
        @test isapprox(cap.point[1], 1/√3, atol=1e-10)
        @test isapprox(cap.point[2], 1/√3, atol=1e-10)
        @test isapprox(cap.point[3], 1/√3, atol=1e-10)
        # The radius should be the angle between the center and any vertex
        expected_radius = acos(1/√3)
        @test isapprox(cap.radius, expected_radius, atol=1e-10)

        # Test with nearly colinear points (small angle between them)
        # Points near the equator with very small angular separation
        ϵ = 1e-6  # Very small angle in radians
        a = UnitSphericalPoint(1.0, 0.0, 0.0)
        b = UnitSphericalPoint(cos(ϵ), sin(ϵ), 0.0)
        c = UnitSphericalPoint(cos(2ϵ), sin(2ϵ), 0.0)
        cap = SphericalCap(a, b, c)
        
        # The circumcenter should be near the north pole
        @test isapprox(cap.point[1], 0.0, atol=1e-6)
        @test isapprox(cap.point[2], 0.0, atol=1e-6)
        @test isapprox(cap.point[3], 1.0, atol=1e-6)
        # The radius should be approximately π/2
        @test isapprox(cap.radius, π/2, atol=1e-6)

        # # Test with nearly identical points
        # # Points very close to each other in the northern hemisphere
        # ϵ = 1e-8  # Extremely small angle
        # base = UnitSphericalPoint(1/√3, 1/√3, 1/√3)
        # a = base
        # b = UnitSphericalPoint(cos(ϵ), sin(ϵ), 0.0) * (1/√3)  # Small perturbation
        # c = UnitSphericalPoint(cos(2ϵ), sin(2ϵ), 0.0) * (1/√3)  # Another small perturbation
        # cap = SphericalCap(a, b, c)
        
        # # The circumcenter should be very close to the base point
        # @test isapprox(cap.point[1], 1/√3, atol=1e-6)
        # @test isapprox(cap.point[2], 1/√3, atol=1e-6)
        # @test isapprox(cap.point[3], 1/√3, atol=1e-6)
        # # The radius should be very small
        # @test isapprox(cap.radius, ϵ, atol=1e-6)

        # Test with points forming a very thin triangle
        # Points near the north pole with small angular separation
        ϵ = 1e-6
        a = UnitSphericalPoint(0.0, 0.0, 1.0)  # North pole
        b = UnitSphericalPoint(sin(ϵ), 0.0, cos(ϵ))
        c = UnitSphericalPoint(0.0, sin(ϵ), cos(ϵ))
        cap = SphericalCap(a, b, c)
        
        # The circumcenter should be very close to the north pole
        @test isapprox(cap.point[1], 0.0, atol=1e-6)
        @test isapprox(cap.point[2], 0.0, atol=1e-6)
        @test isapprox(cap.point[3], 1.0, atol=1e-6)
        # The radius should be approximately ϵ
        @test isapprox(cap.radius, ϵ, atol=1e-6)
    end

    @testset "Circumcenter winding order independence" begin
        # The circumcenter should be the same regardless of point order
        # This tests the fix for the antipodal circumcenter bug

        using LinearAlgebra: norm

        # Test with a triangle in the first octant
        a = UnitSphericalPoint(1.0, 0.0, 0.0)
        b = UnitSphericalPoint(0.0, 1.0, 0.0)
        c = UnitSphericalPoint(0.0, 0.0, 1.0)

        # All 6 permutations should give the same circumcenter
        centers = [
            UnitSpherical.circumcenter_on_unit_sphere(a, b, c),
            UnitSpherical.circumcenter_on_unit_sphere(a, c, b),
            UnitSpherical.circumcenter_on_unit_sphere(b, a, c),
            UnitSpherical.circumcenter_on_unit_sphere(b, c, a),
            UnitSpherical.circumcenter_on_unit_sphere(c, a, b),
            UnitSpherical.circumcenter_on_unit_sphere(c, b, a),
        ]

        # All centers should be approximately equal
        for i in 2:6
            @test isapprox(centers[1][1], centers[i][1], atol=1e-10)
            @test isapprox(centers[1][2], centers[i][2], atol=1e-10)
            @test isapprox(centers[1][3], centers[i][3], atol=1e-10)
        end

        # The center should be at (1/√3, 1/√3, 1/√3) - the smaller circumcircle
        expected = 1/√3
        @test isapprox(centers[1][1], expected, atol=1e-10)
        @test isapprox(centers[1][2], expected, atol=1e-10)
        @test isapprox(centers[1][3], expected, atol=1e-10)

        # Test with random triangles - all should give radius < π/2
        # (the smaller circumcircle, not the antipodal one)
        using Random
        Random.seed!(42)
        for _ in 1:100
            # Generate 3 random points on the unit sphere
            points = [UnitSphericalPoint(randn(), randn(), randn()) for _ in 1:3]
            points = [p / norm(p) for p in points]  # Normalize

            cap = SphericalCap(points[1], points[2], points[3])

            # The circumcircle should always be the smaller one (radius ≤ π/2)
            # unless points are nearly collinear (on a great circle)
            @test cap.radius <= π/2 + 1e-10

            # All 3 points should be equidistant from the center
            d1 = UnitSpherical.spherical_distance(cap.point, points[1])
            d2 = UnitSpherical.spherical_distance(cap.point, points[2])
            d3 = UnitSpherical.spherical_distance(cap.point, points[3])
            @test isapprox(d1, d2, atol=1e-10)
            @test isapprox(d2, d3, atol=1e-10)
            @test isapprox(d1, cap.radius, atol=1e-10)
        end
    end

    @testset "Merging of SphericalCaps" begin
        function test_merge(p1, p2, r1, r2, pmerged, rmerged)
            r1 = deg2rad(r1)
            r2 = deg2rad(r2)
            cap1 = SphericalCap(p1, r1)
            cap2 = SphericalCap(p2, r2)
            capmerged = UnitSpherical._merge(cap1, cap2)
            @test all(isapprox.(GeographicFromUnitSphere()(capmerged.point), pmerged))
            @test isapprox(capmerged.radius, deg2rad(rmerged))
        end

        test_merge((10.0, 0.0), (30.0, 0.0), 5, 10, (22.5, 0.0), 17.5)
        test_merge((10.0, 0.0), (30.0, 0.0), 15, 10, (17.5, 0.0), 22.5)
        test_merge((10.0, 0.0), (30.0, 0.0), 40, 5, (10.0, 0.0), 40.0)
        test_merge((10.0, 0.0), (30.0, 0.0), 5, 50, (30.0, 0.0), 50.0)
    end
end

@testset "Spherical orientation predicate" begin
    using GeometryOps.UnitSpherical: spherical_orient

    # Points on the equator
    a = UnitSphericalPoint(1.0, 0.0, 0.0)  # (0°, 0°)
    b = UnitSphericalPoint(0.0, 1.0, 0.0)  # (90°, 0°)

    # Point in northern hemisphere - should be "left" of a→b (positive)
    c_north = UnitSphericalPoint(0.0, 0.0, 1.0)  # North pole
    @test spherical_orient(a, b, c_north) == 1

    # Point in southern hemisphere - should be "right" of a→b (negative)
    c_south = UnitSphericalPoint(0.0, 0.0, -1.0)  # South pole
    @test spherical_orient(a, b, c_south) == -1

    # Point on the great circle - should return 0
    c_on = UnitSphericalPoint(-1.0, 0.0, 0.0)  # (180°, 0°)
    @test spherical_orient(a, b, c_on) == 0

    # Test with nearly collinear points (numerical stability)
    ε = 1e-10
    a_near = UnitSphericalPoint(1.0, 0.0, 0.0)
    b_near = UnitSphericalPoint(1.0 - ε, ε, 0.0)
    b_near = b_near / norm(b_near)  # Normalize
    c_near = UnitSphericalPoint(1.0 - 2ε, 2ε, 0.0)
    c_near = c_near / norm(c_near)
    # Should handle near-collinear gracefully (uses robust_cross_product)
    @test spherical_orient(a_near, b_near, c_near) in (-1, 0, 1)
end

@testset "Point on spherical arc" begin
    using GeometryOps.UnitSpherical: point_on_spherical_arc

    # Arc from (0°,0°) to (90°,0°) along equator
    a = UnitSphericalPoint(1.0, 0.0, 0.0)
    b = UnitSphericalPoint(0.0, 1.0, 0.0)

    # Point at (45°,0°) - midpoint of arc
    mid = UnitSphericalPoint(1/√2, 1/√2, 0.0)
    @test point_on_spherical_arc(mid, a, b) == true

    # Endpoints should be on the arc
    @test point_on_spherical_arc(a, a, b) == true
    @test point_on_spherical_arc(b, a, b) == true

    # Point on great circle but outside arc (at 180°,0°)
    outside = UnitSphericalPoint(-1.0, 0.0, 0.0)
    @test point_on_spherical_arc(outside, a, b) == false

    # Point not on great circle (north pole)
    off_circle = UnitSphericalPoint(0.0, 0.0, 1.0)
    @test point_on_spherical_arc(off_circle, a, b) == false

    # Test arc crossing the antimeridian
    # Arc from (170°,0°) to (-170°,0°) - the SHORT way
    a2 = UnitSphereFromGeographic()((170.0, 0.0))
    b2 = UnitSphereFromGeographic()((-170.0, 0.0))
    mid2 = UnitSphereFromGeographic()((180.0, 0.0))
    @test point_on_spherical_arc(mid2, a2, b2) == true

    # Point on the LONG way around (at 0°,0°)
    far = UnitSphericalPoint(1.0, 0.0, 0.0)
    @test point_on_spherical_arc(far, a2, b2) == false
end

@testset "Great circle arc intersection" begin
    using GeometryOps.UnitSpherical: spherical_arc_intersection, ArcIntersectionResult
    using GeometryOps.UnitSpherical: arc_cross, arc_hinge, arc_disjoint, arc_overlap

    # Two arcs that cross: equator vs prime meridian
    # Arc 1: (−45°,0°) to (45°,0°) along equator
    a1 = UnitSphereFromGeographic()((-45.0, 0.0))
    b1 = UnitSphereFromGeographic()((45.0, 0.0))
    # Arc 2: (0°,−45°) to (0°,45°) along prime meridian
    a2 = UnitSphereFromGeographic()((0.0, -45.0))
    b2 = UnitSphereFromGeographic()((0.0, 45.0))

    result = spherical_arc_intersection(a1, b1, a2, b2)
    @test result.type == arc_cross
    @test length(result.points) == 1
    # Intersection should be at (0°,0°)
    expected = UnitSphericalPoint(1.0, 0.0, 0.0)
    @test isapprox(result.points[1], expected, atol=1e-10)
    # Check fractions (α along arc1, β along arc2)
    @test isapprox(result.fracs[1][1], 0.5, atol=1e-10)  # midpoint of arc1
    @test isapprox(result.fracs[1][2], 0.5, atol=1e-10)  # midpoint of arc2

    # Two arcs that share an endpoint (hinge)
    a3 = UnitSphereFromGeographic()((0.0, 0.0))
    b3 = UnitSphereFromGeographic()((45.0, 0.0))
    a4 = UnitSphereFromGeographic()((0.0, 0.0))
    b4 = UnitSphereFromGeographic()((0.0, 45.0))

    result = spherical_arc_intersection(a3, b3, a4, b4)
    @test result.type == arc_hinge
    @test length(result.points) == 1
    @test isapprox(result.fracs[1][1], 0.0, atol=1e-10)  # start of arc1
    @test isapprox(result.fracs[1][2], 0.0, atol=1e-10)  # start of arc2

    # Two disjoint arcs
    a5 = UnitSphereFromGeographic()((0.0, 0.0))
    b5 = UnitSphereFromGeographic()((10.0, 0.0))
    a6 = UnitSphereFromGeographic()((20.0, 0.0))
    b6 = UnitSphereFromGeographic()((30.0, 0.0))

    result = spherical_arc_intersection(a5, b5, a6, b6)
    @test result.type == arc_disjoint
    @test isempty(result.points)

    # Overlapping collinear arcs
    a7 = UnitSphereFromGeographic()((0.0, 0.0))
    b7 = UnitSphereFromGeographic()((30.0, 0.0))
    a8 = UnitSphereFromGeographic()((20.0, 0.0))
    b8 = UnitSphereFromGeographic()((50.0, 0.0))

    result = spherical_arc_intersection(a7, b7, a8, b8)
    @test result.type == arc_overlap
    @test length(result.points) == 2
end