using Test

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
        @test UnitSpherical._contains(big_cap, small_cap)
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
end