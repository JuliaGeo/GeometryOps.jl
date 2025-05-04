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