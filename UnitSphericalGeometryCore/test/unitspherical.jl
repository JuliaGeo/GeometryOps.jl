using Test
using LinearAlgebra
using UnitSphericalGeometryCore
import GeoInterface as GI
import Extents

@testset "spherical_distance" begin
    @testset "Basic correctness" begin
        # Orthogonal axes → π/2
        @test spherical_distance(UnitSphericalPoint(1.0, 0.0, 0.0),
                                 UnitSphericalPoint(0.0, 1.0, 0.0)) ≈ π/2
        # Identical points → 0
        p = UnitSphericalPoint(1.0, 0.0, 0.0)
        @test spherical_distance(p, p) == 0.0
        # Antipodal → π
        @test spherical_distance(UnitSphericalPoint(1.0, 0.0, 0.0),
                                 UnitSphericalPoint(-1.0, 0.0, 0.0)) ≈ π
        # Known π/4 separation
        @test spherical_distance(UnitSphericalPoint(1.0, 0.0, 0.0),
                                 UnitSphericalPoint(cos(π/4), sin(π/4), 0.0)) ≈ π/4
        # Known π/3 separation
        @test spherical_distance(UnitSphericalPoint(1.0, 0.0, 0.0),
                                 UnitSphericalPoint(cos(π/3), sin(π/3), 0.0)) ≈ π/3
    end

    @testset "Near-identical points (atan2 stability vs acos)" begin
        # For each ε, check that the atan2 form recovers ε with very small
        # relative error. This is the case where `acos(clamp(x⋅y, -1, 1))`
        # degrades — x⋅y = cos(ε) ≈ 1 - ε²/2, and double precision can't
        # distinguish that from 1 once ε ≲ 1e-8.
        a = UnitSphericalPoint(1.0, 0.0, 0.0)
        for ε in (1e-4, 1e-6, 1e-8, 1e-10, 1e-12)
            b = UnitSphericalPoint(cos(ε), sin(ε), 0.0)
            d = spherical_distance(a, b)
            @test d ≈ ε rtol=1e-12

            # Contrast: the naive acos form collapses to 0 around ε ≈ 1e-8.
            d_naive = acos(clamp(a ⋅ b, -1.0, 1.0))
            if ε ≤ 1e-10
                @test d_naive == 0.0 || abs(d_naive - ε) / ε > 1e-2
            end
        end
    end

    @testset "Near-antipodal points" begin
        # Should return ≈ π - ε. The atan2 form stays well-conditioned here
        # (cross norm ≈ sin(ε) ≈ ε, dot ≈ -1 + ε²/2).
        a = UnitSphericalPoint(1.0, 0.0, 0.0)
        for ε in (1e-2, 1e-4, 1e-6, 1e-8)
            # Point at angle (π - ε) from a, rotated in the xy-plane
            b = UnitSphericalPoint(cos(π - ε), sin(π - ε), 0.0)
            @test spherical_distance(a, b) ≈ π - ε rtol=1e-10
        end
    end

    # Ported from Google's S2 geometry library, TEST(S1Angle, ConstructorsThatMeasureAngles):
    # https://github.com/google/s2geometry/blob/a4f0cf58a9cfc214585c39de6e3682384fac0917/src/s2/s1angle_test.cc#L144-L152
    # S2's Vector3::Angle accepts non-unit vectors (the atan2 form is scale-invariant);
    # UnitSphericalPoint is typed for unit vectors, so we normalize the (0,0,2) case.
    @testset "s2geometry ConstructorsThatMeasureAngles" begin
        @test spherical_distance(UnitSphericalPoint(1.0, 0.0, 0.0),
                                 UnitSphericalPoint(0.0, 0.0, 1.0)) == π/2
        @test spherical_distance(UnitSphericalPoint(1.0, 0.0, 0.0),
                                 UnitSphericalPoint(1.0, 0.0, 0.0)) == 0.0
    end

    @testset "Symmetry and range" begin
        using Random
        Random.seed!(0xA7A7)
        for _ in 1:50
            a = UnitSphericalPoint(randn(), randn(), randn())
            b = UnitSphericalPoint(randn(), randn(), randn())
            a = a / norm(a); b = b / norm(b)
            d_ab = spherical_distance(a, b)
            d_ba = spherical_distance(b, a)
            @test d_ab == d_ba
            @test 0 ≤ d_ab ≤ π + 1e-12
        end
    end
end

@testset "slerp" begin
    @testset "Basic correctness" begin
        # Quarter-turn along the equator
        @test slerp(UnitSphericalPoint(1.0, 0.0, 0.0),
                    UnitSphericalPoint(0.0, 1.0, 0.0), 0.5) ≈
              UnitSphericalPoint(sqrt(2)/2, sqrt(2)/2, 0.0)
        # Endpoints
        a = UnitSphericalPoint(1.0, 0.0, 0.0)
        b = UnitSphericalPoint(0.0, 1.0, 0.0)
        @test slerp(a, b, 0.0) ≈ a
        @test slerp(a, b, 1.0) ≈ b
        # Identical inputs — fast path
        @test slerp(a, a, 0.5) == a
    end

    @testset "Random well-separated points stay on the unit sphere" begin
        using Random
        Random.seed!(0x5E5E)
        for _ in 1:50
            va = randn(3); va ./= norm(va)
            vb = randn(3); vb ./= norm(vb)
            # Skip anything within 1° of antipodal — that's the degenerate regime
            # covered by the dedicated testsets below.
            (va ⋅ vb) < -0.9998 && continue
            a = UnitSphericalPoint(va[1], va[2], va[3])
            b = UnitSphericalPoint(vb[1], vb[2], vb[3])
            for t in (0.0, 0.25, 0.5, 0.75, 1.0)
                @test isapprox(norm(slerp(a, b, t)), 1.0, atol=1e-12)
            end
        end
    end

    @testset "Well-conditioned near-antipodal geometry" begin
        # When the perpendicular direction is captured cleanly (b lives
        # exactly in the xz-plane, a on the x-axis), slerp stays accurate
        # right up against π.
        a = UnitSphericalPoint(1.0, 0.0, 0.0)
        for θ_deg in (179.0, 179.99, 179.9999, 179.99999999)
            θ = deg2rad(θ_deg)
            b = UnitSphericalPoint(cos(θ), 0.0, sin(θ))
            @test isapprox(norm(slerp(a, b, 0.5)), 1.0, atol=1e-10)
        end
    end

    # Previously broken under the divide-by-sin(Ω) form; now pass via the
    # tangent-vector form in src/utils/UnitSpherical/slerp.jl.
    @testset "Cancellation-prone near-antipodal" begin
        a = UnitSphericalPoint(1.0, 0.0, 0.0)
        for ε in (1e-10, 1e-12, 1e-14, 1e-15, 1e-16)
            v = [-1.0, 0.0, ε]; v ./= norm(v)
            b = UnitSphericalPoint(v[1], v[2], v[3])
            m = slerp(a, b, 0.5)
            @test isapprox(norm(m), 1.0, atol=1e-10)
        end
    end

    @testset "Exactly antipodal returns a unit vector" begin
        a = UnitSphericalPoint(1.0, 0.0, 0.0)
        b = UnitSphericalPoint(-1.0, 0.0, 0.0)
        m = slerp(a, b, 0.5)
        @test all(isfinite, m)
        @test isapprox(norm(m), 1.0, atol=1e-10)
    end

    # Regression test for ConservativeRegridding.jl#83. A half-sphere
    # equatorial band has a truly antipodal diagonal between
    # (0°, −15°) and (180°, +15°); under the old slerp this produced
    # a (0,0,0) midpoint that poisoned the SphericalCap radius.
    @testset "ConservativeRegridding#83 antipodal diagonal" begin
        to_sphere = UnitSphereFromGeographic()
        p1 = to_sphere((0.0, -15.0))
        p3 = to_sphere((180.0, +15.0))
        @test p1 ⋅ p3 ≈ -1.0
        m = slerp(p1, p3, 0.5)
        @test isapprox(norm(m), 1.0, atol=1e-10)
    end
end

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
        @test Extents.intersects(cap1, cap2)
        @test Extents.intersects(cap2, cap1)
        @test !UnitSphericalGeometryCore._disjoint(cap1, cap2)

        # Create two caps that don't intersect
        cap3 = SphericalCap(UnitSphericalPoint(1.0, 0.0, 0.0), π/8)
        cap4 = SphericalCap(UnitSphericalPoint(0.0, 0.0, 1.0), π/8)
        @test !Extents.intersects(cap3, cap4)
        @test UnitSphericalGeometryCore._disjoint(cap3, cap4)

        # Test containment
        big_cap = SphericalCap(UnitSphericalPoint(1.0, 0.0, 0.0), π/2)
        small_cap = SphericalCap(UnitSphericalPoint(1/√2, 1/√2, 0.0), π/4)
        @test Extents.contains(big_cap, small_cap)
        @test Extents.within(small_cap, big_cap)
        @test Extents.contains(big_cap, small_cap; strict=true)
        @test !Extents.contains(small_cap, big_cap)
        @test !Extents.within(big_cap, small_cap)
    end

    @testset "ULP-near contact and coincident centres" begin
        centre = UnitSphericalPoint(1.0, 0.0, 0.0)
        for θ in (1e-10, 0.75, π / 2)
            other = UnitSphericalPoint(cos(θ), sin(θ), 0.0)
            contact = spherical_distance(centre, other)
            pointcap = SphericalCap(other, 0.0)
            @test !Extents.intersects(SphericalCap(centre, prevfloat(contact)), pointcap)
            @test Extents.intersects(SphericalCap(centre, contact), pointcap)
            @test Extents.intersects(SphericalCap(centre, nextfloat(contact)), pointcap)
        end

        pointcap = SphericalCap(centre, 0.0)
        same = SphericalCap(centre, 0.0)
        wider = SphericalCap(centre, 0.25)
        @test Extents.intersects(pointcap, same)
        @test Extents.contains(pointcap, same)
        @test Extents.contains(wider, pointcap)
        @test !Extents.contains(pointcap, wider)
        @test Extents.union(pointcap, same) === pointcap
        @test Extents.union(wider, pointcap) === wider
        @test Extents.union(pointcap, wider) === wider
        @test Extents.union(pointcap, wider; strict=true) === wider

        whole = SphericalCap(centre, π)
        antipode = UnitSphericalPoint(-1.0, 0.0, 0.0)
        @test UnitSphericalGeometryCore._contains(whole, antipode)
        @test Extents.contains(whole, SphericalCap(antipode, π / 3))
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
            UnitSphericalGeometryCore.circumcenter_on_unit_sphere(a, b, c),
            UnitSphericalGeometryCore.circumcenter_on_unit_sphere(a, c, b),
            UnitSphericalGeometryCore.circumcenter_on_unit_sphere(b, a, c),
            UnitSphericalGeometryCore.circumcenter_on_unit_sphere(b, c, a),
            UnitSphericalGeometryCore.circumcenter_on_unit_sphere(c, a, b),
            UnitSphericalGeometryCore.circumcenter_on_unit_sphere(c, b, a),
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
            d1 = UnitSphericalGeometryCore.spherical_distance(cap.point, points[1])
            d2 = UnitSphericalGeometryCore.spherical_distance(cap.point, points[2])
            d3 = UnitSphericalGeometryCore.spherical_distance(cap.point, points[3])
            @test isapprox(d1, d2, atol=1e-10)
            @test isapprox(d2, d3, atol=1e-10)
            @test isapprox(d1, cap.radius, atol=1e-10)
        end
    end

    @testset "Union and growth of SphericalCaps" begin
        function test_merge(p1, p2, r1, r2, pmerged, rmerged)
            r1 = deg2rad(r1)
            r2 = deg2rad(r2)
            cap1 = SphericalCap(p1, r1)
            cap2 = SphericalCap(p2, r2)
            capmerged = Extents.union(cap1, cap2)
            @test all(isapprox.(GeographicFromUnitSphere()(capmerged.point), pmerged))
            @test isapprox(capmerged.radius, deg2rad(rmerged))
            @test Extents.contains(capmerged, cap1)
            @test Extents.contains(capmerged, cap2)
        end

        test_merge((10.0, 0.0), (30.0, 0.0), 5, 10, (22.5, 0.0), 17.5)
        test_merge((10.0, 0.0), (30.0, 0.0), 15, 10, (17.5, 0.0), 22.5)
        test_merge((10.0, 0.0), (30.0, 0.0), 40, 5, (10.0, 0.0), 40.0)
        test_merge((10.0, 0.0), (30.0, 0.0), 5, 50, (30.0, 0.0), 50.0)

        east = UnitSphericalPoint(1.0, 0.0, 0.0)
        west = UnitSphericalPoint(-1.0, 0.0, 0.0)
        antipodal_merge = Extents.union(SphericalCap(east, 0.0), SphericalCap(west, 0.0))
        @test UnitSphericalGeometryCore._contains(antipodal_merge, east)
        @test UnitSphericalGeometryCore._contains(antipodal_merge, west)
        wide_merge = Extents.union(SphericalCap(east, 3π / 4), SphericalCap(west, 3π / 4))
        @test wide_merge.radius == Float64(π)
        @test Extents.contains(wide_merge, SphericalCap(east, 3π / 4))
        @test Extents.contains(wide_merge, SphericalCap(west, 3π / 4))

        using Random: Xoshiro
        rng = Xoshiro(0xcafecafe)
        for _ in 1:100
            a = SphericalCap(rand(rng, UnitSphericalPoint{Float64}), 0.4 * rand(rng))
            b = SphericalCap(rand(rng, UnitSphericalPoint{Float64}), 0.4 * rand(rng))
            merged = Extents.union(a, b)
            @test Extents.contains(merged, a)
            @test Extents.contains(merged, b)
        end

        cap = SphericalCap(UnitSphericalPoint(1.0, 0.0, 0.0), 0.2)
        @test Extents.grow(cap, 0.0) === cap
        dilated = Extents.grow(cap, 0.75)
        @test dilated.radius > 0.5
        @test Extents.contains(dilated, cap)
        @test Extents.grow(cap, 0.5).radius > 2 * cap.radius
        @test Extents.grow(cap, 10.0).radius == Float64(π)
        @test Extents.grow(cap, -0.25).radius == cap.radius / 2
        @test_throws DomainError Extents.grow(cap, -0.75)
    end
end


@testset "Spherical cap XYZ bounds" begin
    in_extent(p, ext) =
        ext.X[1] <= p[1] <= ext.X[2] &&
        ext.Y[1] <= p[2] <= ext.Y[2] &&
        ext.Z[1] <= p[3] <= ext.Z[2]

    function boundary_points(cap, n = 72)
        c = collect(cap.point)
        axis = abs(c[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
        u = normalize(cross(c, axis))
        v = cross(c, u)
        return [UnitSphericalPoint(Tuple(
            cos(cap.radius) .* c .+
            sin(cap.radius) .* (cos(ϕ) .* u .+ sin(ϕ) .* v)))
            for ϕ in range(0.0, 2π; length = n + 1)[1:end-1]]
    end

    tiny = SphericalCap(UnitSphericalPoint(1.0, 0.0, 0.0), 1e-12)
    @test which(Extents.extent, (typeof(tiny),)) == which(Extents.extent, (Any,))
    @test Extents.extent(tiny) === nothing
    tiny_ext = convert(Extents.Extent, tiny)
    @test tiny_ext isa Extents.Extent{(:X, :Y, :Z)}
    @test convert(Extents.Extent{(:X, :Y, :Z)}, tiny) === tiny_ext
    @test convert(typeof(tiny_ext), tiny) === tiny_ext
    @test_throws ArgumentError convert(Extents.Extent{(:X, :Y)}, tiny)
    @test tiny_ext.X[1] <= cos(tiny.radius) <= tiny_ext.X[2]
    @test tiny_ext.Y[1] <= -sin(tiny.radius) < sin(tiny.radius) <= tiny_ext.Y[2]
    @test tiny_ext.Z[1] <= -sin(tiny.radius) < sin(tiny.radius) <= tiny_ext.Z[2]

    polar = SphericalCap(UnitSphericalPoint(0.0, 0.0, 1.0), 0.2)
    polar_ext = convert(Extents.Extent, polar)
    @test polar_ext.X[1] <= -sin(polar.radius) < sin(polar.radius) <= polar_ext.X[2]
    @test polar_ext.Y[1] <= -sin(polar.radius) < sin(polar.radius) <= polar_ext.Y[2]
    @test polar_ext.Z[1] <= cos(polar.radius) <= polar_ext.Z[2] == 1.0

    antimeridian = SphericalCap(UnitSphereFromGeographic()((180.0, 0.0)), 0.3)
    antimeridian_ext = convert(Extents.Extent, antimeridian)
    @test antimeridian_ext.X[1] == -1.0
    @test antimeridian_ext.X[1] <= -cos(antimeridian.radius) <= antimeridian_ext.X[2]
    @test antimeridian_ext.Y[1] <= -sin(antimeridian.radius) < sin(antimeridian.radius) <= antimeridian_ext.Y[2]

    hemisphere = SphericalCap(UnitSphericalPoint(0.0, 0.0, 1.0), π / 2)
    hemisphere_ext = convert(Extents.Extent, hemisphere)
    @test hemisphere_ext.X == (-1.0, 1.0)
    @test hemisphere_ext.Y == (-1.0, 1.0)
    @test hemisphere_ext.Z[1] <= 0.0 <= hemisphere_ext.Z[2] == 1.0

    whole = SphericalCap(UnitSphericalPoint(0.0, 0.0, 1.0), nextfloat(Float64(π)))
    whole_ext = convert(Extents.Extent, whole)
    @test whole_ext == Extents.Extent(X = (-1.0, 1.0), Y = (-1.0, 1.0), Z = (-1.0, 1.0))

    arbitrary = SphericalCap(normalize(UnitSphericalPoint(1.0, -2.0, 3.0)), 0.7)
    for cap in (tiny, polar, antimeridian, hemisphere, arbitrary)
        ext = convert(Extents.Extent, cap)
        @test in_extent(cap.point, ext)
        @test all(p -> in_extent(p, ext), boundary_points(cap))
        @test Extents.intersects(cap, ext)
        @test Extents.intersects(cap, ext) == Extents.intersects(ext, cap)
    end

    cap32 = SphericalCap(UnitSphericalPoint(0.0f0, 0.0f0, 1.0f0), 0.2f0)
    ext32 = convert(Extents.Extent, cap32)
    @test ext32.X isa Tuple{Float32, Float32}
    @test all(p -> in_extent(p, ext32), boundary_points(cap32))
end

@testset "Cap–extent intersection" begin
    cap = SphericalCap(UnitSphericalPoint(1.0, 0.0, 0.0), π/4)
    around_center = Extents.Extent(X = (0.9, 1.1), Y = (-0.1, 0.1), Z = (-0.1, 0.1))
    antipodal = Extents.Extent(X = (-1.1, -0.9), Y = (-0.1, 0.1), Z = (-0.1, 0.1))
    @test Extents.intersects(cap, around_center)
    @test Extents.intersects(around_center, cap)
    @test !Extents.intersects(cap, antipodal)
    @test !Extents.intersects(antipodal, cap)

    # Degenerate (point) extents straddling the cap boundary
    just_inside = UnitSphericalPoint(cos(π/4 - 1e-3), sin(π/4 - 1e-3), 0.0)
    just_outside = UnitSphericalPoint(cos(π/4 + 1e-3), sin(π/4 + 1e-3), 0.0)
    @test Extents.intersects(cap, GI.extent(just_inside))
    @test !Extents.intersects(cap, GI.extent(just_outside))

    # A whole-sphere cap reaches the antipodal box; radii past π clamp to full
    @test Extents.intersects(SphericalCap(UnitSphericalPoint(1.0, 0.0, 0.0), π), antipodal)
    @test Extents.intersects(SphericalCap(UnitSphericalPoint(1.0, 0.0, 0.0), 3π/2), antipodal)

    # Tiny radii keep full precision (the chord radius is computed as
    # `2sin(radius/2)`, not from `radiuslike = cos(radius)`, which rounds
    # to 1 below radius ≈ 1.5e-8 and would empty the cap)
    tinycap = SphericalCap(UnitSphericalPoint(1.0, 0.0, 0.0), 1e-9)
    @test Extents.intersects(tinycap, GI.extent(UnitSphericalPoint(cos(0.9e-9), sin(0.9e-9), 0.0)))
    @test !Extents.intersects(tinycap, GI.extent(UnitSphericalPoint(cos(1.1e-9), sin(1.1e-9), 0.0)))
end

@testset "Points from geometries" begin
    @testset "3D GeoInterface points are unit xyz, 2D are lon/lat" begin
        @test UnitSphericalPoint(GI.Point(0.0, 0.0, 1.0)) == UnitSphericalPoint(0.0, 0.0, 1.0)
        @test UnitSphericalPoint(GI.Point(0.0, 90.0)) ≈ UnitSphericalPoint(0.0, 0.0, 1.0)
        @test UnitSphericalPoint((0.0, 0.0, 1.0)) == UnitSphericalPoint(0.0, 0.0, 1.0)
        @test UnitSphericalPoint([0.0, 90.0]) ≈ UnitSphericalPoint(0.0, 0.0, 1.0)
    end

    @testset "to_unit_spherical_points" begin
        ring = GI.LinearRing([(0.0, 0.0), (90.0, 0.0), (0.0, 90.0), (0.0, 0.0)])
        pts = to_unit_spherical_points(ring)
        @test pts isa Vector{UnitSphericalPoint{Float64}}
        @test length(pts) == 4
        @test pts[1] ≈ UnitSphericalPoint(1.0, 0.0, 0.0)
        @test pts[2] ≈ UnitSphericalPoint(0.0, 1.0, 0.0)
        @test pts[3] ≈ UnitSphericalPoint(0.0, 0.0, 1.0)
        # existing unit spherical points pass through
        @test to_unit_spherical_points(GI.LineString(pts)) == pts
    end
end

@testset "spherical_bounding_cap" begin
    inside(cap, p) = spherical_distance(cap.point, p) <= cap.radius

    @testset "Point" begin
        cap = spherical_bounding_cap(GI.Point(10.0, 20.0))
        @test cap isa SphericalCap{Float64}
        @test cap.point ≈ UnitSphereFromGeographic()((10.0, 20.0))
        @test 0 < cap.radius < 1e-6
        @test SphericalCap(GI.Point(10.0, 20.0)).radius == 0
    end

    @testset "Small geometries contain their vertices and edges" begin
        line = GI.LineString([(0.0, 0.0), (10.0, 0.0), (10.0, 10.0)])
        cap = spherical_bounding_cap(line)
        @test cap.radius < π / 2
        pts = to_unit_spherical_points(line)
        @test all(p -> inside(cap, p), pts)
        for i in 1:length(pts)-1, t in range(0, 1; length = 51)
            @test inside(cap, slerp(pts[i], pts[i+1], t))
        end
        # a polygon's cap is built from the vertices of all of its rings
        poly = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
                           [(2.0, 2.0), (3.0, 2.0), (3.0, 3.0), (2.0, 2.0)]])
        pcap = spherical_bounding_cap(poly)
        @test all(p -> inside(pcap, p), to_unit_spherical_points(poly))
        @test SphericalCap(poly) == pcap
        # 3D (unit xyz) input is taken as-is
        xyz = GI.LineString([GI.Point(1.0, 0.0, 0.0), GI.Point(0.0, 1.0, 0.0)])
        @test spherical_bounding_cap(xyz).point ≈ UnitSphericalPoint(normalize([1.0, 1.0, 0.0]))
    end

    @testset "Degenerate and wide geometries give the whole sphere" begin
        # vertices spread over the equator: the vertex mass is ~0
        equator = GI.LineString([(0.0, 0.0), (90.0, 0.0), (180.0, 0.0), (-90.0, 0.0)])
        @test spherical_bounding_cap(equator).radius == Float64(π)
        # vertices more than a quarter turn from their mass center
        wide = GI.LineString([(0.0, 0.0), (0.0, 80.0), (0.0, -80.0), (170.0, 0.0)])
        @test spherical_bounding_cap(wide).radius == Float64(π)
        whole = spherical_bounding_cap(equator)
        @test Extents.intersects(whole, spherical_bounding_cap(GI.Point(0.0, 0.0)))
    end

    @testset "Disjoint caps reject, overlapping caps do not" begin
        a = spherical_bounding_cap(GI.LineString([(0.0, 0.0), (1.0, 1.0)]))
        b = spherical_bounding_cap(GI.LineString([(10.0, 10.0), (11.0, 11.0)]))
        c = spherical_bounding_cap(GI.LineString([(0.5, 0.5), (12.0, 12.0)]))
        @test !Extents.intersects(a, b)
        @test Extents.intersects(a, c) && Extents.intersects(b, c)
    end
end
