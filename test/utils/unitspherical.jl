using Test
using LinearAlgebra

using GeometryOps.UnitSpherical

import GeometryOps as GO
import GeoInterface as GI
import Extents
import UnitSphericalGeometryCore

@testset "GeometryOps.UnitSpherical builds on UnitSphericalGeometryCore" begin
    @test GO.UnitSpherical.UnitSphericalPoint === UnitSphericalGeometryCore.UnitSphericalPoint
    @test GO.UnitSpherical.SphericalCap === UnitSphericalGeometryCore.SphericalCap
    @test GO.UnitSpherical.slerp === UnitSphericalGeometryCore.slerp
    @test GO.UnitSpherical.spherical_bounding_cap === UnitSphericalGeometryCore.spherical_bounding_cap
    @test GO.UnitSpherical.robust_cross_product === UnitSphericalGeometryCore.robust_cross_product
    @test GO.UnitSpherical.RobustCrossProduct === UnitSphericalGeometryCore.RobustCrossProduct
    # UnitSphericalGeometryCore has no predicates; those are GeometryOps' own
    @test !isdefined(UnitSphericalGeometryCore, :spherical_orient)
    @test isdefined(GO.UnitSpherical, :spherical_orient)
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


@testset "Spatial tree pruning with SphericalCap against brute force" begin
    using Random: Xoshiro
    rng = Xoshiro(1312)
    points = rand(rng, UnitSphericalPoint{Float64}, 1000)
    tree = GO.FlexibleRTrees.RTree(GO.FlexibleRTrees.HPR(), points)
    for radius in (0.05, 0.5, 2.0, π)
        querycap = SphericalCap(rand(rng, UnitSphericalPoint{Float64}), radius)
        hits = sort!(GO.SpatialTreeInterface.depth_first_search(
            Base.Fix1(Extents.intersects, querycap), tree))
        @test hits == findall(p -> spherical_distance(querycap.point, p) <= querycap.radius, points)
    end
end

@testset "spherical_arc_extent" begin
    # Quarter arc along the equator: x and y peak at the endpoints, z is flat
    ext = spherical_arc_extent(UnitSphericalPoint(1.0, 0.0, 0.0), UnitSphericalPoint(0.0, 1.0, 0.0))
    @test ext isa Extents.Extent{(:X, :Y, :Z)}
    @test ext.X[1] ≈ 0 atol = 1e-14
    @test ext.X[2] ≈ 1 atol = 1e-14
    @test ext.Y[1] ≈ 0 atol = 1e-14
    @test ext.Y[2] ≈ 1 atol = 1e-14
    @test ext.Z[1] ≈ 0 atol = 1e-14
    @test ext.Z[2] ≈ 0 atol = 1e-14

    # Two points at z = 0.9 joined over the pole: the arc reaches z = 1
    z = 0.9; s = sqrt(1 - z^2)
    bulging = spherical_arc_extent(UnitSphericalPoint(0.0, -s, z), UnitSphericalPoint(0.0, s, z))
    @test bulging.Z[2] ≈ 1 atol = 1e-14
    @test bulging.Z[2] > z
    @test bulging.Z[1] ≈ z atol = 1e-14
    @test bulging.Y[1] ≈ -s atol = 1e-14
    @test bulging.Y[2] ≈ s atol = 1e-14

    # Geographic (lon, lat) input
    @test spherical_arc_extent((0.0, 0.0), (90.0, 0.0)) ==
        spherical_arc_extent(UnitSphericalPoint(1.0, 0.0, 0.0), UnitSphericalPoint(0.0, 1.0, 0.0))

    # Degenerate arc: a point
    p = UnitSphericalPoint(1.0, 0.0, 0.0)
    degenerate = spherical_arc_extent(p, p)
    @test degenerate.X[1] <= 1 <= degenerate.X[2]
    @test degenerate.X[2] - degenerate.X[1] < 1e-14

    @testset "Containment of dense arc samples" begin
        using Random: Xoshiro
        rng = Xoshiro(999)
        inext(q, e) = e.X[1] <= q[1] <= e.X[2] && e.Y[1] <= q[2] <= e.Y[2] && e.Z[1] <= q[3] <= e.Z[2]
        for _ in 1:100
            a, b = rand(rng, UnitSphericalPoint{Float64}), rand(rng, UnitSphericalPoint{Float64})
            e = spherical_arc_extent(a, b)
            @test all(t -> inext(slerp(a, b, t), e), range(0.0, 1.0, length = 101))
            # nearly-degenerate arc
            b2 = UnitSphericalPoint(normalize(a + 1e-12 * rand(rng, UnitSphericalPoint{Float64})))
            e2 = spherical_arc_extent(a, b2)
            @test all(t -> inext(slerp(a, b2, t), e2), range(0.0, 1.0, length = 11))
        end
    end
end
