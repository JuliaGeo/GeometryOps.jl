using Test
using LinearAlgebra

import GeometryOps as GO, GeoInterface as GI
import Extents
using GeometryOps.UnitSpherical
using Random: Xoshiro

# k vertices of the z = z₀ circle, CCW seen from +z (the S2 interior-on-left
# convention: the enclosed region is the cap containing the north pole)
polar_ring(z, k) = [UnitSphericalPoint(sqrt(1 - z^2) * cos(t), sqrt(1 - z^2) * sin(t), z)
                    for t in range(0, 2π; length = k + 1)[1:(end - 1)]]

@testset "extent(Planar(), ...)" begin
    poly = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
    @test GO.extent(GO.Planar(), poly) == GI.extent(poly)
    @test GO.extent(GO.Planar(), GI.Point(1.0, 2.0)) == GI.extent(GI.Point(1.0, 2.0))
end

@testset "extent(Spherical(), ...)" begin
    m = GO.Spherical()
    z = 0.9; s = sqrt(1 - z^2)

    @testset "CCW polar cap ring encloses the pole" begin
        ext = GO.extent(m, GI.LinearRing(polar_ring(z, 8)))
        @test ext isa Extents.Extent{(:X, :Y, :Z)}
        @test ext.Z[2] == 1
        @test ext.Z[1] ≈ z atol = 1e-12
        @test -1 < ext.X[1] && ext.X[2] < 1
        @test -1 < ext.Y[1] && ext.Y[2] < 1
    end

    @testset "CW ring encloses the same cap (default is winding-independent)" begin
        rext = GO.extent(m, GI.LinearRing(reverse(polar_ring(z, 8))))
        @test rext == GO.extent(m, GI.LinearRing(polar_ring(z, 8)))
        @test rext.Z[2] == 1        # still the polar cap, not the complement
        @test rext.X[2] < 1 && rext.Y[2] < 1
    end

    @testset "Geographic polygon around the pole" begin
        cap = GI.Polygon([[(lon, 60.0) for lon in 0.0:30.0:360.0]])
        ext = GO.extent(m, cap)
        @test ext.Z[2] == 1
        @test ext.Z[1] ≈ sind(60) atol = 1e-12
    end

    @testset "No enclosure: region extent equals curve extent" begin
        pts = [(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 15.0), (5.0, 5.0)]
        @test GO.extent(m, GI.Polygon([pts])) == GO.extent(m, GI.LineString(pts))
    end

    @testset "Axis point on the boundary clamps, and extends nothing else" begin
        pts = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]  # corner at +eₓ
        ext = GO.extent(m, GI.Polygon([pts]))
        @test ext.X[2] == 1.0
        @test ext.X[1] > 0.9 && ext.Y[2] < 0.2 && ext.Z[2] < 0.2
        @test ext.Y[1] > -0.1 && ext.Z[1] > -0.1
    end

    @testset "Vertex exactly at the pole" begin
        ring = GI.LinearRing([UnitSphericalPoint(0.0, 0.0, 1.0),
                              UnitSphericalPoint(s, 0.0, z),
                              UnitSphericalPoint(0.0, s, z),
                              UnitSphericalPoint(0.0, 0.0, 1.0)])
        ext = GO.extent(m, ring)
        @test ext.Z[2] >= 1
        @test all(b -> all(isfinite, b), values(ext))
    end

    @testset "Ring below the equator encloses the southern cap" begin
        # the CCW ring at z = cosd(100) has the >hemisphere northern region
        # on its left; the default mode bounds the ENCLOSED southern cap
        # (the left region is `oriented = true` behavior, tested below)
        ext = GO.extent(m, GI.LinearRing(polar_ring(cosd(100), 16)))
        @test ext.Z[1] == -1                       # south pole enclosed
        @test ext.Z[2] ≈ cosd(100) atol = 1e-12    # rim on top; arcs bulge south
        @test -1 < ext.X[1] && ext.X[2] < 1
        @test -1 < ext.Y[1] && ext.Y[2] < 1
    end

    @testset "Dumbbell: both poles through a thin corridor" begin
        # two polar caps (above ±85°) joined by a corridor over lon ∈ (355°, 5°):
        # small area, zero net winding about the polar axis, and both poles and
        # (1, 0, 0) interior
        north = [(lon, 85.0) for lon in range(5.0, 355.0; length = 15)]    # eastward, pole on the left
        south = [(lon, -85.0) for lon in range(355.0, 5.0; length = 15)]   # westward, pole on the left
        ext = GO.extent(m, GI.Polygon([vcat(north, south, [north[1]])]))
        @test ext.Z == (-1.0, 1.0)
        @test ext.X[2] == 1.0       # (1, 0, 0) is inside the corridor
        @test ext.X[1] > -1         # (-1, 0, 0) is outside
        @test -0.2 < ext.Y[1] && ext.Y[2] < 0.2
    end

    @testset "Lonlat polar cells: exact pole vertex, no far-pole leak" begin
        for lon0 in 0.0:30.0:330.0
            cell = GI.Polygon([[(lon0, 80.0), (lon0 + 30.0, 80.0),
                                (lon0 + 30.0, 90.0), (lon0, 90.0), (lon0, 80.0)]])
            ext = GO.extent(m, cell)
            @test ext.Z[2] == 1.0   # the pole is a vertex
            @test ext.Z[1] > 0.9    # the south pole must not leak in
            @test ext.X[1] > -1 && ext.X[2] < 1 && ext.Y[1] > -1 && ext.Y[2] < 1
        end
    end

    @testset "Points, multis, and collections" begin
        p = GI.Point(0.0, 0.0)  # lon/lat → (1, 0, 0)
        ep = GO.extent(m, p)
        @test ep.X[1] == ep.X[2] == 1.0
        emp = GO.extent(m, GI.MultiPoint([(0.0, 0.0), (90.0, 0.0)]))
        @test emp.X == (0.0, 1.0) && emp.Y == (0.0, 1.0)

        north = GI.Polygon([[(lon, 60.0) for lon in 0.0:30.0:360.0]])
        south = GI.Polygon([[(lon, -60.0) for lon in 360.0:-30.0:0.0]])  # CCW around the south pole
        ext = GO.extent(m, GI.MultiPolygon([north, south]))
        @test ext.Z == (-1.0, 1.0)
    end

    @testset "LineStrings are curves, not regions" begin
        ext = GO.extent(m, GI.LineString(vcat(polar_ring(z, 8), [polar_ring(z, 8)[1]])))
        @test ext.Z[2] < 1  # no interior, no pole
    end

    @testset "Random cells contain their samples and enclosed axis points" begin
        rng = Xoshiro(2026)
        inext(q, e) = e.X[1] <= q[1] <= e.X[2] && e.Y[1] <= q[2] <= e.Y[2] && e.Z[1] <= q[3] <= e.Z[2]
        axispoints = [UnitSphericalPoint(v...) for v in
            ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))]
        for _ in 1:50
            c = rand(rng, UnitSphericalPoint{Float64})
            r = 0.05 + 0.95 * rand(rng)
            u = normalize(cross(c, abs(c.z) < 0.9 ? UnitSphericalPoint(0.0, 0.0, 1.0) : UnitSphericalPoint(1.0, 0.0, 0.0)))
            v = cross(c, u)
            k = rand(rng, 3:8)
            ring = [UnitSphericalPoint(cos(r) * c + sin(r) * (cos(t) * u + sin(t) * v))
                    for t in range(0, 2π; length = k + 1)[1:(end - 1)]]
            ext = GO.extent(m, GI.LinearRing(ring))
            # boundary and interior samples lie inside
            samples = [slerp(ring[i], ring[mod1(i + 1, k)], t)
                       for i in 1:k for t in range(0.0, 1.0; length = 33)]
            @test all(q -> inext(q, ext), samples)
            @test inext(c, ext)
            # the polygon is star-shaped around c, so it contains the cap
            # whose radius is the boundary's least distance to c — any axis
            # point in that cap must be covered
            r_in = 0.98 * minimum(q -> spherical_distance(c, q), samples)
            for a in axispoints
                spherical_distance(c, a) < r_in && @test inext(a, ext)
            end
            # the reversed ring encloses the same cap by default…
            @test GO.extent(m, GI.LinearRing(reverse(ring))) == ext
            # …and bounds the complement under `oriented = true`: it shares
            # the boundary, covers axis points outside the cap, and between
            # them the two boxes cover every axis point
            om = GO.Spherical(oriented = true)
            @test GO.extent(om, GI.LinearRing(ring)) == ext   # CCW: same region
            rext = GO.extent(om, GI.LinearRing(reverse(ring)))
            @test all(q -> inext(q, rext), samples)
            for a in axispoints
                spherical_distance(c, a) > r && @test inext(a, rext)
                @test inext(a, ext) || inext(a, rext)
            end
        end
    end
end

@testset "extent(Spherical(oriented = true), ...)" begin
    m = GO.Spherical()
    om = GO.Spherical(oriented = true)
    z = 0.9

    @testset "CCW polar cap matches the default mode" begin
        ring = GI.LinearRing(polar_ring(z, 8))
        @test GO.extent(om, ring) == GO.extent(m, ring)
        @test GO.extent(om, ring).Z[2] == 1
    end

    @testset "CW ring is the complement (left-of-ring region)" begin
        ext = GO.extent(om, GI.LinearRing(reverse(polar_ring(z, 8))))
        @test ext.Z[1] == -1
        @test ext.Z[2] < 1          # sup z of the complement is on the boundary
        @test ext.X == (-1.0, 1.0)  # complement contains ±eₓ and ±e_y
        @test ext.Y == (-1.0, 1.0)
    end

    @testset "Cap larger than a hemisphere" begin
        # the region on the CCW ring's left at z = cosd(100) reaches over
        # the north pole and past the equator
        ext = GO.extent(om, GI.LinearRing(polar_ring(cosd(100), 16)))
        @test ext.Z[2] == 1
        @test ext.X == (-1.0, 1.0)  # ±eₓ and ±e_y are interior
        @test ext.Y == (-1.0, 1.0)
        @test ext.Z[1] < cosd(100)  # arcs bulge below the vertex circle
    end

    @testset "Geographic CW polygon around the pole is the complement" begin
        cap_cw = GI.Polygon([[(lon, 60.0) for lon in 360.0:-30.0:0.0]])
        ext = GO.extent(om, cap_cw)
        @test ext.Z[1] == -1 && ext.Z[2] < 1
        # while the default mode reads it as the cap
        dext = GO.extent(m, cap_cw)
        @test dext.Z[2] == 1
        @test dext.Z[1] ≈ sind(60) atol = 1e-12
    end
end
