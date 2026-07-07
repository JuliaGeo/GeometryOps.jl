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

    @testset "CW ring is the complement (S2 convention)" begin
        ext = GO.extent(m, GI.LinearRing(reverse(polar_ring(z, 8))))
        @test ext.Z[1] == -1
        @test ext.Z[2] < 1          # sup z of the complement is on the boundary
        @test ext.X == (-1.0, 1.0)  # complement contains ±eₓ and ±e_y
        @test ext.Y == (-1.0, 1.0)
    end

    @testset "Geographic polygon around the pole" begin
        cap = GI.Polygon([[(lon, 60.0) for lon in 0.0:30.0:360.0]])
        ext = GO.extent(m, cap)
        @test ext.Z[2] == 1
        @test ext.Z[1] ≈ sind(60) atol = 1e-12
    end

    @testset "No enclosure: region extent equals curve extent" begin
        pts = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]
        @test GO.extent(m, GI.Polygon([pts])) == GO.extent(m, GI.LineString(pts))
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

    @testset "Cap larger than a hemisphere" begin
        ext = GO.extent(m, GI.LinearRing(polar_ring(cosd(100), 16)))
        @test ext.Z[2] == 1
        @test ext.X == (-1.0, 1.0)  # ±eₓ and ±e_y are interior
        @test ext.Y == (-1.0, 1.0)
        @test ext.Z[1] < cosd(100)  # arcs bulge below the vertex circle
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
        end
    end
end
