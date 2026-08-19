using Test
import GeometryOps as GO
import GeoInterface as GI
using GeometryOps.UnitSpherical: UnitSphereFromGeographic, UnitSphericalPoint

# Probe for the cache allocation test.  This has to be a top-level function: measuring
# `@allocated` directly on testset-local variables boxes them into the keyword `NamedTuple`
# and charges the algorithm bytes that it never actually allocates.
_area_bytes(alg, a, b, cache) = @allocated GO.intersection_area(alg, a, b; cache)

function _spherical_polygon(coords)
    transform = UnitSphereFromGeographic()
    points = [transform((lon, lat)) for (lon, lat) in coords]
    push!(points, points[1])
    return GI.Polygon([points])
end

square_a = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
square_b = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
small = GI.Polygon([[(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5), (0.5, 0.5)]])
far = GI.Polygon([[(9.0, 9.0), (10.0, 9.0), (10.0, 10.0), (9.0, 10.0), (9.0, 9.0)]])
edge_neighbour = GI.Polygon([[(2.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 0.0)]])

@testset "intersection_area" begin
    @testset "Planar - $(nameof(typeof(alg)))" for alg in
            (GO.ConvexConvexSutherlandHodgman(), GO.OverlayNG())
        @test GO.intersection_area(alg, square_a, square_b) ≈ 1.0
        @test GO.intersection_area(alg, square_a, small) ≈ 1.0
        @test GO.intersection_area(alg, small, square_a) ≈ 1.0
        @test GO.intersection_area(alg, square_a, square_a) ≈ 4.0
        @test GO.intersection_area(alg, square_a, far) == 0.0
        # sharing only an edge is a zero-area intersection
        @test GO.intersection_area(alg, square_a, edge_neighbour) == 0.0
        # and it agrees with building the polygon and measuring it
        for (a, b) in ((square_a, square_b), (square_a, small), (square_a, far))
            @test GO.intersection_area(alg, a, b) ≈ GO.area(GO.intersection(alg, a, b))
        end
    end

    @testset "Planar OverlayNG - non-convex inputs" begin
        alg = GO.OverlayNG()
        donut = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
                            [(3.0, 3.0), (3.0, 7.0), (7.0, 7.0), (7.0, 3.0), (3.0, 3.0)]])
        box = GI.Polygon([[(2.0, 2.0), (8.0, 2.0), (8.0, 8.0), (2.0, 8.0), (2.0, 2.0)]])
        # a hole survives into the result, and is subtracted
        @test GO.intersection_area(alg, donut, box) ≈ 20.0
        @test GO.intersection_area(alg, donut, box) ≈ GO.area(GO.intersection(alg, donut, box))

        # two disconnected result polygons
        bar = GI.Polygon([[(-1.0, 4.0), (11.0, 4.0), (11.0, 6.0), (-1.0, 6.0), (-1.0, 4.0)]])
        @test GO.intersection_area(alg, donut, bar) ≈ GO.area(GO.intersection(alg, donut, bar))

        # multipolygon input
        mp = GI.MultiPolygon([square_a, far])
        @test GO.intersection_area(alg, mp, square_b) ≈ GO.area(GO.intersection(alg, mp, square_b))

        # anything below dimension 2 has no area
        line = GI.LineString([(0.0, 0.0), (2.0, 2.0)])
        @test GO.intersection_area(alg, line, square_a) == 0.0
        @test GO.intersection_area(alg, GI.Point(1.0, 1.0), square_a) == 0.0
    end

    @testset "Sutherland-Hodgman rejects what `intersection` rejects" begin
        line = GI.LineString([(0.0, 0.0), (2.0, 2.0)])
        for alg in (GO.ConvexConvexSutherlandHodgman(), GO.ConvexConvexSutherlandHodgman(GO.Spherical()))
            @test_throws ArgumentError GO.intersection_area(alg, line, square_a)
            @test_throws ArgumentError GO.intersection_area(alg, square_a, GI.Point(1.0, 1.0))
        end
    end

    @testset "Float type" begin
        for alg in (GO.ConvexConvexSutherlandHodgman(), GO.OverlayNG())
            @test GO.intersection_area(alg, square_a, square_b, Float32) isa Float32
        end
        # the spherical radius must not silently widen the result back to Float64
        sph = GO.ConvexConvexSutherlandHodgman(GO.Spherical())
        sa = _spherical_polygon([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)])
        sb = _spherical_polygon([(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)])
        f32_cache = GO.SutherlandHodgmanCache(sph, Float32)
        @test GO.intersection_area(sph, sa, sb, Float32; cache = f32_cache) isa Float32
        @test GO.intersection_area(GO.OverlayNG(GO.Spherical()), square_a, square_b, Float32) isa Float32
    end

    @testset "Spherical" begin
        sh = GO.ConvexConvexSutherlandHodgman(GO.Spherical())
        ov = GO.OverlayNG(GO.Spherical())
        coords_a = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]
        coords_b = [(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)]
        coords_small = [(0.5, 0.5), (1.0, 0.5), (1.0, 1.0), (0.5, 1.0)]
        coords_far = [(50.0, 50.0), (52.0, 50.0), (52.0, 52.0), (50.0, 52.0)]

        for (ca, cb) in ((coords_a, coords_b), (coords_a, coords_small),
                         (coords_small, coords_a), (coords_a, coords_far))
            sa, sb = _spherical_polygon(ca), _spherical_polygon(cb)
            @test GO.intersection_area(sh, sa, sb) ≈
                GO.area(GO.Spherical(), GO.intersection(sh, sa, sb))

            la, lb = GI.Polygon([push!(collect(ca), first(ca))]), GI.Polygon([push!(collect(cb), first(cb))])
            @test GO.intersection_area(ov, la, lb) ≈
                GO.area(GO.Spherical(), GO.intersection(ov, la, lb))
            # both engines measure the same lune, whichever chart they work in
            @test GO.intersection_area(sh, sa, sb) ≈ GO.intersection_area(ov, la, lb) rtol=1e-8
        end

        # the lon/lat output chart reaches the same area as the xyz one it rounds from
        la = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
        lb = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
        lonlat = GO.OverlayNG(GO.Spherical(); point_type = Tuple{Float64,Float64})
        @test GO.intersection_area(lonlat, la, lb) ≈ GO.intersection_area(ov, la, lb) rtol=1e-9

        # a hole in the spherical result is subtracted, as it is on the plane
        donut = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
                            [(3.0, 3.0), (3.0, 7.0), (7.0, 7.0), (7.0, 3.0), (3.0, 3.0)]])
        box = GI.Polygon([[(2.0, 2.0), (8.0, 2.0), (8.0, 8.0), (2.0, 8.0), (2.0, 2.0)]])
        @test GO.intersection_area(ov, donut, box) ≈
            GO.area(GO.Spherical(), GO.intersection(ov, donut, box))

        # the radius of the manifold scales the result
        sa, sb = _spherical_polygon(coords_a), _spherical_polygon(coords_b)
        unit = GO.ConvexConvexSutherlandHodgman(GO.Spherical(; radius = 1.0))
        @test GO.intersection_area(sh, sa, sb) ≈
            GO.intersection_area(unit, sa, sb) * GO.Spherical().radius^2
    end

    @testset "Cache" begin
        planar = GO.ConvexConvexSutherlandHodgman()
        spherical = GO.ConvexConvexSutherlandHodgman(GO.Spherical())
        planar_cache = GO.SutherlandHodgmanCache(planar)
        spherical_cache = GO.SutherlandHodgmanCache(spherical)
        sa = _spherical_polygon([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)])
        sb = _spherical_polygon([(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)])

        @test GO.intersection_area(planar, square_a, square_b; cache = planar_cache) ≈
            GO.intersection_area(planar, square_a, square_b)
        @test GO.intersection_area(spherical, sa, sb; cache = spherical_cache) ≈
            GO.intersection_area(spherical, sa, sb)

        # a reused cache is the whole point: with one, nothing is allocated at all
        _area_bytes(planar, square_a, square_b, planar_cache)
        _area_bytes(spherical, sa, sb, spherical_cache)
        @test _area_bytes(planar, square_a, square_b, planar_cache) == 0
        @test _area_bytes(spherical, sa, sb, spherical_cache) == 0
        # Guard against passing vacuously: without a cache the buffers cost something
        _area_bytes(planar, square_a, square_b, nothing)
        @test _area_bytes(planar, square_a, square_b, nothing) > 0

        @test_throws ArgumentError GO.intersection_area(planar, square_a, square_b; cache = spherical_cache)
    end
end
