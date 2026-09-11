using Test
import GeometryOps as GO
import GeoInterface as GI

@testset "spherical proper crossing nodes" begin
    # The meridian has zero XY displacement in kernel coordinates: applying
    # the planar crossing solve here divides by zero.
    fixtures = (
        ([(0., 0.), (10., 0.)], [(5., -1.), (5., 1.)]),
        ([(170., 0.), (-170., 0.)], [(180., -5.), (180., 5.)]),
        ([(0., 35.), (10., 45.)], [(0., 45.), (10., 35.)]),
    )
    for (a, b) in fixtures, reverse_a in (false, true), reverse_b in (false, true),
            swap in (false, true), exact in (GO.True(), GO.False()),
            accelerator in (GO.NestedLoop(), GO.DoubleNaturalTree())
        A = GI.LineString(reverse_a ? reverse(a) : a)
        B = GI.LineString(reverse_b ? reverse(b) : b)
        A, B = swap ? (B, A) : (A, B)
        alg = GO.RelateNG(; manifold=GO.Spherical(), exact, accelerator)
        prepared = GO.prepare(alg, A)
        @test GO.relate_predicate(prepared, GO.pred_intersects(), B)
        @test GO.relate(alg, A, B, "0F1FF0102")
        @test GO.relate_predicate(prepared, GO.pred_crosses(), B)
        @test !GO.relate_predicate(prepared, GO.pred_touches(), B)
        @test !GO.relate_predicate(prepared, GO.pred_disjoint(), B)
    end
end

@testset "spherical crossing representatives retain direction and scale" begin
    m = GO.Spherical()
    # Nonzero Z and both antipodal candidates. Tiny arcs have a rational
    # crossing direction too small to convert to Float64 before scaling.
    for t in (0.5, 1e-200), sign in (-1., 1.), axis in (1, 3)
        xyz = ((sign, t, 0.), (sign, -t, 0.), (sign, 0., t), (sign, 0., -t))
        pts = map(xyz) do p
            q = axis == 1 ? p : (p[3], p[2], p[1])
            GO.rk_normalize_usp(GO.UnitSphericalPoint(q...))
        end
        expected = axis == 1 ? (sign, 0., 0.) : (0., 0., sign)
        for (a, b, c, d) in (pts, (pts[2], pts[1], pts[3], pts[4]),
                (pts[3], pts[4], pts[1], pts[2]))
            key = GO.crossing_node(a, b, c, d)
            p = GO._crossing_locate_point(m, key)
            @test p isa GO.UnitSphericalPoint{Float64}
            @test (GI.x(p), GI.y(p), GI.z(p)) == expected
        end
    end
end

@testset "crossing node location in mixed spherical collections" begin
    m = GO.Spherical()
    equator = GI.LineString([(0., 0.), (10., 0.)])
    meridian = GI.LineString([(5., -1.), (5., 1.)])
    pts = map(p -> GO._to_kernel_point(m, p), ((0., 0.), (10., 0.), (5., -1.), (5., 1.)))
    key = GO.crossing_node(pts...)
    for (polygon, covered) in (
        (GI.Polygon([[(3., -2.), (7., -2.), (7., 2.), (3., 2.), (3., -2.)]]), true),
        (GI.Polygon([[(20., 20.), (25., 20.), (25., 25.), (20., 25.), (20., 20.)]]), false),
    )
        collection = GI.GeometryCollection([equator, polygon])
        rg = GO.RelateGeometry(m, collection; exact=GO.True())
        @test GO.locate_node(rg, key, nothing) == GO.LOC_INTERIOR
        @test GO.is_node_in_area(rg, key, nothing) == covered
        alg = GO.RelateNG(; manifold=m)
        @test GO.relate_predicate(GO.prepare(alg, collection), GO.pred_intersects(), meridian)
    end
    # A line crossing a polygon boundary also materializes a symbolic node.
    box = GI.Polygon([[(3., -2.), (7., -2.), (7., 2.), (3., 2.), (3., -2.)]])
    @test GO.relate(GO.RelateNG(; manifold=m), equator, box, "T*T******")

    @test GO.relate(equator, meridian, "0F1FF0102")
end
