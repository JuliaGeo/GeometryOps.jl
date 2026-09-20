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

@testset "exact spherical coincidence classes" begin
    m = GO.Spherical()
    rg = GO.RelateGeometry(m, GI.LineString([(0., 0.), (10., 0.)]); exact=GO.True())
    for axis in 1:3, t in (0.5, 1e-200), sign in (-1., 1.)
        point(x, y, z) = GO.UnitSphericalPoint(circshift([x, y, z], axis - 1)...)
        a, b = point(sign, t, 0.), point(sign, -t, 0.)
        c, d = point(sign, 0., t), point(sign, 0., -t)
        e, f = point(sign, t, t), point(sign, -t, -t)
        crossing1 = GO.crossing_node(a, b, c, d)
        crossing2 = GO.crossing_node(a, b, e, f)
        vertex = GO.vertex_node(point(sign, 0., 0.))
        antipode = GO.vertex_node(point(-sign, 0., 0.))
        near = GO.vertex_node(point(sign, eps(), 0.))
        identity = GO._spherical_node_identity(crossing1)
        @test identity == GO._spherical_node_identity(vertex)
        @test identity == GO._spherical_node_identity(crossing2)
        @test identity != GO._spherical_node_identity(antipode)
        @test identity != GO._spherical_node_identity(near)
        for keys in ((crossing1, crossing2, vertex, antipode, near),
                (near, antipode, vertex, crossing2, crossing1))
            tc = GO.TopologyComputer(GO.RelateMatrixPredicate(), rg, rg)
            for k in keys
                sections = GO.NodeSections(k)
                section = GO.NodeSection(true, GO.DIM_L, Int32(1), Int32(0),
                    nothing, !k.is_crossing, a, k, b)
                GO.add_node_section!(sections, section)
                tc.node_sections[k] = sections
            end
            GO._merge_coincident_nodes!(tc)
            @test length(tc.node_sections) == 3
            @test haskey(tc.node_sections, vertex)
            @test length(tc.node_sections[vertex].sections) == 3
            @test all(s -> s.node == vertex, tc.node_sections[vertex].sections)
        end
    end

    # Direction identity also applies when there are no crossing keys.
    vertex = GO.vertex_node(GO.UnitSphericalPoint(1., 1., 1.))
    scaled = GO.vertex_node(GO.UnitSphericalPoint(2., 2., 2.))
    different_z = GO.vertex_node(GO.UnitSphericalPoint(1., 1., -1.))
    @test GO._spherical_node_identity(vertex) == GO._spherical_node_identity(scaled)
    @test GO._spherical_node_identity(vertex) != GO._spherical_node_identity(different_z)
    tc = GO.TopologyComputer(GO.RelateMatrixPredicate(), rg, rg)
    for k in (vertex, scaled, different_z)
        tc.node_sections[k] = GO.NodeSections(k)
    end
    GO._merge_coincident_nodes!(tc)
    @test length(tc.node_sections) == 2
    @test haskey(tc.node_sections, different_z)
end

@testset "spherical polygon crossing a hole near an axis" begin
    outer = [(-35., -30.), (35., -30.), (35., 30.), (-35., 30.), (-35., -30.)]
    hole = [(-0.001, -0.001), (-0.001, 0.001), (0.001, 0.001),
        (0.001, -0.001), (-0.001, -0.001)]
    polygon = GI.Polygon([outer, hole])
    for sign in (-1., 1.), exact in (GO.True(), GO.False()),
            accelerator in (GO.NestedLoop(), GO.DoubleNaturalTree())
        # Four vertices of the HEALPix level-2 cell adjacent to the origin.
        # The small hole crosses both edges incident to the X-axis vertex.
        points = [
            GO.UnitSphericalPoint(0.9670673281592066, sign * 0.19236165165966837, 1/6),
            GO.UnitSphericalPoint(1., 0., 0.),
            GO.UnitSphericalPoint(0.9670673281592066, sign * 0.19236165165966837, -1/6),
            GO.UnitSphericalPoint(0.9238795325112867, sign * 0.3826834323650898, 0.),
        ]
        push!(points, first(points))
        cell = GI.Polygon([GI.LinearRing(points)])
        alg = GO.RelateNG(; manifold=GO.Spherical(), exact, accelerator)
        prepared = GO.prepare(alg, polygon)
        @test !GO.relate_predicate(prepared, GO.pred_contains(), cell)
        @test !GO.relate_predicate(prepared, GO.pred_covers(), cell)
        @test GO.relate_predicate(prepared, GO.pred_intersects(), cell)
        @test GO.relate_predicate(prepared, GO.pred_overlaps(), cell)
    end
end
