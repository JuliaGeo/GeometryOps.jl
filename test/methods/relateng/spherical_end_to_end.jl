# End-to-end `relate(RelateNG(; manifold = Spherical()), …)` smoke tests
# (Task 14). With the spherical kernel and the manifold-derived point type
# in place, the engine runs on the unindexed NestedLoop path. Away from the
# poles and the antimeridian, spherical and planar topology agree (the spike
# measured 2000/2000), so the patch test compares the two DE-9IM matrices
# directly; the pole-containing ring is a case only the spherical kernel
# gets right.

using Test
import GeometryOps as GO
import GeometryOps: Spherical, RelateNG
import GeoInterface as GI

alg = RelateNG(; manifold = Spherical())

@testset "spherical relate agrees with planar on a small mid-latitude patch" begin
    # two overlapping boxes near (10°, 45°): away from poles/antimeridian the
    # spherical and planar topology agree
    A = GI.Polygon([GI.LinearRing([(0., 40.), (10., 40.), (10., 50.), (0., 50.), (0., 40.)])])
    B = GI.Polygon([GI.LinearRing([(5., 45.), (15., 45.), (15., 55.), (5., 55.), (5., 45.)])])
    @test GO.relate(alg, A, B) == GO.relate(A, B)           # same DE-9IM
    @test GO.relate(alg, A, B, "T*T***T**")                  # overlaps
end

@testset "spherical handles a pole-containing ring" begin
    cap = GI.Polygon([GI.LinearRing([(0., 80.), (120., 80.), (240., 80.), (0., 80.)])])
    pt  = GI.Point(0., 89.)                                  # near north pole
    @test GO.relate(alg, cap, pt, "T*****FF*")               # contains
end

# Task 15: the tree accelerator must build 3D great-circle arc extents, not a
# 2D coordinate box. A long near-equatorial arc (lon 0 → 170) bulges to y ≈ 1
# at lon 90 while its endpoint box has y ∈ [0, 0.17]; a polygon straddling the
# equator at lon 90 crosses it there. A 2D endpoint box prunes that pair away
# (wrong DE-9IM); the bulge-aware `spherical_arc_extent` keeps it.
@testset "spherical tree accelerator agrees with NestedLoop (arc bulge)" begin
    A = GI.Polygon([GI.LinearRing([(0., 0.), (170., 0.), (85., 40.), (0., 0.)])])
    B = GI.Polygon([GI.LinearRing([(88., -2.), (92., -2.), (92., 2.), (88., 2.), (88., -2.)])])
    tree = RelateNG(; manifold = Spherical(), accelerator = GO.DoubleNaturalTree())
    loop = RelateNG(; manifold = Spherical(), accelerator = GO.NestedLoop())
    @test GO.relate(tree, A, B) == GO.relate(loop, A, B)
end

# Task 15: above the 32-segment threshold AutoAccelerator picks the tree path;
# it must give the same DE-9IM as the unindexed nested loop on a real ring.
@testset "spherical AutoAccelerator picks the tree above threshold" begin
    n = 48                                                   # 48 segments > threshold
    ringpts = [(10.0 + 8cosd(t), 45.0 + 5sind(t)) for t in range(0, 360; length = n + 1)]
    A = GI.Polygon([GI.LinearRing(ringpts)])
    B = GI.Polygon([GI.LinearRing([(8., 43.), (20., 43.), (20., 52.), (8., 52.), (8., 43.)])])
    loop = RelateNG(; manifold = Spherical(), accelerator = GO.NestedLoop())
    @test GO.relate(alg, A, B) == GO.relate(loop, A, B)
end

# Task 17: prepared spherical relate (A indexed once, in 3D) must agree with
# the unprepared nested-loop relate over several B geometries. The prepared
# edge index is the dimension-generic natural-order `RTree` over 3D arc extents.
@testset "spherical prepared relate agrees with unprepared" begin
    n = 48
    ringpts = [(10.0 + 8cosd(t), 45.0 + 5sind(t)) for t in range(0, 360; length = n + 1)]
    A = GI.Polygon([GI.LinearRing(ringpts)])
    prep = GO.prepare(alg, A)
    loop = RelateNG(; manifold = Spherical(), accelerator = GO.NestedLoop())
    Bs = (
        GI.Polygon([GI.LinearRing([(8., 43.), (20., 43.), (20., 52.), (8., 52.), (8., 43.)])]),
        GI.Polygon([GI.LinearRing([(11., 45.), (13., 45.), (13., 47.), (11., 47.), (11., 45.)])]),
        GI.Point(10., 45.),
        GI.Point(40., 80.),
    )
    for B in Bs
        @test GO.relate(prep, B) == GO.relate(loop, A, B)
    end
end

# Ring winding must not change which region a polygon bounds: real-world data
# arrives in either winding (Natural Earth ships shapefile-convention CW
# shells), and the kernel contract for `rk_point_in_ring` is the area enclosed
# by the ring — the region smaller than a hemisphere — like the planar
# ray-crossing parity. Under the previous S2 left-of-ring reading, a CW
# continent meant "the whole sphere minus the continent" and intersected
# everything.
@testset "polygon interior is winding-independent" begin
    sq(x0, y0, s) = [(x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s), (x0, y0)]
    ccw(pts) = GI.Polygon([GI.LinearRing(pts)])
    cw(pts) = GI.Polygon([GI.LinearRing(reverse(pts))])
    canada_ish = sq(-100.0, 50.0, 10.0)
    australia_ish = sq(140.0, -35.0, 15.0)
    continent = sq(-140.0, 20.0, 60.0)                       # continent-sized
    far_pt = GI.Point(150.0, -25.0)
    @testset "disjoint stays disjoint in all four winding combinations" begin
        for A in (ccw(canada_ish), cw(canada_ish)), B in (ccw(australia_ish), cw(australia_ish))
            @test !GO.intersects(alg, A, B)
            @test GO.disjoint(alg, A, B)
        end
        for A in (ccw(continent), cw(continent))
            @test !GO.intersects(alg, A, far_pt)
            @test !GO.intersects(alg, A, ccw(australia_ish))
        end
    end
    @testset "containment in all four winding combinations" begin
        inner = sq(-97.0, 53.0, 2.0)
        for A in (ccw(canada_ish), cw(canada_ish)), B in (ccw(inner), cw(inner))
            @test GO.contains(alg, A, B)
            @test GO.within(alg, B, A)
            @test !GO.touches(alg, A, B)
        end
    end
    @testset "shared-edge neighbors touch in all four winding combinations" begin
        left = sq(0.0, 40.0, 10.0)
        right = sq(10.0, 40.0, 10.0)
        for A in (ccw(left), cw(left)), B in (ccw(right), cw(right))
            @test GO.touches(alg, A, B)
            @test GO.intersects(alg, A, B)
        end
    end
    @testset "equator-symmetric ring (planar isCCW tie case)" begin
        # two vertices tie for the extreme xyz-y coordinate, which broke the
        # planar extreme-vertex orientation test on the sphere: both windings
        # read CW
        eq = sq(-5.0, -5.0, 10.0)
        inside = GI.Point(0.0, 0.0)
        for A in (ccw(eq), cw(eq))
            @test GO.contains(alg, A, inside)
            @test !GO.intersects(alg, A, far_pt)
        end
    end
end

# A degenerate zero-area sliver ring with a repeated vertex (NE 110m North
# Korea ships an `[A, A, B, A]` first polygon) must not swallow the sphere:
# the zero-length arc `A → A` used to read every query point as on-boundary.
@testset "degenerate sliver polygon in a multipolygon stays local" begin
    a = (130.78, 42.22)
    sliver = GI.LinearRing([a, a, (130.78002, 42.22), a])
    body = GI.LinearRing([(124., 38.), (130., 38.), (130., 43.), (124., 43.), (124., 38.)])
    mp = GI.MultiPolygon([GI.Polygon([sliver]), GI.Polygon([body])])
    far = GI.Polygon([GI.LinearRing([(-60., -25.), (-55., -25.), (-55., -20.), (-60., -20.), (-60., -25.)])])
    @test !GO.intersects(alg, mp, far)
    @test !GO.intersects(alg, mp, GI.Point(-60.0, -20.0))
    @test GO.intersects(alg, mp, GI.Point(127.0, 40.0))
end

# The winding-independent region also drives the kernel interaction bounds: a
# CW polar cap must still bound the cap (its enclosed pole included), not the
# complement — which would under-cover the pole and prune away a contained
# point.
@testset "CW polar cap still means the cap" begin
    cap_cw = GI.Polygon([GI.LinearRing([(0., 80.), (240., 80.), (120., 80.), (0., 80.)])])
    pt = GI.Point(0., 89.)                                   # near north pole
    @test GO.relate(alg, cap_cw, pt, "T*****FF*")            # contains
    e = GO.rk_interaction_bounds(Spherical(), cap_cw)
    @test e.Z[2] >= 1.0                                      # bounds reach the enclosed pole
end

# Prepared point location goes through the indexed spherical locator (the
# longitude-interval edge index with meridian-arc crossing parity against a
# pole anchor); unprepared location is the exact ring scan. End to end the
# two must produce the same DE-9IM for every point, including the awkward
# geometries: an antimeridian-straddling box (split index intervals), polar
# caps over either pole (enclosed pole, polar-axis queries), and a shell
# with the south pole on its boundary (the parity anchor falls back to the
# north pole). The grid deliberately includes both poles, the ±180° seam,
# and points on ring boundaries.
@testset "prepared point location agrees with unprepared" begin
    rings = (
        [(170., -10.), (-170., -10.), (-170., 10.), (170., 10.), (170., -10.)],  # antimeridian box
        [(0., 80.), (120., 80.), (240., 80.), (0., 80.)],                        # north polar cap
        [(0., -80.), (120., -80.), (240., -80.), (0., -80.)],                    # south polar cap
        [(0., -90.), (20., -60.), (-20., -60.), (0., -90.)],                     # south pole on boundary
    )
    qs = [GI.Point(lon, lat) for lon in -180.0:30.0:180.0 for lat in -90.0:15.0:90.0]
    push!(qs, GI.Point(180.0, 5.0), GI.Point(175.0, 0.0), GI.Point(0.0, -70.0))
    for pts in rings, w in (pts, reverse(pts))
        A = GI.Polygon([GI.LinearRing(w)])
        prep = GO.prepare(alg, A)
        n_mismatch = count(q -> GO.relate(prep, q) != GO.relate(alg, A, q), qs)
        @test n_mismatch == 0
    end
end

# A ring may carry antipodal VERTEX pairs even though antipodal edges are
# rejected at ingest — the AntipodalEdgeSplit output is exactly such a ring.
# The Girard orientation fan must not run a chord through a vertex antipodal
# to its apex (the fan triangle has no defined geodesic and its excess
# degenerates to zero, flipping the shared orientation bit and with it the
# polygon interior — CI regression on the fixed ring below).
@testset "antipodal vertex pair does not flip the interior" begin
    # AntipodalEdgeSplit's output for the (0,0) → (180,0) edge, verbatim:
    # the first vertex (0,0) is antipodal to the mid-ring vertex (180,0)
    split_pts = [(0., 0.), (90., 0.), (180., 0.), (90., 80.), (0., 0.)]
    inside = GI.Point(10., 10.)      # under the (0,0)→(90,80) arc (lat 44.6° at lon 10)
    outside = GI.Point(10., -10.)    # southern hemisphere
    for pts in (split_pts, reverse(split_pts))
        poly = GI.Polygon([GI.LinearRing(pts)])
        @test GO.relate(alg, poly, inside, "T*****FF*")      # contains
        @test GO.contains(alg, poly, inside)
        @test !GO.intersects(alg, poly, outside)
    end
    # a different configuration: the first vertex's antipode mid-ring, with
    # the ring elsewhere entirely off the first vertex's great circles
    pts2 = [(20., 30.), (100., 10.), (-160., -30.), (-100., 10.), (20., 30.)]
    inside2 = GI.Point(0., 40.)
    for pts in (pts2, reverse(pts2))
        poly = GI.Polygon([GI.LinearRing(pts)])
        @test GO.contains(alg, poly, inside2)
        @test !GO.intersects(alg, poly, GI.Point(0., -80.))
    end
end

# Task 18: an exactly-antipodal edge has no unique great-circle arc; the kernel
# refuses it at ingest with a message pointing at the AntipodalEdgeSplit remedy.
@testset "spherical antipodal edge throws informatively" begin
    p0   = GO._to_kernel_point(Spherical(), (0., 0.))       # (1, 0, 0)
    p180 = GO._to_kernel_point(Spherical(), (180., 0.))     # (-1, 0, 0): antipodal
    p90  = GO._to_kernel_point(Spherical(), (90., 0.))      # (0, 1, 0)
    err = try; GO._validate_relate_edges(Spherical(), GI.LineString([p0, p180])); nothing; catch e; e; end
    @test err isa ArgumentError
    @test occursin("AntipodalEdgeSplit", err.msg)
    #-- a normal edge and a repeated (zero-length) vertex do NOT throw
    @test (GO._validate_relate_edges(Spherical(), GI.LineString([p0, p90])); true)
    @test (GO._validate_relate_edges(Spherical(), GI.LineString([p0, p0])); true)
    #-- the whole relate rejects a polygon carrying an antipodal edge at ingest
    bad = GI.Polygon([GI.LinearRing([(0., 0.), (180., 0.), (90., 80.), (0., 0.)])])
    @test_throws ArgumentError GO.relate(alg, bad, GI.Point(10., 10.))
end
