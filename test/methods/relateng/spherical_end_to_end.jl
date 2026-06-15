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
# (wrong DE-9IM); the bulge-aware `arc_extent` keeps it.
@testset "spherical tree accelerator agrees with NestedLoop (arc bulge)" begin
    A = GI.Polygon([GI.LinearRing([(0., 0.), (170., 0.), (85., 40.), (0., 0.)])])
    B = GI.Polygon([GI.LinearRing([(88., -2.), (92., -2.), (92., 2.), (88., 2.), (88., -2.)])])
    tree = RelateNG(; manifold = Spherical(), accelerator = GO.DoubleSTRtree())
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
# edge index is the dimension-generic `NaturalIndex` over 3D arc extents.
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

# Task 18: an exactly-antipodal edge has no unique great-circle arc; the kernel
# refuses it at ingest with a message pointing at the AntipodalEdgeSplit remedy.
@testset "spherical antipodal edge throws informatively" begin
    p0   = GO._to_kernel_point(Spherical(), (0., 0.))       # (1, 0, 0)
    p180 = GO._to_kernel_point(Spherical(), (180., 0.))     # (-1, 0, 0): antipodal
    p90  = GO._to_kernel_point(Spherical(), (90., 0.))      # (0, 1, 0)
    err = try; GO.arc_extent(p0, p180); nothing; catch e; e; end
    @test err isa ArgumentError
    @test occursin("AntipodalEdgeSplit", err.msg)
    #-- a normal edge and a repeated (zero-length) vertex do NOT throw
    @test (GO.arc_extent(p0, p90); true)
    @test (GO.arc_extent(p0, p0); true)
    #-- the whole relate rejects a polygon carrying an antipodal edge at ingest
    bad = GI.Polygon([GI.LinearRing([(0., 0.), (180., 0.), (90., 80.), (0., 0.)])])
    @test_throws ArgumentError GO.relate(alg, bad, GI.Point(10., 10.))
end
