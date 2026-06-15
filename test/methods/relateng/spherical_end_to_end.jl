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
