# Spherical edge-crossing conformance, ported from Google S2.
#
# Source: s2geometry `src/s2/s2edge_crosser_test.cc` (`TEST(S2, Crossings)`,
# `CollinearEdgesThatDontTouch`, `CoincidentZeroLengthEdgesThatDontTouch`) and
# `src/s2/s2edge_crossings.h` (the `CrossingSign` contract). These are S2's
# hand-built numerically-extreme edge-crossing cases — several require exact
# arithmetic to resolve (floating-point underflow; determinants needing >2000
# bits of precision) — and are the canonical battery for a robust spherical
# crossing predicate.
#
# Mapping S2 `CrossingSign(a,b,c,d)` onto GeometryOps `rk_classify_intersection`
# (validated empirically against the spherical kernel before this file landed):
#
#   CrossingSign == +1  (interior crossing of both edges) <-> SS_PROPER
#   CrossingSign == -1  (no crossing)                     <-> SS_DISJOINT
#   CrossingSign ==  0  (a vertex is shared)              <-> SS_TOUCH / SS_COLLINEAR
#                                                             with the incidence
#                                                             flag(s) set
#
# Two of S2's edge-crossing tests are deliberately NOT ported, because they have
# no faithful analog in this kernel:
#
#   * `GetIntersection` / `GrazingIntersections` compute an intersection *point*.
#     This kernel is symbolic: a proper crossing's node has no coordinate
#     anywhere in the engine (design D2), so there is nothing to compare.
#
#   * `CoincidentZeroLengthEdgesThatDontTouch` asserts that four points which are
#     *exactly proportional* on the sphere never cross. That holds only under
#     S2's symbolic-perturbation model of `Sign`. This kernel instead computes
#     the exact sign of the *actual* Float64 coordinates (via ExactPredicates /
#     Rational), and the rounding of `(1+k*eps)*p` makes the four points genuinely
#     non-collinear, so a tiny real interior crossing is the correct answer for
#     the given inputs. The S2 test's own comment notes it "depends on the
#     particular symbolic perturbations used by s2pred::Sign()". See the
#     `CoincidentProportionalEdges` testset below, which asserts the property
#     this kernel *does* guarantee (permutation-consistency) on such inputs.
#
# `RobustCrossProd` (also in s2edge_crossings_test.cc) is already ported, in
# test/utils/robustcrossproduct.jl.

using Test
import GeometryOps as GO
import GeometryOps: Spherical, True, False
import GeometryOps.UnitSpherical: UnitSphericalPoint, slerp
import GeoInterface as GI
using LinearAlgebra: normalize
using Random

# Build a unit-length spherical kernel point from raw components, normalizing
# exactly as S2's `S2Point::Normalize()` does (the S2 test harness normalizes
# every vertex before testing).
_nv(x, y, z) = (v = normalize([Float64(x), Float64(y), Float64(z)]);
                UnitSphericalPoint(v[1], v[2], v[3]))

# S2::Origin() — the conventional reference point, copied verbatim from
# s2pointutil.h. Used by S2 crossing cases 4 and 5 to exercise an edge whose
# endpoint is the origin reference.
const _S2_ORIGIN = (-0.0099994664350250197, 0.0025924542609324121, 0.99994664350250195)

# Port of S2's `TestCrossings(a, b, c, d, crossing_sign, ...)`. We assert the
# `CrossingSign`->`SegSegClass.kind` mapping plus the symmetry that S2's harness
# checks: `kind` is invariant under reversing either edge and under swapping the
# two edges.
function check_s2_crossing(m, a, b, c, d, crossing_sign; exact)
    r = GO.rk_classify_intersection(m, a, b, c, d; exact)

    if exact === True()
        # The exact path is authoritative — S2's `CrossingSign` is robust, so the
        # faithful analog is this kernel's exact path. Assert the full mapping.
        if crossing_sign == 1
            @test r.kind == GO.SS_PROPER
        elseif crossing_sign == -1
            @test r.kind == GO.SS_DISJOINT
        else # crossing_sign == 0: at least one vertex is shared between the edges
            @test r.kind != GO.SS_PROPER
            @test r.kind != GO.SS_DISJOINT
            @test r.a0_on_b || r.a1_on_b || r.b0_on_a || r.b1_on_a
        end
    else
        # The Float64 fast path is a conservative accelerator, not an oracle: it
        # falls back to SS_TOUCH on cases that need exact arithmetic (9, 11) and
        # cannot register an incidence whose coplanarity determinant does not
        # round to *exactly* zero (so a shared vertex of two skew normalized
        # edges, as in case 6, reads as SS_DISJOINT). The contract it must honor
        # is one-directional: never deny a real crossing, never invent one.
        if crossing_sign == 1
            @test r.kind != GO.SS_DISJOINT
        else
            @test r.kind != GO.SS_PROPER
        end
    end

    # S2 TestCrossings symmetry: classification is invariant under reversal of
    # either edge and under swapping the two edges. (Holds on both paths.)
    for (p, q, s, t) in ((b, a, c, d), (a, b, d, c), (b, a, d, c), (c, d, a, b))
        @test GO.rk_classify_intersection(m, p, q, s, t; exact).kind == r.kind
    end
    return r
end

function s2_crossings_suite(m; exact)
    @testset "S2 Crossings: 12 hand-built cases" begin
        # 1. Two regular edges that cross.
        check_s2_crossing(m, _nv(1, 2, 1), _nv(1, -3, 0.5),
                             _nv(1, -0.5, -3), _nv(0.1, 0.5, 3), 1; exact)
        # 2. Two regular edges that intersect at antipodal points (so the arcs
        #    themselves do not meet).
        check_s2_crossing(m, _nv(1, 2, 1), _nv(1, -3, 0.5),
                             _nv(-1, 0.5, 3), _nv(-0.1, -0.5, -3), -1; exact)
        # 3. Two edges on the same great circle that start at antipodal points.
        check_s2_crossing(m, _nv(0, 0, -1), _nv(0, 1, 0),
                             _nv(0, 0, 1), _nv(0, 1, 1), -1; exact)
        # 4. Two edges that cross where one vertex is S2::Origin().
        check_s2_crossing(m, _nv(1, 0, 0), _nv(_S2_ORIGIN...),
                             _nv(1, -0.1, 1), _nv(1, 1, -0.1), 1; exact)
        # 5. Antipodal intersection where one vertex is S2::Origin().
        check_s2_crossing(m, _nv(1, 0, 0), _nv(_S2_ORIGIN...),
                             _nv(-1, 0.1, -1), _nv(-1, -1, 0.1), -1; exact)
        # 6. Two edges that share an endpoint (2,3,4).
        check_s2_crossing(m, _nv(7, -2, 3), _nv(2, 3, 4),
                             _nv(2, 3, 4), _nv(-1, 2, 5), 0; exact)
        # 7. Edges that barely cross near the middle of one edge. AB is ~ in the
        #    x=y plane; CD is ~ perpendicular and ends exactly at the x=y plane.
        check_s2_crossing(m, _nv(1, 1, 1), _nv(1, prevfloat(1.0), -1),
                             _nv(11, -12, -1), _nv(10, 10, 1), 1; exact)
        # 8. As (7) but the edges are separated by a distance of about 1e-15.
        check_s2_crossing(m, _nv(1, 1, 1), _nv(1, nextfloat(1.0), -1),
                             _nv(1, -1, 0), _nv(1, 1, 0), -1; exact)
        # 9. Barely cross near the end of both edges — cannot be handled in plain
        #    double precision due to floating-point underflow.
        check_s2_crossing(m, _nv(0, 0, 1), _nv(2, -1e-323, 1),
                             _nv(1, -1, 1), _nv(1e-323, 0, 1), 1; exact)
        # 10. As (9) but separated by a distance of about 1e-640.
        check_s2_crossing(m, _nv(0, 0, 1), _nv(2, 1e-323, 1),
                             _nv(1, -1, 1), _nv(1e-323, 0, 1), -1; exact)
        # 11. Barely cross near the middle of one edge — the exact determinant of
        #     some triangles here needs more than 2000 bits of precision.
        check_s2_crossing(m, _nv(1, -1e-323, -1e-323), _nv(1e-323, 1, 1e-323),
                             _nv(1, -1, 1e-323), _nv(1, 1, 0), 1; exact)
        # 12. As (11) but separated by a distance of about 1e-640.
        check_s2_crossing(m, _nv(1, 1e-323, -1e-323), _nv(-1e-323, 1, 1e-323),
                             _nv(1, -1, 1e-323), _nv(1, 1, 0), -1; exact)
    end

    @testset "S2 CollinearEdgesThatDontTouch" begin
        # Two disjoint sub-arcs of one minor arc: a..b == arc[0.00,0.05] and
        # c..d == arc[0.95,1.00]. They share a great circle but never overlap, so
        # they must never be reported as a proper crossing.
        rng = MersenneTwister(0x5e2c011)
        rp() = (v = normalize(randn(rng, 3)); UnitSphericalPoint(v[1], v[2], v[3]))
        for _ in 1:500
            a = rp(); d = rp()
            b = slerp(a, d, 0.05)
            c = slerp(a, d, 0.95)
            r = GO.rk_classify_intersection(m, a, b, c, d; exact)
            @test r.kind != GO.SS_PROPER
            # On the exact path the kernel resolves them as fully disjoint.
            exact === True() && @test r.kind == GO.SS_DISJOINT
        end
    end

    @testset "CoincidentProportionalEdges (kernel convention, not S2's)" begin
        # S2's `CoincidentZeroLengthEdgesThatDontTouch` asserts non-crossing for
        # exactly-proportional points under symbolic perturbation. This kernel
        # uses exact signs of the actual Float64 coordinates instead, so it does
        # not promise that. What it *does* promise — and what we assert here — is
        # that the classification is self-consistent: invariant under reversing
        # either edge and under swapping the two edges (the same symmetry the S2
        # harness relies on).
        rng = MersenneTwister(0x5e2c022)
        for _ in 1:300
            p = normalize(randn(rng, 3))
            a = UnitSphericalPoint(((1 - 3e-16) .* p)...)
            b = UnitSphericalPoint(((1 - 1e-16) .* p)...)
            c = UnitSphericalPoint(p...)
            d = UnitSphericalPoint(((1 + 2e-16) .* p)...)
            r = GO.rk_classify_intersection(m, a, b, c, d; exact)
            for (w, x, y, z) in ((b, a, c, d), (a, b, d, c), (b, a, d, c), (c, d, a, b))
                @test GO.rk_classify_intersection(m, w, x, y, z; exact).kind == r.kind
            end
        end
    end
end

@testset "S2 edge-crossing conformance: Spherical" begin
    @testset "exact = $E" for E in (True(), False())
        s2_crossings_suite(Spherical(); exact = E)
    end
end
