# Spherical RelateKernel conformance suite (design layer contract, Task 9).
# Standalone counterpart to kernel_conformance.jl. Inputs are exactly-
# representable integer xyz vectors in *general position*: every predicate is
# a sign of an integer determinant, so the exact answer is unambiguous and the
# suite is deterministic. Points are NOT unit length — sign predicates do not
# require it; where unit length matters (arc membership) we use exact-on-grid
# configurations.

using Test
import GeometryOps as GO
import GeometryOps: Spherical, True, False
import GeometryOps.UnitSpherical: UnitSphericalPoint, slerp
import GeometryOps.Extents as Extents
import GeoInterface as GI
using Random
using LinearAlgebra: ⋅, cross

const USPt = UnitSphericalPoint{Float64}
_sgn(x) = Int(sign(x))
_usp(x, y, z) = UnitSphericalPoint(Float64(x), Float64(y), Float64(z))

function kernel_conformance_suite_spherical(m; exact)
    rng = Random.MersenneTwister(0x5e1a7e)
    # random integer-component direction (a point on the sphere only up to
    # scaling; sign predicates are scale-invariant in each argument)
    rpt() = _usp(rand(rng, -8:8), rand(rng, -8:8), rand(rng, -8:8))
    function nonzero()
        p = rpt()
        while iszero(GI.x(p)) && iszero(GI.y(p)) && iszero(GI.z(p))
            p = rpt()
        end
        return p
    end

    @testset "rk_orient: antisymmetry / cyclic invariance / degeneracy" begin
        for _ in 1:500
            a, b, c = nonzero(), nonzero(), nonzero()
            o = _sgn(GO.rk_orient(m, a, b, c; exact))
            @test o == -_sgn(GO.rk_orient(m, b, a, c; exact))
            @test o == -_sgn(GO.rk_orient(m, a, c, b; exact))
            @test o == _sgn(GO.rk_orient(m, b, c, a; exact))
            @test o == _sgn(GO.rk_orient(m, c, a, b; exact))
            @test GO.rk_orient(m, a, a, b; exact) == 0
            @test GO.rk_orient(m, a, b, a; exact) == 0
            @test GO.rk_orient(m, a, b, b; exact) == 0
        end
    end
    @testset "rk_point_on_segment: endpoints, midpoint, off-arc" begin
        a = _usp(1, 0, 0); b = _usp(0, 1, 0)
        @test GO.rk_point_on_segment(m, a, a, b; exact)              # endpoint
        @test GO.rk_point_on_segment(m, b, a, b; exact)              # endpoint
        @test GO.rk_point_on_segment(m, _usp(1, 1, 0), a, b; exact)  # interior (same great circle, within span)
        @test !GO.rk_point_on_segment(m, _usp(0, 0, 1), a, b; exact) # pole, off the circle
        @test !GO.rk_point_on_segment(m, _usp(-1, 1, 0), a, b; exact) # on circle, outside the minor-arc span
    end
    @testset "arc_extent contains the arc (bulge captured)" begin
        rng2 = Random.MersenneTwister(42)
        for _ in 1:500
            p = GO.rk_normalize_usp(_usp(randn(rng2), randn(rng2), randn(rng2)))
            q = GO.rk_normalize_usp(_usp(randn(rng2), randn(rng2), randn(rng2)))
            e = GO.arc_extent(p, q)
            for s in 0:0.05:1
                u = slerp(p, q, s)
                @test e.X[1] <= GI.x(u) <= e.X[2]
                @test e.Y[1] <= GI.y(u) <= e.Y[2]
                @test e.Z[1] <= GI.z(u) <= e.Z[2]
            end
        end
    end

    @testset "rk_interaction_bounds is 3D and contains the converted vertices" begin
        ring = GI.LinearRing([(0.,0.), (10.,0.), (10.,10.), (0.,10.), (0.,0.)])
        e = GO.rk_interaction_bounds(m, ring)
        @test hasproperty(e, :Z)
        for p in GI.getpoint(ring)
            u = GO.rk_normalize_usp(UnitSphericalPoint((Float64(GI.x(p)), Float64(GI.y(p)))))
            @test e.X[1] <= GI.x(u) <= e.X[2]
            @test e.Y[1] <= GI.y(u) <= e.Y[2]
            @test e.Z[1] <= GI.z(u) <= e.Z[2]
        end
    end
    @testset "rk_bounds_covers respects Z" begin
        big = Extents.Extent(X = (0., 2.), Y = (0., 2.), Z = (0., 2.))
        inside = Extents.Extent(X = (0.5, 1.), Y = (0.5, 1.), Z = (0.5, 1.))
        outsideZ = Extents.Extent(X = (0.5, 1.), Y = (0.5, 1.), Z = (0.5, 3.))
        @test GO.rk_bounds_covers(big, inside)
        @test !GO.rk_bounds_covers(big, outsideZ)
        @test !GO.rk_bounds_disjoint(big, inside)
        # the 2D covering relation is unchanged
        big2 = Extents.Extent(X = (0., 2.), Y = (0., 2.))
        @test GO.rk_bounds_covers(big2, Extents.Extent(X = (0.5, 1.), Y = (0.5, 1.)))
        @test !GO.rk_bounds_covers(big2, Extents.Extent(X = (0.5, 3.), Y = (0.5, 1.)))
    end
    @testset "rk_classify_intersection: symmetry and incidence consistency" begin
        n_proper = 0; n_touch = 0; n_collinear = 0
        for _ in 1:2000
            a0, a1, b0, b1 = nonzero(), nonzero(), nonzero(), nonzero()
            r = GO.rk_classify_intersection(m, a0, a1, b0, b1; exact)
            # swapping A and B: kind invariant, flag pairs permuted
            s = GO.rk_classify_intersection(m, b0, b1, a0, a1; exact)
            @test s.kind == r.kind
            @test (s.a0_on_b, s.a1_on_b, s.b0_on_a, s.b1_on_a) ==
                  (r.b0_on_a, r.b1_on_a, r.a0_on_b, r.a1_on_b)
            # reversing a segment: kind invariant, its two flags swapped
            v = GO.rk_classify_intersection(m, a1, a0, b0, b1; exact)
            @test v.kind == r.kind
            @test (v.a0_on_b, v.a1_on_b, v.b0_on_a, v.b1_on_a) ==
                  (r.a1_on_b, r.a0_on_b, r.b0_on_a, r.b1_on_a)
            w = GO.rk_classify_intersection(m, a0, a1, b1, b0; exact)
            @test w.kind == r.kind
            @test (w.a0_on_b, w.a1_on_b, w.b0_on_a, w.b1_on_a) ==
                  (r.a0_on_b, r.a1_on_b, r.b1_on_a, r.b0_on_a)
            if r.kind == GO.SS_PROPER
                n_proper += 1
                @test !(r.a0_on_b || r.a1_on_b || r.b0_on_a || r.b1_on_a)
                # proper: each endpoint strictly off the other arc's great circle
                @test GO.rk_orient(m, b0, b1, a0; exact) != 0
                @test GO.rk_orient(m, b0, b1, a1; exact) != 0
                @test GO.rk_orient(m, a0, a1, b0; exact) != 0
                @test GO.rk_orient(m, a0, a1, b1; exact) != 0
            end
            r.kind == GO.SS_TOUCH && (n_touch += 1)
            r.kind == GO.SS_COLLINEAR && (n_collinear += 1)
            # each incidence flag agrees exactly with rk_point_on_segment
            @test r.a0_on_b == GO.rk_point_on_segment(m, a0, b0, b1; exact)
            @test r.a1_on_b == GO.rk_point_on_segment(m, a1, b0, b1; exact)
            @test r.b0_on_a == GO.rk_point_on_segment(m, b0, a0, a1; exact)
            @test r.b1_on_a == GO.rk_point_on_segment(m, b1, a0, a1; exact)
        end
        # proper crossings are common for two random great circles; touch and
        # collinear (two coplanar great circles) are rare/measure-zero on the
        # sphere, so they are exercised by the explicit cases below.
        @test n_proper > 20

        # --- hand-built decidable configurations ---
        # proper crossing: two axis-aligned great circles meeting at +x=(1,0,0),
        # strictly interior to both minor arcs
        r = GO.rk_classify_intersection(m, _usp(1,0,1), _usp(1,0,-1), _usp(1,1,0), _usp(1,-1,0); exact)
        @test r.kind == GO.SS_PROPER
        @test !(r.a0_on_b || r.a1_on_b || r.b0_on_a || r.b1_on_a)
        # shared-endpoint touch (non-collinear arcs sharing (1,0,0))
        r = GO.rk_classify_intersection(m, _usp(1,0,0), _usp(0,1,0), _usp(1,0,0), _usp(0,0,1); exact)
        @test r.kind == GO.SS_TOUCH && r.a0_on_b && r.b0_on_a
        # collinear overlap on the equator: arcs [+x,(1,1,0)] and [(1,1,0)... ] overlapping
        r = GO.rk_classify_intersection(m, _usp(1,0,0), _usp(0,1,0), _usp(1,1,0), _usp(-1,1,0); exact)
        @test r.kind == GO.SS_COLLINEAR
        # collinear disjoint on the equator: [+x, (2,1,0)] vs [(-1,2,0), -y]
        r = GO.rk_classify_intersection(m, _usp(1,0,0), _usp(2,1,0), _usp(-1,2,0), _usp(0,-1,0); exact)
        @test r.kind == GO.SS_DISJOINT
        # T-touch: b0 on the interior of arc a (a on equator, b dips to the pole)
        r = GO.rk_classify_intersection(m, _usp(1,0,0), _usp(0,1,0), _usp(1,1,0), _usp(0,0,1); exact)
        @test r.kind == GO.SS_TOUCH && r.b0_on_a && !r.a0_on_b
    end
    # testsets added task-by-task below
end

@testset "Kernel conformance: Spherical" begin
    @testset "exact = $E" for E in (True(), False())
        kernel_conformance_suite_spherical(Spherical(); exact = E)
    end
end
