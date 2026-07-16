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
        # zero-length arc (repeated ring vertex): only its endpoint direction
        # is on it — the degenerate span test must not accept every point
        @test GO.rk_point_on_segment(m, a, a, a; exact)
        @test GO.rk_point_on_segment(m, _usp(2, 0, 0), a, a; exact)  # same direction, different scale
        @test !GO.rk_point_on_segment(m, b, a, a; exact)
        @test !GO.rk_point_on_segment(m, _usp(-1, 0, 0), a, a; exact) # antipode of the endpoint
    end
    # arc containment (bulge capture) is the shared `spherical_arc_extent`'s
    # contract, tested exhaustively in test/utils/unitspherical.jl

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
    @testset "angle ordering: vertex-apex fan reproduces planar order" begin
        # apex +x; reference axis is +y (argmin |component|), so the tangent
        # frame is (u,v) = (+y,+z) and a direction-point (0,dx,dy) has tangent
        # coords exactly (dx,dy) — the spherical comparator reproduces the planar
        # fan order verbatim.
        apex = _usp(1, 0, 0)
        node = GO.vertex_node(apex)
        dirs = [(1,0),(2,1),(1,1),(1,2),(0,1),
                (-1,2),(-1,1),(-2,1),(-1,0),
                (-2,-1),(-1,-1),(-1,-2),
                (0,-1),(1,-2),(1,-1),(2,-1)]
        fan = [_usp(0, dx, dy) for (dx, dy) in dirs]
        vcmp(p, q) = GO.rk_compare_edge_dir(m, node, p, q; exact)
        for p in fan
            @test vcmp(p, p) == 0
        end
        for i in eachindex(fan), j in eachindex(fan)
            @test vcmp(fan[i], fan[j]) == (i == j ? 0 : (i < j ? -1 : 1))
        end
        for i in eachindex(fan), j in eachindex(fan), k in eachindex(fan)
            if vcmp(fan[i], fan[j]) < 0 && vcmp(fan[j], fan[k]) < 0
                @test vcmp(fan[i], fan[k]) < 0
            end
        end
        @test_throws ArgumentError GO.rk_quadrant(m, apex, apex)  # zero-length direction
    end

    @testset "isCrossing / isInteriorSegment on a spherical corner" begin
        # apex +x, direction-point (0,dx,dy) -> planar direction (dx,dy); the
        # planar truth table maps over verbatim.
        nd = GO.vertex_node(_usp(1, 0, 0))
        d(dx, dy) = _usp(0, dx, dy)
        # X cross: a arms (1,1)/(-1,-1), b arms (-1,1)/(1,-1)
        @test GO.rk_is_crossing(m, nd, d(-1,-1), d(1,1), d(-1,1), d(1,-1); exact)
        # collinear b arm -> not crossing
        @test !GO.rk_is_crossing(m, nd, d(-1,-1), d(1,1), d(-1,-1), d(1,-1); exact)
        # both b arms same side of the a-line -> not crossing
        @test !GO.rk_is_crossing(m, nd, d(-1,-1), d(1,1), d(-1,1), d(0,1); exact)
        # crossing-node apex is rejected
        cn = GO.crossing_node(_usp(1,0,1), _usp(1,0,-1), _usp(1,1,0), _usp(1,-1,0))
        @test_throws ArgumentError GO.rk_is_crossing(m, cn, d(-1,-1), d(1,1), d(-1,1), d(1,-1); exact)

        # isInteriorSegment: corner a0->node->a1, ring interior on the right
        @test !GO.rk_is_interior_segment(m, nd, d(0,1), d(1,0), d(1,1); exact)
        @test GO.rk_is_interior_segment(m, nd, d(0,1), d(1,0), d(-1,-1); exact)
        @test GO.rk_is_interior_segment(m, nd, d(1,0), d(0,1), d(1,1); exact)
        @test !GO.rk_is_interior_segment(m, nd, d(1,0), d(0,1), d(-1,-1); exact)
    end

    @testset "rk_crossing_dirs_ccw and crossing-apex comparison" begin
        # proper crossing at +x=(1,0,0): a-arc (1,0,1)->(1,0,-1), b-arc (1,1,0)->(1,-1,0)
        a0, a1, b0, b1 = _usp(1,0,1), _usp(1,0,-1), _usp(1,1,0), _usp(1,-1,0)
        @test GO.rk_classify_intersection(m, a0, a1, b0, b1; exact).kind == GO.SS_PROPER
        ccw = GO.rk_crossing_dirs_ccw(m, a0, a1, b0, b1; exact)
        @test Set(ccw) == Set((a0, a1, b0, b1))   # a permutation of the four endpoints
        @test ccw[1] == a1                          # starts from a1 (contract)
        # the crossing-apex comparator is a strict weak order on the four arms
        cn = GO.crossing_node(a0, a1, b0, b1)
        ccmp(p, q) = GO.rk_compare_edge_dir(m, cn, p, q; exact)
        arms = (a0, a1, b0, b1)
        for x in arms
            @test ccmp(x, x) == 0
        end
        for x in arms, y in arms
            @test ccmp(x, y) == -ccmp(y, x)         # antisymmetry
            x === y || @test ccmp(x, y) != 0        # distinct arms have distinct angles
        end
    end
    @testset "rk_nodes_coincide: reflexive, symmetric, cross-kind" begin
        # crossing at +x via two different segment pairs
        cx_a = GO.crossing_node(_usp(1,0,1), _usp(1,0,-1), _usp(1,1,0), _usp(1,-1,0))
        cx_b = GO.crossing_node(_usp(2,1,0), _usp(2,-1,0), _usp(2,0,1), _usp(2,0,-1))
        # crossing at +y
        cy = GO.crossing_node(_usp(0,1,1), _usp(0,1,-1), _usp(1,1,0), _usp(-1,1,0))
        vx = GO.vertex_node(_usp(1,0,0))
        vx2 = GO.vertex_node(_usp(2,0,0))     # parallel to +x -> same sphere point
        vy = GO.vertex_node(_usp(0,1,0))
        keys = (cx_a, cx_b, cy, vx, vx2, vy)
        for k in keys
            @test GO.rk_nodes_coincide(m, k, k; exact)            # reflexive
        end
        for k1 in keys, k2 in keys
            @test GO.rk_nodes_coincide(m, k1, k2; exact) ==
                  GO.rk_nodes_coincide(m, k2, k1; exact)          # symmetric
        end
        # vertex coincidence by parallelism (not exact bit equality)
        @test vx != vx2 && GO.rk_nodes_coincide(m, vx, vx2; exact)
        @test !GO.rk_nodes_coincide(m, vx, vy; exact)
        # cross-kind: vertex lies on a proper crossing
        @test GO.rk_nodes_coincide(m, cx_a, vx; exact)
        @test GO.rk_nodes_coincide(m, cx_a, vx2; exact)
        @test !GO.rk_nodes_coincide(m, cy, vx; exact)
        # two distinct crossing pairs meeting at the same apex (+x)
        @test cx_a != cx_b && GO.rk_nodes_coincide(m, cx_a, cx_b; exact)
        @test !GO.rk_nodes_coincide(m, cx_a, cy; exact)
        # antipodal directions are NOT the same point
        @test !GO.rk_nodes_coincide(m, GO.vertex_node(_usp(1,1,0)), GO.vertex_node(_usp(-1,-1,0)); exact)
    end
    @testset "rk_point_in_ring: parity, boundary, winding independence" begin
        # CCW-from-above diamond at z=1 encircling the north pole; integer
        # vertices (membership decidable: boundary via exact coplanarity, the
        # meridian-parity crossings are scale-invariant signs).
        verts = [_usp(2,0,1), _usp(0,2,1), _usp(-2,0,1), _usp(0,-2,1), _usp(2,0,1)]
        ring = GI.LinearRing(verts)
        @test GO.rk_point_in_ring(m, _usp(2,0,1), ring; exact) == GO.LOC_BOUNDARY     # vertex
        @test GO.rk_point_in_ring(m, _usp(1,1,1), ring; exact) == GO.LOC_BOUNDARY     # edge midpoint
        @test GO.rk_point_in_ring(m, _usp(-1,1,1), ring; exact) == GO.LOC_BOUNDARY    # edge midpoint
        @test GO.rk_point_in_ring(m, _usp(1,1,10), ring; exact) == GO.LOC_INTERIOR    # polar cap
        @test GO.rk_point_in_ring(m, _usp(1,0,8), ring; exact) == GO.LOC_INTERIOR
        @test GO.rk_point_in_ring(m, _usp(3,1,0), ring; exact) == GO.LOC_EXTERIOR     # equatorial
        @test GO.rk_point_in_ring(m, _usp(3,-1,0), ring; exact) == GO.LOC_EXTERIOR

        # Winding independence (kernel contract): a ring encloses the same
        # region — the one smaller than a hemisphere — in either winding, like
        # the planar ray-crossing parity. Boundary is unchanged too.
        rev = GI.LinearRing(reverse(verts))
        @test GO.rk_point_in_ring(m, _usp(1,1,10), rev; exact) == GO.LOC_INTERIOR
        @test GO.rk_point_in_ring(m, _usp(3,1,0), rev; exact) == GO.LOC_EXTERIOR
        @test GO.rk_point_in_ring(m, _usp(1,1,1), rev; exact) == GO.LOC_BOUNDARY
    end
    @testset "area interaction bounds reach the enclosed pole" begin
        verts = [_usp(2,0,1), _usp(0,2,1), _usp(-2,0,1), _usp(0,-2,1), _usp(2,0,1)]
        poly = GI.Polygon([GI.LinearRing(verts)])
        e = GO.rk_interaction_bounds(m, poly)
        @test e.Z[2] >= 1.0                # interior reaches the enclosed north pole
        # the boundary ring (curve) bound tops out well below the pole
        eb = GO.rk_interaction_bounds(m, GI.LinearRing(verts))
        @test eb.Z[2] < 0.99
        # only the enclosed +z axis is extended; the equatorial axes are not
        @test e.X == eb.X && e.Y == eb.Y && e.Z[1] == eb.Z[1]
    end
    # testsets added task-by-task below
end

@testset "Kernel conformance: Spherical" begin
    @testset "exact = $E" for E in (True(), False())
        kernel_conformance_suite_spherical(Spherical(); exact = E)
    end
end

# `_ring_is_ccw(::Spherical)`: the S2 turning-angle (`GetCurvature`)
# formulation. Exact-independent (turn signs are always exact), so tested
# outside the exact-parameterized suite.
@testset "spherical ring orientation: turning-angle curvature" begin
    m = Spherical()
    curvature(pts) = GO._spherical_loop_curvature(
        GO._prune_loop_degeneracies([GO.rk_normalize_usp(p) for p in pts]))

    @testset "antipodal vertex pair (AntipodalEdgeSplit output)" begin
        # the 80° lune [(0,0), (90,0), (180,0), (90,80)]: its first and third
        # vertices are antipodal — the fan formulation this replaces NaN'd on
        # the chord between them; every turn angle here is between adjacent
        # vertices, so the pair never meets in one term
        lune = [GO._to_kernel_point(m, p) for p in [(0., 0.), (90., 0.), (180., 0.), (90., 80.)]]
        @test curvature(lune) ≈ deg2rad(200) atol = 1e-12       # area 2π − curvature ≈ 2.7925268 sr
        @test GO._ring_is_ccw(m, lune; exact = True())
        @test !GO._ring_is_ccw(m, reverse(lune); exact = True())
    end

    @testset "exactly negated under reversal, invariant under rotation" begin
        rng = Random.Xoshiro(0xc0ffee)
        for _ in 1:50
            k = rand(rng, 3:9)
            ring = [GO.rk_normalize_usp(_usp(randn(rng), randn(rng), randn(rng))) for _ in 1:k]
            c = curvature(ring)
            @test curvature(reverse(ring)) == -c                # exact, not ≈
            for r in 1:(k - 1)
                @test curvature(vcat(ring[(r + 1):end], ring[1:r])) == c
            end
        end
    end

    @testset "exact hemisphere is normalized in both windings" begin
        # equator ring: curvature 0, intrinsically winding-ambiguous — the
        # `>= -maxError` rule (S2 IsNormalized) reads it CCW in both windings
        eq = [_usp(1, 0, 0), _usp(0, 1, 0), _usp(-1, 0, 0), _usp(0, -1, 0)]
        @test GO._ring_is_ccw(m, eq; exact = True())
        @test GO._ring_is_ccw(m, reverse(eq); exact = True())
    end

    @testset "degeneracy pruning" begin
        cap = [_usp(2, 0, 1), _usp(0, 2, 1), _usp(-2, 0, 1), _usp(0, -2, 1)]
        # duplicate run (AA) plus a retraced whisker (ABA) leave the bit alone
        spiked = [cap[1], cap[1], cap[2], _usp(5, 5, 9), cap[2], cap[3], cap[4]]
        @test GO._ring_is_ccw(m, spiked; exact = True()) == GO._ring_is_ccw(m, cap; exact = True()) == true
        @test curvature(spiked) == curvature(cap)
        # whisker straddling the closure: [B, A, X, Y, A] prunes to [A, X, Y]
        wrap = [cap[2], cap[1], _usp(3, 1, 2), _usp(1, 3, 2), cap[1]]
        @test length(GO._prune_loop_degeneracies(wrap)) == 3
        # a ring degenerate after pruning bounds no area
        @test !GO._ring_is_ccw(m, [cap[1], cap[2], cap[1], cap[2]]; exact = True())
        @test !GO._ring_is_ccw(m, [cap[1], cap[2]]; exact = True())
    end

    @testset "closed and open forms agree" begin
        cap = [_usp(2, 0, 1), _usp(0, 2, 1), _usp(-2, 0, 1), _usp(0, -2, 1)]
        closed = vcat(cap, [cap[1]])
        @test GO._ring_is_ccw(m, closed; exact = True()) == GO._ring_is_ccw(m, cap; exact = True())
        @test !GO._ring_is_ccw(m, reverse(closed); exact = True())
    end
end
