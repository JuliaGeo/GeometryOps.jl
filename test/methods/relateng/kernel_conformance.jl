# RelateKernel conformance testset (design layer contract, Task 9).
#
# This suite is the *specification* a kernel implementation must satisfy:
# it is written as a function over a manifold so that a future `Spherical`
# kernel can be instantiated against the very same property checks. Inputs
# are drawn from a small integer grid with a fixed-seed RNG, so every
# orientation/classification has an exactly representable answer and the
# suite is deterministic across runs.

using Test
import GeometryOps as GO
import GeometryOps: Planar, True, False
import GeoInterface as GI
using Random

# `rk_orient` may return Float64 (plain determinant) or Int (exact
# predicate); all properties are stated on signs.
_sgn(x) = Int(sign(x))

# --- Ground-truth reference for the crossing-apex differential test ---
#
# A direct reimplementation of JTS PolygonNodeTopology.compareAngle,
# evaluated at the *exact rational apex* of a proper crossing using
# Rational{BigInt} arithmetic throughout (Float64 values are dyadic
# rationals, so the conversion is exact). This is intentionally independent
# of the kernel's symbolic-apex derivation in `rk_compare_edge_dir`.

# Quadrant convention identical to rk_quadrant: NE=0, NW=1, SW=2, SE=3,
# axis directions on the `>= 0` side.
function _ref_quadrant(ox::Rational{BigInt}, oy::Rational{BigInt}, p)
    px, py = Rational{BigInt}(p[1]), Rational{BigInt}(p[2])
    @assert !(px == ox && py == oy)
    if px >= ox
        return py >= oy ? 0 : 3
    else
        return py >= oy ? 1 : 2
    end
end

# Sign of orient(origin, q, p) over rationals (positive when p is CCW of q).
function _ref_orient_sign(ox::Rational{BigInt}, oy::Rational{BigInt}, q, p)
    R = Rational{BigInt}
    qx, qy = R(q[1]), R(q[2])
    px, py = R(p[1]), R(p[2])
    return Int(sign((qx - ox) * (py - oy) - (qy - oy) * (px - ox)))
end

# compareAngle anchored at the positive X-axis of the rational apex.
function _ref_compare_angle(apex::Tuple{Rational{BigInt}, Rational{BigInt}}, p, q)
    ox, oy = apex
    quadrant_p = _ref_quadrant(ox, oy, p)
    quadrant_q = _ref_quadrant(ox, oy, q)
    quadrant_p > quadrant_q && return 1
    quadrant_p < quadrant_q && return -1
    return _ref_orient_sign(ox, oy, q, p)
end

function kernel_conformance_suite(m; exact)
    # Fixed seed: same property sample within a Julia version (Julia only
    # guarantees within-version RNG stream reproducibility). The property
    # checks below are sample-independent, and the count thresholds
    # (`n_proper > 20`/`> 50`, etc.) are robust to any reasonable uniform
    # sample from this grid.
    rng = Random.MersenneTwister(0x5e1a7e)
    # Random point on a small integer grid: all predicates exactly decidable.
    rpt() = (Float64(rand(rng, -8:8)), Float64(rand(rng, -8:8)))

    @testset "rk_orient: antisymmetry / cyclic invariance / degeneracy" begin
        for _ in 1:500
            a, b, c = rpt(), rpt(), rpt()
            o = _sgn(GO.rk_orient(m, a, b, c; exact))
            # antisymmetry: swapping any two arguments flips the sign
            @test o == -_sgn(GO.rk_orient(m, b, a, c; exact))
            @test o == -_sgn(GO.rk_orient(m, a, c, b; exact))
            # cyclic invariance
            @test o == _sgn(GO.rk_orient(m, b, c, a; exact))
            @test o == _sgn(GO.rk_orient(m, c, a, b; exact))
            # degeneracy: any repeated point is collinear
            @test GO.rk_orient(m, a, a, b; exact) == 0
            @test GO.rk_orient(m, a, b, a; exact) == 0
            @test GO.rk_orient(m, a, b, b; exact) == 0
        end
    end

    @testset "rk_classify_intersection: symmetry and incidence consistency" begin
        n_proper = 0
        n_touch = 0
        n_collinear = 0
        for _ in 1:1000
            a0, a1, b0, b1 = rpt(), rpt(), rpt(), rpt()
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
            # SS_PROPER: crossing strictly interior to both segments,
            # so no endpoint incidence whatsoever
            if r.kind == GO.SS_PROPER
                n_proper += 1
                @test !(r.a0_on_b || r.a1_on_b || r.b0_on_a || r.b1_on_a)
                # a proper crossing means each endpoint is strictly off the
                # other segment's line: all four orientations are nonzero
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
        # the random sample must actually exercise the interesting kinds
        @test n_proper > 20
        @test n_touch > 0
        @test n_collinear > 0

        # shared-endpoint configurations (non-collinear) classify as touch
        r = GO.rk_classify_intersection(m, (0., 0.), (1., 0.), (0., 0.), (0., 1.); exact)
        @test r.kind == GO.SS_TOUCH && r.a0_on_b && r.b0_on_a
    end

    @testset "rk_compare_edge_dir: strict weak order on a 16-direction fan" begin
        # fan of 16 directions in strictly increasing angular order around
        # a non-origin apex, covering all four quadrants and both axes
        apex = (3.0, -2.0)
        node = GO.vertex_node(apex)
        dirs = [(1., 0.), (2., 1.), (1., 1.), (1., 2.), (0., 1.),     # NE
                (-1., 2.), (-1., 1.), (-2., 1.), (-1., 0.),           # NW
                (-2., -1.), (-1., -1.), (-1., -2.),                   # SW
                (0., -1.), (1., -2.), (1., -1.), (2., -1.)]           # SE
        fan = [(apex[1] + dx, apex[2] + dy) for (dx, dy) in dirs]
        cmp(p, q) = GO.rk_compare_edge_dir(m, node, p, q; exact)
        for p in fan
            @test cmp(p, p) == 0                       # irreflexivity of <
        end
        for p in fan, q in fan
            @test cmp(p, q) == -cmp(q, p)              # antisymmetry
        end
        for p in fan, q in fan, r in fan               # transitivity
            if cmp(p, q) < 0 && cmp(q, r) < 0
                @test cmp(p, r) < 0
            elseif cmp(p, q) == 0 && cmp(q, r) == 0
                @test cmp(p, r) == 0
            end
        end
        # the comparator reproduces the fan's angular order exactly
        for i in eachindex(fan), j in eachindex(fan)
            @test cmp(fan[i], fan[j]) == (i == j ? 0 : (i < j ? -1 : 1))
        end
    end

    # The crossing-apex differential test is planar-specific: it uses the
    # planar internal `GO._exact_crossing_point` and a planar-Cartesian
    # rational reference (`_ref_compare_angle`). A Spherical kernel must
    # supply its own differential reference for this property.
    if m isa GO.Planar
        @testset "rk_compare_edge_dir: crossing apex vs exact rational reference" begin
            # Differential test (Task 8 review item): for proper crossings on an
            # integer grid, the symbolic crossing-apex comparison must agree in
            # sign with compareAngle evaluated at the exact rational apex.
            n_proper = 0
            for _ in 1:2000
                a0, a1, b0, b1 = rpt(), rpt(), rpt(), rpt()
                r = GO.rk_classify_intersection(m, a0, a1, b0, b1; exact)
                r.kind == GO.SS_PROPER || continue
                n_proper += 1
                node = GO.crossing_node(a0, a1, b0, b1)
                apex = GO._exact_crossing_point(a0, a1, b0, b1)
                endpoints = (a0, a1, b0, b1)
                for p in endpoints, q in endpoints
                    @test _sgn(GO.rk_compare_edge_dir(m, node, p, q; exact)) ==
                          _ref_compare_angle(apex, p, q)
                end
            end
            @test n_proper > 50    # the filter kept a meaningful sample
        end
    end

    @testset "rk_nodes_coincide: reflexive, symmetric, consistent with ==" begin
        v1 = GO.vertex_node((2.0, 3.0))
        v2 = GO.vertex_node((2.0, 3.0))
        v3 = GO.vertex_node((2.0, 4.0))
        c1 = GO.crossing_node((0., 0.), (2., 2.), (0., 2.), (2., 0.))
        c2 = GO.crossing_node((2., 0.), (0., 2.), (2., 2.), (0., 0.))  # same pair, permuted
        c3 = GO.crossing_node((0., 0.), (3., 1.), (0., 1.), (3., 0.))  # crosses elsewhere
        keys = (v1, v2, v3, c1, c2, c3)
        for k in keys
            @test GO.rk_nodes_coincide(m, k, k; exact)            # reflexive
        end
        for k1 in keys, k2 in keys
            @test GO.rk_nodes_coincide(m, k1, k2; exact) ==
                  GO.rk_nodes_coincide(m, k2, k1; exact)          # symmetric
            if k1 == k2
                @test GO.rk_nodes_coincide(m, k1, k2; exact)      # == implies coincide
            end
        end
        @test v1 == v2 && GO.rk_nodes_coincide(m, v1, v2; exact)
        @test c1 == c2 && GO.rk_nodes_coincide(m, c1, c2; exact)
        @test !GO.rk_nodes_coincide(m, v1, v3; exact)
        @test !GO.rk_nodes_coincide(m, c1, c3; exact)
        # cross-kind coincidence: c1 crosses at exactly (1, 1)
        @test GO.rk_nodes_coincide(m, c1, GO.vertex_node((1.0, 1.0)); exact)
        @test !GO.rk_nodes_coincide(m, c3, GO.vertex_node((1.0, 1.0)); exact)
        # slow path: two *distinct* crossing pairs sharing the same apex —
        # keys differ, yet the nodes coincide at (1, 1)
        c_a = GO.crossing_node((0., 0.), (2., 2.), (0., 2.), (2., 0.))
        c_b = GO.crossing_node((1., 0.), (1., 2.), (0., 1.), (2., 1.))
        @test c_a != c_b
        @test GO.rk_nodes_coincide(m, c_a, c_b; exact)
    end

    @testset "rk_point_in_ring agrees with rk_point_on_segment on edges" begin
        # non-convex ring with a reflex vertex; integer coordinates so edge
        # midpoints are exactly representable
        ring_pts = [(0., 0.), (8., 0.), (8., 8.), (4., 4.), (0., 8.), (0., 0.)]
        ring = GI.LinearRing(ring_pts)
        edges = [(ring_pts[i], ring_pts[i + 1]) for i in 1:length(ring_pts)-1]
        on_edge(p) = any(GO.rk_point_on_segment(m, p, q0, q1; exact) for (q0, q1) in edges)
        # deterministic probes: vertices, edge midpoints, interior, exterior
        midpoints = [((q0[1] + q1[1]) / 2, (q0[2] + q1[2]) / 2) for (q0, q1) in edges]
        probes = vcat(ring_pts, midpoints, [(2.0, 2.0), (4.0, 5.0), (9.0, 1.0), (4.0, 7.0)])
        # randomized probes on the grid covering inside/outside/boundary
        for _ in 1:300
            push!(probes, (Float64(rand(rng, -1:9)), Float64(rand(rng, -1:9))))
        end
        for p in probes
            loc = GO.rk_point_in_ring(m, p, ring; exact)
            @test (loc == GO.LOC_BOUNDARY) == on_edge(p)
        end
    end
end

@testset "Kernel conformance: Planar" begin
    @testset "exact = $E" for E in (True(), False())
        # exact = False() is included because all suite inputs are small
        # integer-grid coordinates, for which plain Float64 determinants are
        # already exact — both code paths must satisfy the contract here.
        kernel_conformance_suite(Planar(); exact = E)
    end
end
