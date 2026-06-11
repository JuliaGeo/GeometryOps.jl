# Tests for the planar RelateKernel: orientation, point-on-segment,
# point-in-ring and interaction bounds (design doc D1/D2).

using Test
import GeometryOps as GO
import GeometryOps: Planar, True, False
import GeoInterface as GI

const PT = Tuple{Float64, Float64}
m = Planar()

@testset "rk_orient" begin
    @test GO.rk_orient(m, (0.0,0.0), (1.0,0.0), (0.0,1.0); exact = True()) > 0
    @test GO.rk_orient(m, (0.0,0.0), (1.0,0.0), (0.0,-1.0); exact = True()) < 0
    @test GO.rk_orient(m, (0.0,0.0), (1.0,0.0), (2.0,0.0); exact = True()) == 0
    # adversarial near-collinear: exact must get the sign right.
    # NOTE: the plan used `0.5 + 1e-17`, but that rounds to exactly 0.5 in
    # Float64 (eps(0.5)/2 ≈ 5.6e-17), making the point exactly collinear.
    # Use nextfloat(0.5) — one ulp above the line — so the test is meaningful.
    a, b = (0.0, 0.0), (1.0, 1.0)
    c = (0.5, nextfloat(0.5))   # above the line by one ulp
    @test c[2] != 0.5           # guard: perturbation survives rounding
    @test GO.rk_orient(m, a, b, c; exact = True()) == GO.rk_orient(m, a, b, (0.5, 0.6); exact = True())
    # plain FP determinant evaluates to exactly 0.0 here; true sign is +1
    @test GO.rk_orient(m, (12.0,12.0), (24.0,24.0), (0.5, nextfloat(0.5)); exact = True()) > 0
end

@testset "rk_point_on_segment" begin
    @test GO.rk_point_on_segment(m, (0.5,0.5), (0.0,0.0), (1.0,1.0); exact = True()) == true
    @test GO.rk_point_on_segment(m, (2.0,2.0), (0.0,0.0), (1.0,1.0); exact = True()) == false   # collinear, outside
    @test GO.rk_point_on_segment(m, (0.5,0.6), (0.0,0.0), (1.0,1.0); exact = True()) == false
    @test GO.rk_point_on_segment(m, (0.0,0.0), (0.0,0.0), (1.0,1.0); exact = True()) == true    # endpoint inclusive
end

@testset "rk_point_in_ring" begin
    ring = GI.LinearRing([(0.0,0.0), (10.0,0.0), (10.0,10.0), (0.0,10.0), (0.0,0.0)])
    @test GO.rk_point_in_ring(m, (5.0,5.0), ring; exact = True()) == GO.LOC_INTERIOR
    @test GO.rk_point_in_ring(m, (5.0,0.0), ring; exact = True()) == GO.LOC_BOUNDARY
    @test GO.rk_point_in_ring(m, (15.0,5.0), ring; exact = True()) == GO.LOC_EXTERIOR
end

@testset "bounds" begin
    pa = GI.Polygon([[(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,0.0)]])
    ea = GO.rk_interaction_bounds(m, pa)
    @test !GO.rk_bounds_disjoint(ea, ea)
    @test GO.rk_bounds_covers(ea, ea)
end

@testset "rk_classify_intersection" begin
    cl(a0,a1,b0,b1) = GO.rk_classify_intersection(m, a0, a1, b0, b1; exact = True())
    # disjoint
    @test cl((0.,0.),(1.,0.),(0.,1.),(1.,1.)).kind == GO.SS_DISJOINT
    # proper crossing
    r = cl((0.,0.),(2.,2.),(0.,2.),(2.,0.))
    @test r.kind == GO.SS_PROPER
    @test !(r.a0_on_b || r.a1_on_b || r.b0_on_a || r.b1_on_a)
    # touch: b0 on interior of a
    r = cl((0.,0.),(2.,0.),(1.,0.),(1.,1.))
    @test r.kind == GO.SS_TOUCH && r.b0_on_a && !r.b1_on_a && !r.a0_on_b && !r.a1_on_b
    # touch: shared endpoint
    r = cl((0.,0.),(1.,0.),(1.,0.),(1.,1.))
    @test r.kind == GO.SS_TOUCH && r.a1_on_b && r.b0_on_a
    # collinear overlap
    r = cl((0.,0.),(2.,0.),(1.,0.),(3.,0.))
    @test r.kind == GO.SS_COLLINEAR && r.b0_on_a && r.a1_on_b
    # collinear disjoint
    @test cl((0.,0.),(1.,0.),(2.,0.),(3.,0.)).kind == GO.SS_DISJOINT
    # collinear, touching only at one shared endpoint -> SS_TOUCH
    r = cl((0.,0.),(1.,0.),(1.,0.),(2.,0.))
    @test r.kind == GO.SS_TOUCH && r.a1_on_b && r.b0_on_a
    # containment: b inside a (collinear)
    r = cl((0.,0.),(3.,0.),(1.,0.),(2.,0.))
    @test r.kind == GO.SS_COLLINEAR && r.b0_on_a && r.b1_on_a
    # degenerate: zero-length b on a
    r = cl((0.,0.),(2.,0.),(1.,0.),(1.,0.))
    @test r.kind == GO.SS_TOUCH && r.b0_on_a && r.b1_on_a

    # --- additional edge cases (beyond the plan) ---
    # T-touch: a0 on b's interior, non-collinear segments
    r = cl((1.,0.),(2.,1.),(0.,0.),(2.,0.))
    @test r.kind == GO.SS_TOUCH && r.a0_on_b && !r.a1_on_b && !r.b0_on_a && !r.b1_on_a
    # proper crossing with swapped sides: orientations (+,-)/(-,+) vs (-,+)/(+,-)
    r = cl((0.,2.),(2.,0.),(0.,0.),(2.,2.))
    @test r.kind == GO.SS_PROPER
    # almost-crossing: b straddles line(a) but a does not reach line(b)
    r = cl((0.,0.),(1.,1.),(3.,0.),(3.,4.))
    @test r.kind == GO.SS_DISJOINT
    @test !(r.a0_on_b || r.a1_on_b || r.b0_on_a || r.b1_on_a)
    # almost-crossing the other way: a straddles line(b) but b does not reach line(a)
    @test cl((3.,0.),(3.,4.),(0.,0.),(1.,1.)).kind == GO.SS_DISJOINT
    # both zero-length at the same point: degenerate touch, all flags true
    r = cl((1.,1.),(1.,1.),(1.,1.),(1.,1.))
    @test r.kind == GO.SS_TOUCH && r.a0_on_b && r.a1_on_b && r.b0_on_a && r.b1_on_a
    # both zero-length at different points (collinear trivially): disjoint
    r = cl((0.,0.),(0.,0.),(1.,0.),(1.,0.))
    @test r.kind == GO.SS_DISJOINT
    @test !(r.a0_on_b || r.a1_on_b || r.b0_on_a || r.b1_on_a)
    # collinear, shared endpoint, opposite directions (no interior overlap)
    r = cl((1.,0.),(0.,0.),(1.,0.),(2.,0.))
    @test r.kind == GO.SS_TOUCH && r.a0_on_b && r.b0_on_a && !r.a1_on_b && !r.b1_on_a
    # identical segments: collinear, all four flags
    r = cl((0.,0.),(1.,0.),(0.,0.),(1.,0.))
    @test r.kind == GO.SS_COLLINEAR && r.a0_on_b && r.a1_on_b && r.b0_on_a && r.b1_on_a

    # --- systematic symmetry checks over representative configurations ---
    seg_pairs = [
        ((0.,0.),(1.,0.),(0.,1.),(1.,1.)),   # disjoint
        ((0.,0.),(2.,2.),(0.,2.),(2.,0.)),   # proper crossing
        ((0.,0.),(2.,0.),(1.,0.),(1.,1.)),   # touch: b0 on interior of a
        ((0.,0.),(1.,0.),(1.,0.),(1.,1.)),   # touch: shared endpoint
        ((0.,0.),(2.,0.),(1.,0.),(3.,0.)),   # collinear overlap
        ((0.,0.),(1.,0.),(2.,0.),(3.,0.)),   # collinear disjoint
        ((0.,0.),(1.,0.),(1.,0.),(2.,0.)),   # collinear abutment
        ((0.,0.),(3.,0.),(1.,0.),(2.,0.)),   # collinear containment
        ((0.,0.),(2.,0.),(1.,0.),(1.,0.)),   # zero-length b on a
        ((1.,0.),(2.,1.),(0.,0.),(2.,0.)),   # T-touch: a0 on b's interior
        ((0.,0.),(1.,1.),(3.,0.),(3.,4.)),   # almost-crossing
        ((0.,0.),(1.,0.),(0.,0.),(1.,0.)),   # identical segments
    ]
    for (a0, a1, b0, b1) in seg_pairs
        r = cl(a0, a1, b0, b1)
        # swapping A and B: same kind, flags permuted
        s = cl(b0, b1, a0, a1)
        @test s.kind == r.kind
        @test (s.a0_on_b, s.a1_on_b, s.b0_on_a, s.b1_on_a) == (r.b0_on_a, r.b1_on_a, r.a0_on_b, r.a1_on_b)
        # reversing a's endpoints: kind invariant, a's flags swapped, b's unchanged
        v = cl(a1, a0, b0, b1)
        @test v.kind == r.kind
        @test (v.a0_on_b, v.a1_on_b, v.b0_on_a, v.b1_on_a) == (r.a1_on_b, r.a0_on_b, r.b0_on_a, r.b1_on_a)
    end
end

@testset "NodeKey" begin
    v = GO.vertex_node((1.0, 2.0))
    v2 = GO.vertex_node((1.0, 2.0))
    @test v == v2 && hash(v) == hash(v2)
    c1 = GO.crossing_node((0.,0.), (2.,2.), (0.,2.), (2.,0.))
    c2 = GO.crossing_node((0.,2.), (2.,0.), (2.,2.), (0.,0.))  # same pair, swapped & reversed
    @test c1 == c2 && hash(c1) == hash(c2)
    @test v != c1
    # signed zeros: -0.0 == 0.0 numerically but has different bits; constructors
    # must normalize so the default bit-pattern ==/hash agrees with coordinate equality
    z1 = GO.vertex_node((-0.0, 0.0))
    z2 = GO.vertex_node((0.0, -0.0))
    @test z1 == z2 && hash(z1) == hash(z2)
    cz1 = GO.crossing_node((-0.0,-2.), (0.,2.), (-2.,0.), (2.,-0.0))
    cz2 = GO.crossing_node((0.0,-2.), (-0.0,2.), (-2.,-0.0), (2.,0.))
    @test cz1 == cz2 && hash(cz1) == hash(cz2)
end

@testset "exact crossing coincidence (rational slow path)" begin
    # X crossing at exactly (1,1); a vertex node placed there must coincide
    c = GO.crossing_node((0.,0.), (2.,2.), (0.,2.), (2.,0.))
    @test GO.rk_nodes_coincide(m, c, GO.vertex_node((1.0, 1.0)); exact = True()) == true
    @test GO.rk_nodes_coincide(m, c, GO.vertex_node((1.0, 1.0 + eps(1.0))); exact = True()) == false
    # two crossings meeting at the same point
    c2 = GO.crossing_node((1.,0.), (1.,2.), (0.,1.), (2.,1.))
    @test GO.rk_nodes_coincide(m, c, c2; exact = True()) == true
    # crossing point with non-representable rational coordinates
    c3 = GO.crossing_node((0.,0.), (3.,1.), (0.,1.), (3.,0.))  # crosses at (1.5, 0.5)
    @test GO.rk_nodes_coincide(m, c3, GO.vertex_node((1.5, 0.5)); exact = True()) == true
end
