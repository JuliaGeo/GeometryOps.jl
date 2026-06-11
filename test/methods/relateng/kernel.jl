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
end
