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
