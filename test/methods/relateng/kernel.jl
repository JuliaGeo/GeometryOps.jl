# Tests for the planar RelateKernel: orientation, point-on-segment,
# point-in-ring and interaction bounds (design doc D1/D2).

using Test
import GeometryOps as GO
import GeometryOps: Planar, True, False
import GeometryOps.UnitSpherical: UnitSphericalPoint
import GeoInterface as GI
import Extents

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
    @test Extents.intersects(ea, ea)
    @test Extents.covers(ea, ea)
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

@testset "rk_quadrant" begin
    # JTS Quadrant convention: NE=0, NW=1, SW=2, SE=3, numbered CCW from the
    # positive X-axis; axis directions belong to the `>= 0` side (so +x and
    # +y are NE, -x is NW, -y is SE).
    o = (0.0, 0.0)
    @test GO.rk_quadrant(m, o, (1.,0.)) == 0    # +x axis is NE
    @test GO.rk_quadrant(m, o, (1.,1.)) == 0
    @test GO.rk_quadrant(m, o, (0.,1.)) == 0    # +y axis is NE
    @test GO.rk_quadrant(m, o, (-1.,1.)) == 1
    @test GO.rk_quadrant(m, o, (-1.,0.)) == 1   # -x axis is NW
    @test GO.rk_quadrant(m, o, (-1.,-1.)) == 2
    @test GO.rk_quadrant(m, o, (0.,-1.)) == 3   # -y axis is SE
    @test GO.rk_quadrant(m, o, (1.,-1.)) == 3
    @test GO.rk_quadrant(m, (2.0, 3.0), (3.0, 3.0)) == 0   # non-origin apex
    @test_throws ArgumentError GO.rk_quadrant(m, o, (0.0, 0.0))   # zero-length direction
end

@testset "edge ordering around a vertex node" begin
    # Sign expectations verified against the Java contract
    # (PolygonNodeTopology.compareAngle): returns negative / zero / positive
    # as angle(P) is less than / equal to / greater than angle(Q), with
    # angles increasing CCW from the positive X-axis. Different quadrants
    # decide the comparison; same-quadrant ties are resolved by
    # Orientation.index(origin, q, p) (CCW -> P greater).
    origin = GO.vertex_node((0.0, 0.0))
    east, north, west, south = (1.,0.), (0.,1.), (-1.,0.), (0.,-1.)
    vcmp(p, q) = GO.rk_compare_edge_dir(m, origin, p, q; exact = True())
    @test vcmp(east, north) < 0   # same quadrant (NE holds both + axes), orientation resolves
    @test vcmp(north, west) < 0   # NE=0 < NW=1
    @test vcmp(west, south) < 0   # NW=1 < SE=3 (south lies on the dx>=0 side)
    @test vcmp(east, east) == 0
    @test vcmp(north, east) > 0
    # same quadrant resolved by orientation
    @test vcmp((2.0, 1.0), (1.0, 2.0)) < 0
    # full fan in strictly increasing angular order: comparator must be a
    # total order matching position in the fan (antisymmetry included)
    fan = [(1.,0.), (2.,1.), (1.,1.), (1.,2.), (0.,1.),            # NE: 0..90
           (-1.,2.), (-1.,1.), (-2.,1.), (-1.,0.),                  # NW: (90)..180
           (-2.,-1.), (-1.,-1.), (-1.,-2.),                         # SW: (180)..(270)
           (0.,-1.), (1.,-2.), (1.,-1.), (2.,-1.)]                  # SE: 270..(360)
    for i in eachindex(fan), j in eachindex(fan)
        @test vcmp(fan[i], fan[j]) == (i == j ? 0 : (i < j ? -1 : 1))
    end
end

@testset "crossing-node incident edge order" begin
    # a: (0,0)->(2,2), b: (0,2)->(2,0); crossing at symbolic (1,1)
    dirs = GO.rk_crossing_dirs_ccw(m, (0.,0.), (2.,2.), (0.,2.), (2.,0.); exact = True())
    # CCW order starting from direction toward a1=(2,2):
    @test dirs == ((2.,2.), (0.,2.), (0.,0.), (2.,0.))
    # reversing b gives the same cyclic order (other orientation branch)
    dirs2 = GO.rk_crossing_dirs_ccw(m, (0.,0.), (2.,2.), (2.,0.), (0.,2.); exact = True())
    @test dirs2 == ((2.,2.), (0.,2.), (0.,0.), (2.,0.))

    # rk_compare_edge_dir with a crossing apex reproduces compareAngle
    # anchored at the positive X-axis of the (symbolic) crossing point (1,1):
    # angles from the apex: (2,2)=45deg, (0,2)=135deg, (0,0)=225deg, (2,0)=315deg
    cn = GO.crossing_node((0.,0.), (2.,2.), (0.,2.), (2.,0.))
    ccmp(p, q) = GO.rk_compare_edge_dir(m, cn, p, q; exact = True())
    order = [(2.,2.), (0.,2.), (0.,0.), (2.,0.)]
    for i in eachindex(order), j in eachindex(order)
        @test ccmp(order[i], order[j]) == (i == j ? 0 : (i < j ? -1 : 1))
    end
    # directions which are not incident endpoints (D3 coincidence-merged
    # nodes) fall back to the exact rational apex (1,1): (5,5) lies on the
    # same ray from the apex as (2,2)
    @test ccmp((5.,5.), (2.,2.)) == 0
    @test ccmp((5.,5.), (0.,2.)) < 0
    @test ccmp((0.,2.), (5.,5.)) > 0
    # same-quadrant pair around a crossing apex, resolved by orientation:
    # a: (0,0)->(4,1) x b: (1,-1)->(2,4) cross properly at (24/19, 6/19);
    # both a1=(4,1) and b1=(2,4) are NE of the apex, with angle(4,1) smaller
    cn2 = GO.crossing_node((0.,0.), (4.,1.), (1.,-1.), (2.,4.))
    @test GO.rk_compare_edge_dir(m, cn2, (4.,1.), (2.,4.); exact = True()) < 0
    @test GO.rk_compare_edge_dir(m, cn2, (2.,4.), (4.,1.); exact = True()) > 0
    @test GO.rk_compare_edge_dir(m, cn2, (2.,4.), (2.,4.); exact = True()) == 0
end

@testset "isCrossing / isInteriorSegment" begin
    n = (1.0, 1.0)
    @test GO.rk_is_crossing(m, GO.vertex_node(n), (0.,0.), (2.,2.), (0.,2.), (2.,0.); exact = True())
    @test !GO.rk_is_crossing(m, GO.vertex_node(n), (0.,0.), (2.,2.), (2.,0.), (2.,2.); exact = True())
    # both b-arms on the same side of the a-corner: a touch, not a crossing
    @test !GO.rk_is_crossing(m, GO.vertex_node(n), (0.,0.), (2.,2.), (0.,2.), (1.,2.); exact = True())
    # collinear arm -> reported as not crossing (Java contract)
    @test !GO.rk_is_crossing(m, GO.vertex_node(n), (0.,0.), (2.,2.), (0.,0.), (2.,0.); exact = True())
    # crossing-node apexes are rejected: proper crossings are crossings by
    # construction (JTS TopologyComputer.updateAreaAreaCross short-circuits
    # them with `isProper ||` before ever calling isCrossing)
    cn = GO.crossing_node((0.,0.), (2.,2.), (0.,2.), (2.,0.))
    @test_throws ArgumentError GO.rk_is_crossing(m, cn, (0.,0.), (2.,2.), (0.,2.), (2.,0.); exact = True())

    # isInteriorSegment: corner a0 -> node -> a1, ring interior on the right
    nd = GO.vertex_node((0.0, 0.0))
    @test !GO.rk_is_interior_segment(m, nd, (0.,1.), (1.,0.), (1.,1.); exact = True())
    @test GO.rk_is_interior_segment(m, nd, (0.,1.), (1.,0.), (-1.,-1.); exact = True())
    # reversed corner flips the interior side
    @test GO.rk_is_interior_segment(m, nd, (1.,0.), (0.,1.), (1.,1.); exact = True())
    @test !GO.rk_is_interior_segment(m, nd, (1.,0.), (0.,1.), (-1.,-1.); exact = True())
    @test_throws ArgumentError GO.rk_is_interior_segment(m, cn, (0.,0.), (2.,2.), (1.,1.); exact = True())
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

@testset "node helpers are N-dimensional (3D round-trips)" begin
    u = UnitSphericalPoint(0.0, -0.0, 1.0)
    k = GO.vertex_node(u)
    # signed zero normalized away, all three components preserved
    @test GI.x(k.pt) == 0.0 && GI.y(k.pt) == 0.0 && GI.z(k.pt) == 1.0
    @test !signbit(GI.y(k.pt))               # -0.0 -> +0.0
    # planar 2-tuple path unchanged
    k2 = GO.vertex_node((3.0, 4.0))
    @test k2.pt == (3.0, 4.0)
end

@testset "crossing_node canonicalizes in 3D and is order-invariant" begin
    a0 = UnitSphericalPoint(1.0, 0.0, 0.0); a1 = UnitSphericalPoint(0.0, 1.0, 0.0)
    b0 = UnitSphericalPoint(0.0, 0.0, 1.0); b1 = UnitSphericalPoint(1.0, 1.0, 1.0)
    k1 = GO.crossing_node(a0, a1, b0, b1)
    k2 = GO.crossing_node(b1, b0, a1, a0)   # same pair, every order/orientation flipped
    @test k1 == k2
    @test k1.pt isa UnitSphericalPoint       # USP preserved through canonicalization
end
