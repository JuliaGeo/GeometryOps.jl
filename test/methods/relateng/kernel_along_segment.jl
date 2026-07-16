# Tests for the OverlayNG phase-1 kernel additions (design §2.5–§2.7):
# `rk_compare_along_segment` (planar + spherical, float filter vs exact
# fallback), the unified exact-crossing authority NodeKey call shape, and
# `_ring_material_interior_on_left`.

using Test
import GeometryOps as GO
import GeometryOps: Planar, Spherical, True, False
import GeometryOps.UnitSpherical: UnitSphericalPoint, UnitSphereFromGeographic
import GeoInterface as GI
using LinearAlgebra: cross, dot
using Random

# crossing node of horizontal segment (0,0)-(L,0) with the vertical through x.
_planar_xnode(L, x) = GO.crossing_node((0.0, 0.0), (L, 0.0), (x, -1.0), (x, 1.0))

@testset "rk_compare_along_segment (planar)" begin
    m = Planar()
    s0, s1 = (0.0, 0.0), (10.0, 0.0)
    k3 = _planar_xnode(10.0, 3.0)
    k7 = _planar_xnode(10.0, 7.0)
    v5 = GO.vertex_node((5.0, 0.0))

    @test GO.rk_compare_along_segment(m, s0, s1, k3, k7; exact = True()) == -1
    @test GO.rk_compare_along_segment(m, s0, s1, k7, k3; exact = True()) == 1
    @test GO.rk_compare_along_segment(m, s0, s1, k3, v5; exact = True()) == -1
    @test GO.rk_compare_along_segment(m, s0, s1, v5, k7; exact = True()) == -1
    #-- coincidence returns 0 (same node)
    @test GO.rk_compare_along_segment(m, s0, s1, k3, k3; exact = True()) == 0

    #-- reversed segment direction flips every order
    @test GO.rk_compare_along_segment(m, s1, s0, k3, k7; exact = True()) == 1
end

# exact-always along-parameter (rational) for a node key on segment (s0,s1)
function _exact_param(s0, s1, k)
    R = Rational{BigInt}
    dxr = R(s1[1]) - R(s0[1]); dyr = R(s1[2]) - R(s0[2])
    p = GO._exact_node_point(k)
    return (p[1] - R(s0[1])) * dxr + (p[2] - R(s0[2])) * dyr
end

@testset "planar order matches exact-always (dense)" begin
    m = Planar()
    rng = MersenneTwister(42)
    s0, s1 = (0.0, 0.0), (1.0, 0.0)
    for _ in 1:200
        xs = sort!(rand(rng, 6) .* 0.9 .+ 0.05)
        # crossing nodes at those xs, with random (well-conditioned) slopes
        ks = [GO.crossing_node(s0, s1, (x, -rand(rng) - 0.2), (x, rand(rng) + 0.2)) for x in xs]
        # comparator sort
        order_filter = sortperm(collect(1:length(ks));
            lt = (i, j) -> GO.rk_compare_along_segment(m, s0, s1, ks[i], ks[j]; exact = True()) < 0)
        order_exact = sortperm([_exact_param(s0, s1, k) for k in ks])
        @test order_filter == order_exact
        # since xs are strictly increasing, both must equal 1:length
        @test order_filter == collect(1:length(ks))
    end
end

@testset "planar near-parallel escalates soundly" begin
    m = Planar()
    s0, s1 = (0.0, 0.0), (1.0, 0.0)
    # two crossings extremely close together, from near-parallel secondaries:
    # the filter must return the same sign as the exact parameters.
    x1 = 0.5
    x2 = nextfloat(nextfloat(0.5))
    ka = GO.crossing_node(s0, s1, (x1, -1e-9), (x1, 1e-9))
    kb = GO.crossing_node(s0, s1, (x2, -1e-9), (x2, 1e-9))
    got = GO.rk_compare_along_segment(m, s0, s1, ka, kb; exact = True())
    want = _exact_param(s0, s1, ka) < _exact_param(s0, s1, kb) ? -1 : 1
    @test got == want
end

@testset "rk_compare_along_segment (spherical)" begin
    ms = Spherical()
    usp(lon, lat) = UnitSphereFromGeographic()((lon, lat))
    s0, s1 = usp(0.0, 0.0), usp(10.0, 0.0)
    k3 = GO.crossing_node(s0, s1, usp(3.0, -1.0), usp(3.0, 1.0))
    k7 = GO.crossing_node(s0, s1, usp(7.0, -1.0), usp(7.0, 1.0))
    v5 = GO.vertex_node(usp(5.0, 0.0))

    @test GO.rk_compare_along_segment(ms, s0, s1, k3, k7; exact = True()) == -1
    @test GO.rk_compare_along_segment(ms, s0, s1, k7, k3; exact = True()) == 1
    @test GO.rk_compare_along_segment(ms, s0, s1, k3, v5; exact = True()) == -1
    @test GO.rk_compare_along_segment(ms, s0, s1, v5, k7; exact = True()) == -1
    @test GO.rk_compare_along_segment(ms, s0, s1, k3, k3; exact = True()) == 0
    @test GO.rk_compare_along_segment(ms, s1, s0, k3, k7; exact = True()) == 1

    #-- dense monotone check: crossings at increasing longitudes order correctly
    lons = collect(1.0:0.5:9.0)
    ks = [GO.crossing_node(s0, s1, usp(x, -1.0), usp(x, 1.0)) for x in lons]
    order = sortperm(collect(1:length(ks));
        lt = (i, j) -> GO.rk_compare_along_segment(ms, s0, s1, ks[i], ks[j]; exact = True()) < 0)
    @test order == collect(1:length(ks))
end

@testset "unified exact-crossing authority" begin
    k = _planar_xnode(10.0, 4.0)
    @test GO._exact_crossing_point(k) == GO._exact_crossing_point((0.0, 0.0), (10.0, 0.0), (4.0, -1.0), (4.0, 1.0))
    @test GO._exact_node_point(k) == GO._exact_crossing_point(k)
    v = GO.vertex_node((2.0, 3.0))
    @test GO._exact_node_point(v) == (Rational{BigInt}(2), Rational{BigInt}(3))
    @test Float64.(GO._exact_crossing_point(k)) == (4.0, 0.0)
end

@testset "_ring_material_interior_on_left flip" begin
    for m in (Planar(), Spherical())
        # a CW square (shell) and its use as a hole
        ptsll = [(0.0, 0.0), (0.0, 4.0), (4.0, 4.0), (4.0, 0.0), (0.0, 0.0)]
        pts = m isa Planar ? ptsll : [UnitSphereFromGeographic()(p) for p in ptsll]
        il = GO._ring_interior_on_left(m, pts, false; exact = True())
        @test GO._ring_material_interior_on_left(m, pts, false; exact = True()) == il
        @test GO._ring_material_interior_on_left(m, pts, true; exact = True()) == !il
    end
end
