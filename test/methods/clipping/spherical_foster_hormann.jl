#=
# Spherical Foster-Hormann clipping

`FosterHormannClipping(Spherical())` reads polygon edges as great-circle arcs. These tests
cover what that changes relative to the planar algorithm, and pin the two properties the
spherical path is easiest to silently lose: that the planar path is untouched, and that the
answer is still right at the scale of a discrete-global-grid cell rather than only at the
scale of a country.

The oracle is s2geography (Google S2), through the same FFI the OverlayNG differential suite
uses. RelateNG is deliberately *not* used as an oracle here: its spherical predicates are
wrong on exact shared-edge contact, which is the configuration a DGG tiling produces
constantly.
=#

using Test
import GeometryOps as GO
import GeoInterface as GI
import GeometryOpsCore
using GeometryOpsCore: Planar, Spherical, Geodesic
using Random
using LinearAlgebra: norm

include(joinpath(@__DIR__, "..", "..", "external", "s2geography", "s2geography.jl"))
using .S2Geog
using GeometryOpsTestHelpers: write_wkb

const S2_OK = s2_available()
const RADIUS_RATIO = (S2Geog.S2_RADIUS / GeometryOpsCore.WGS84_EARTH_MEAN_RADIUS)^2

const ALG_S = GO.FosterHormannClipping(Spherical())
const ALG_P = GO.FosterHormannClipping(Planar())

ring(pts) = GI.Polygon([GI.LinearRing(vcat(pts, [pts[1]]))])

#-- FH returns a vector of polygons; s2 wants one geometry
function as_mpoly(ps)
    isempty(ps) && return nothing
    GI.MultiPolygon([GI.Polygon([GI.LinearRing(collect(GI.getpoint(r))) for r in GI.getring(p)])
                     for p in ps])
end

fh(op, alg, a, b) = (op === :intersection ? GO.intersection :
                     op === :union ? GO.union : GO.difference)(alg, a, b; target = GI.PolygonTrait())

#= Relative area disagreement with s2 for one op, normalized by the size of the inputs so
that a near-empty result is not judged against its own vanishing area. =#
function s2_area_disagreement(op, a, b)
    wa, wb = write_wkb(a), write_wkb(b)
    scale = s2_area(wa) + s2_area(wb)
    scale > 0 || return 0.0
    ours = as_mpoly(fh(op, ALG_S, a, b))
    a_ours = ours === nothing ? 0.0 : s2_area(write_wkb(ours))
    return abs(a_ours - s2_area(s2_overlay(op, wa, wb))) / scale
end

#= An irregular, non-convex cell. Irregular on purpose: a regular n-gon is symmetric enough
to hide an ordering or ulp bug at a shared vertex. `notch` makes one vertex turn inward, so
the ring is genuinely non-convex — the case a HEALPix or ISEA4R tiling forces, and the case
`ConvexConvexSutherlandHodgman` is documented to give undefined results for. =#
function cell(rng, lon0, lat0, s; dup = false)
    base = [(0.0, 0.0), (1.0, 0.13), (2.0, -0.07), (2.1, 1.0), (1.05, 0.55), (0.0, 1.0)]
    pts = [(lon0 + s * (u + 0.07 * (rand(rng) - 0.5)), lat0 + s * (v + 0.07 * (rand(rng) - 0.5)))
           for (u, v) in base]
    #-- HEALPix rings carry coincident vertices near the polar corners
    dup && insert!(pts, 4, pts[4])
    push!(pts, pts[1])
    GI.Polygon([GI.LinearRing(pts)])
end

@testset "constructing a spherical algorithm" begin
    #= Every one of these threw `MethodError: ... is ambiguous` before the constructor was
    narrowed, which made the whole spherical path unreachable: the `Union{Spherical,
    Geodesic}` method was narrower in the manifold but wider in the accelerator than the
    struct's own outer constructor, so neither won. =#
    for alg in (GO.FosterHormannClipping(Spherical()),
                GO.FosterHormannClipping(; manifold = Spherical()),
                GO.FosterHormannClipping(Spherical(), GO.NestedLoop()),
                GO.FosterHormannClipping(Spherical(), GO.AutoAccelerator()))
        @test alg isa GO.FosterHormannClipping
        #-- STRtrees index planar rectangles; the sphere always falls back to a nested loop
        @test alg.accelerator isa GO.NestedLoop
    end
    #-- the planar side still honours an explicit accelerator choice
    @test GO.FosterHormannClipping(Planar(), GO.AutoAccelerator()).accelerator isa GO.AutoAccelerator
end

@testset "Geodesic manifold is rejected at construction" begin
    #= Foster-Hormann has no geodesic implementation of `_get_side` and the other clipping
    primitives. `FosterHormannClipping(Geodesic())` used to construct fine and only fail
    with a bare `MethodError` deep inside the first clip that reached the missing method —
    far from the cause, and only on inputs that exercised that code path. Every spelling
    that would build a Geodesic algorithm must fail immediately and clearly instead. =#
    @test_throws ArgumentError GO.FosterHormannClipping(Geodesic())
    @test_throws ArgumentError GO.FosterHormannClipping(; manifold = Geodesic())
    @test_throws ArgumentError GO.FosterHormannClipping(Geodesic(), GO.NestedLoop())
    @test_throws ArgumentError GO.FosterHormannClipping(Geodesic(), GO.AutoAccelerator())
    #-- direct parametric construction bypasses none of the outer constructors above
    @test_throws ArgumentError GO.FosterHormannClipping{Geodesic{Float64}, GO.NestedLoop}(Geodesic(), GO.NestedLoop())
    #-- Planar and Spherical are unaffected
    @test GO.FosterHormannClipping(Planar()) isa GO.FosterHormannClipping
    @test GO.FosterHormannClipping(Spherical()) isa GO.FosterHormannClipping
end

@testset "crossings land on the great circle, not the chart line" begin
    #= Two lon/lat squares. Their boundaries are not great circles, so the spherical
    crossings are *not* the planar ones: the top edge of A from (0,10) to (10,10) bulges
    north of the parallel, and meets the meridian lon=5 above latitude 10.

    This is the test that fails if `_intersection_point` silently falls back to the planar
    method — which it did, despite taking a manifold argument. =#
    A = ring([(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
    B = ring([(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 15.0)])

    pl = collect(GI.getpoint(GI.getring(fh(:intersection, ALG_P, A, B)[1], 1)))
    sp = collect(GI.getpoint(GI.getring(fh(:intersection, ALG_S, A, B)[1], 1)))

    @test (5.0, 10.0) in pl                     # planar: the chart crossing
    @test !((5.0, 10.0) in sp)                  # spherical: strictly north of it
    #-- the northern crossing on the lon=5 meridian; (5,5) is the other vertex at that lon
    lat_at_5 = only(unique(p[2] for p in sp if isapprox(p[1], 5.0; atol = 1e-9) && p[2] > 9.0))
    @test lat_at_5 > 10.0
    @test isapprox(lat_at_5, 10.0374230459; atol = 1e-8)

    if S2_OK
        @test s2_area_disagreement(:intersection, A, B) < 1e-12
    end
end

@testset "degenerate rings do not divide by zero" begin
    rng = MersenneTwister(20260823)
    for lat0 in (0.0, 84.9, -84.9), s in (2.3e-4, 1.0)
        P = cell(rng, 30.0, lat0, s; dup = true)
        Q = cell(rng, 30.0 + s * 0.9, lat0 + s * 0.35, s; dup = true)
        for op in (:intersection, :union, :difference)
            res = fh(op, ALG_S, P, Q)
            for p in res, r in GI.getring(p), pt in GI.getpoint(r)
                @test isfinite(GI.x(pt)) && isfinite(GI.y(pt))
            end
        end
    end
end

@testset "adjacent cells sharing an edge" begin
    #= Two cells that share a full edge: the common case in a tiling, and the one where a
    duplicated intersection point used to derail the traversal into a `TracingError`. =#
    shared = [(20.0, 10.0), (20.4, 10.9)]
    P = ring([(19.2, 10.2), shared[1], shared[2], (19.3, 11.1)])
    Q = ring([shared[1], (21.1, 9.9), (21.2, 10.8), shared[2]])
    for alg in (ALG_P, ALG_S)
        @test GO.intersection_area(alg, P, Q) ≈ 0 atol = 1e-6 * GO.area(alg.manifold, P)
    end
    if S2_OK
        #-- the union of two edge-sharing cells is their sum: no sliver, no double count
        wa, wb = write_wkb(P), write_wkb(Q)
        u = as_mpoly(fh(:union, ALG_S, P, Q))
        @test u !== nothing
        @test isapprox(s2_area(write_wkb(u)), s2_area(wa) + s2_area(wb); rtol = 1e-12)
    end
end

@testset "near-antipodal cells do not cross" begin
    #= Two great circles always meet, at an antipodal pair. Each circle separating the
    other's endpoints is therefore not enough to say the two *arcs* meet: each arc can
    reach a different one of the two meeting points.

    A cell near lon 0 and a cell near lon 180 are near-antipodal images of one another, so
    each lies almost on the other's great circle and straddles it. Under the two-way
    straddle test they crossed, and two cells on opposite sides of the globe reported an
    intersection the area of a whole cell. =#
    sq(lon, lat, d) = ring([(lon-d, lat-d), (lon+d, lat-d), (lon+d, lat+d), (lon-d, lat+d)])
    for d in (0.25, 0.5, 1.0, 2.0), off in (178.0, 179.0, 179.9, 179.99, 179.999, 180.0, 180.001, 180.5, 181.0)
        @test GO.intersection_area(ALG_S, sq(0.0, 0.0, d), sq(off, 0.0, d)) == 0.0
    end
    #-- and a genuine crossing is still found
    @test GO.intersection_area(ALG_S, sq(0.0, 0.0, 1.0), sq(0.5, 0.0, 1.0)) > 0.0
end

@testset "collinear-point removal is manifold aware" begin
    #= A run of vertices along a parallel is exactly collinear in the chart but not on a
    great circle. Dropping them as redundant — which the planar rule does — moves the
    boundary poleward by the arc's sagitta and loses the sliver between. Canada and the
    United States meet along the 49th parallel, so this was worth ~0.1% of Canada. =#
    north = ring([(-123.0, 49.0), (-110.0, 49.0), (-100.0, 49.0), (-95.0, 49.0), (-95.0, 60.0), (-123.0, 60.0)])
    south = ring([(-123.0, 40.0), (-95.0, 40.0), (-95.0, 49.0), (-100.0, 49.0), (-110.0, 49.0), (-123.0, 49.0)])

    res = fh(:difference, ALG_S, north, south)
    @test length(res) == 1
    lats = [GI.y(pt) for r in GI.getring(res[1]) for pt in GI.getpoint(r)]
    #-- the intermediate vertices of the shared run must survive on the sphere
    @test count(≈(49.0), lats) >= 4

    if S2_OK
        @test s2_area_disagreement(:difference, north, south) < 1e-12
    end
end

#= Two cells that share a chart edge and overlap only slightly.

A parallel is not a great circle, so the shared horizontal edges of two such cells are arcs
of *different* great circles and genuinely cross. The banded `spherical_orient` reports 0
inside an eps*16 window, which at cell scale is wider than the determinant it is judging, so
the crossing was classified as a hinge, the entry/exit alternation collapsed, and the tracer
returned the whole subject ring -- up to 100,000x the true overlap area.

This is the configuration a DGG tiling produces constantly, and the random-offset cells in
the differential corpus never generate it: they share no edge exactly. =#
@testset "shared edge with a thin overlap" begin
    sq(lon, lat, w, h) = GI.Polygon([GI.LinearRing(
        [(lon, lat), (lon + w, lat), (lon + w, lat + h), (lon, lat + h), (lon, lat)])])

    for s in (2.3e-4, 2.3e-3), frac in (0.05, 0.01, 0.002, 1e-3, 1e-4, 1e-5)
        P = sq(0.0, 12.0, s, s)
        Q = sq(s * (1 - frac), 12.0, s, s)   # shares both horizontal edges exactly
        expected = GO.area(Spherical(), P) * frac
        got = GO.intersection_area(ALG_S, P, Q)
        @test got ≈ expected rtol = 1e-3
    end

    #- and the same through the geometry-returning path
    s, frac = 2.3e-4, 1e-4
    P = sq(0.0, 12.0, s, s); Q = sq(s * (1 - frac), 12.0, s, s)
    polys = GO.intersection(ALG_S, P, Q; target = GI.PolygonTrait())
    @test sum(p -> GO.area(Spherical(), p), polys) ≈ GO.area(Spherical(), P) * frac rtol = 1e-3
end

#= A vertex that passes through the clip untouched must come back exactly as it arrived.

`Spherical` computes in `UnitSphericalPoint`. 3D input therefore never converts at all. 2D
lon/lat input converts once at ingress, and would come back as
`GeographicFromUnitSphere(UnitSphereFromGeographic(p))` -- displaced by a few ulps, ~6e-9 m
-- if egress round-tripped it. It does not: a passthrough vertex is fetched from the input
ring by its source index, so it is returned as the very value that went in.

The test is not "is it close": it asserts that any output vertex landing near an input
vertex is *exactly* that vertex, which is what a round trip would break. =#
@testset "passthrough vertices are bit-exact" begin
    usp(lon, lat) = GO.UnitSpherical.UnitSphereFromGeographic()((lon, lat))

    #-- nearest input vertex, and whether the match is exact
    function check(outpts, inpts, dist)
        near = 0
        for o in outpts
            d, best = Inf, nothing
            for i in inpts
                dd = dist(o, i)
                dd < d && (d = dd; best = i)
            end
            #-- "near" means far closer than any real crossing could be to a vertex
            if d < 1e-9
                near += 1
                @test o === best || o == best
            end
        end
        return near
    end

    for lon0 in (0.0, 120.0)
        pa = [(lon0, 0.0), (lon0 + 10.0, 0.0), (lon0 + 10.0, 10.0), (lon0, 10.0)]
        pb = [(lon0 + 5.0, 5.0), (lon0 + 15.0, 5.0), (lon0 + 15.0, 15.0), (lon0 + 5.0, 15.0)]

        #-- 2D lon/lat in, lon/lat out
        A2, B2 = ring(pa), ring(pb)
        out2 = GO.intersection(ALG_S, A2, B2; target = GI.PolygonTrait())
        pts2 = collect(GI.getpoint(GI.getexterior(out2[1])))
        @test eltype(pts2) <: Tuple
        n2 = check(pts2, vcat(pa, pb), (o, i) -> hypot(o[1] - i[1], o[2] - i[2]))
        @test n2 ≥ 2   # both shared corners survive the clip

        #-- 3D xyz in, xyz out
        A3 = GI.Polygon([GI.LinearRing([usp(p...) for p in vcat(pa, [pa[1]])])])
        B3 = GI.Polygon([GI.LinearRing([usp(p...) for p in vcat(pb, [pb[1]])])])
        out3 = GO.intersection(ALG_S, A3, B3; target = GI.PolygonTrait())
        pts3 = collect(GI.getpoint(GI.getexterior(out3[1])))
        @test eltype(pts3) <: GO.UnitSpherical.UnitSphericalPoint
        ins3 = [usp(p...) for p in vcat(pa, pb)]
        n3 = check(pts3, ins3, (o, i) -> norm(o .- i))
        @test n3 ≥ 2

        #-- and the two representations agree on the answer
        @test GO.intersection_area(ALG_S, A2, B2) ≈ GO.intersection_area(ALG_S, A3, B3) rtol = 1e-12
    end
end

#= No-crossing containment is the *normal* case in conservative regridding -- a fine cell
inside a coarse one has no crossings at all -- yet it is the one answer the tracer does not
build. With no crossings there are no nodes to trace, so the ops rebuild the answer from the
input geometry instead, and that rebuild has to land in the representation the tracer would
have emitted or it does not fit the vector it is pushed into.

A corpus of randomly offset overlapping cells always crosses, so it can never reach this
branch. These cases do, on both input representations. =#
@testset "containment and disjoint, both input representations" begin
    usp(lon, lat) = GO.UnitSpherical.UnitSphereFromGeographic()((lon, lat))
    box(c, h) = [(c[1] - h, c[2] - h), (c[1] + h, c[2] - h), (c[1] + h, c[2] + h), (c[1] - h, c[2] + h)]
    as_usp(pts) = GI.Polygon([GI.LinearRing([usp(p...) for p in vcat(pts, [pts[1]])])])

    small_pts, big_pts, away_pts = box((10.0, 45.0), 0.25), box((10.0, 45.0), 1.0), box((40.0, 45.0), 0.25)
    areas = Float64[]

    for (build, P) in ((ring, Tuple), (as_usp, GO.UnitSpherical.UnitSphericalPoint))
        small, big, away = build(small_pts), build(big_pts), build(away_pts)
        a_small, a_big = GO.area(Spherical(), small), GO.area(Spherical(), big)

        #-- the intersection of a contained pair is exactly the contained polygon, either order
        @test GO.intersection_area(ALG_S, small, big) == a_small
        @test GO.intersection_area(ALG_S, big, small) == a_small
        @test GO.intersection_area(ALG_S, small, away) == 0
        push!(areas, a_small, GO.intersection_area(ALG_S, small, big))

        #-- emitted in the representation the caller supplied, not silently converted
        got = GO.intersection(ALG_S, small, big; target = GI.PolygonTrait())
        @test length(got) == 1
        pts = collect(GI.getpoint(GI.getexterior(got[1])))
        @test eltype(pts) <: P
        #-- rebuilt from the input ring, so bit-exact rather than round-tripped
        @test Set(pts) == Set(GI.getpoint(GI.getexterior(small)))

        #-- union of a contained pair is the container; of a disjoint pair, two pieces
        u = GO.union(ALG_S, small, big; target = GI.PolygonTrait())
        @test length(u) == 1
        @test GO.area(Spherical(), u[1]) == a_big
        @test length(GO.union(ALG_S, small, away; target = GI.PolygonTrait())) == 2

        #-- container minus contained leaves a hole; disjoint leaves the subject whole
        d = GO.difference(ALG_S, big, small; target = GI.PolygonTrait())
        @test length(d) == 1
        @test GI.nhole(d[1]) == 1
        dd = GO.difference(ALG_S, small, away; target = GI.PolygonTrait())
        @test length(dd) == 1
        @test GO.area(Spherical(), dd[1]) == a_small
    end

    #-- and the two representations agree exactly, not just to a tolerance
    @test areas[1] == areas[3]
    @test areas[2] == areas[4]
end

@testset "planar path is untouched" begin
    #= Pinned planar results. The spherical work is purely additive; if any of these move,
    a manifold-aware branch has leaked into the planar path. =#
    A = ring([(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
    B = ring([(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 15.0)])
    pts = collect(GI.getpoint(GI.getring(fh(:intersection, ALG_P, A, B)[1], 1)))
    @test pts == [(10.0, 5.0), (10.0, 10.0), (5.0, 10.0), (5.0, 5.0), (10.0, 5.0)]
    @test GO.intersection_area(ALG_P, A, B) == 25.0

    #-- planar collinear removal still collapses a straight run
    line = ring([(0.0, 0.0), (5.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
    box  = ring([(2.0, -1.0), (8.0, -1.0), (8.0, 11.0), (2.0, 11.0)])
    got = collect(GI.getpoint(GI.getring(fh(:intersection, ALG_P, line, box)[1], 1)))
    @test !any(p -> p == (5.0, 0.0), got)
end

@testset "FosterHormannCache" begin
    rng = MersenneTwister(7)
    prs = [(cell(rng, 120.0 + i * 1e-3, 12.0, 2.3e-4),
            cell(rng, 120.0 + i * 1e-3 + 2.07e-4, 12.0 + 8e-5, 2.3e-4)) for i in 1:40]
    append!(prs, [(cell(rng, 10.0 + i * 0.01, 20.0, 1.0),
                   cell(rng, 10.0 + i * 0.01 + 0.9, 20.35, 1.0)) for i in 1:40])

    for alg in (ALG_P, ALG_S)
        cache = GO.FosterHormannCache(alg)
        #-- a cache must not change the answer, to the last bit
        for (a, b) in prs
            @test GO.intersection_area(alg, a, b; cache) === GO.intersection_area(alg, a, b)
        end
        #-- and it must actually remove the per-call allocation of the working set
        a, b = prs[1]
        GO.intersection_area(alg, a, b; cache)
        uncached = @allocated GO.intersection_area(alg, a, b)
        cached = @allocated GO.intersection_area(alg, a, b; cache)
        #= The cache must remove the per-call working set. Julia < 1.12 does not
        elide those allocations, so only the weaker ordering holds there. =#
        @test cached < uncached
        @test cached < uncached ÷ 4 skip = VERSION < v"1.12"
    end

    #-- the float type has to match, and says so
    @test_throws ArgumentError GO.intersection_area(ALG_S, prs[1]..., Float64;
        cache = GO.FosterHormannCache(Float32))
end

if S2_OK
    @testset "spherical FH vs s2geography — cell scale" begin
        #= Cell scale is the point. A HEALPix level-18 edge is ~4e-6 rad; at that size the
        conditioning of the crossing construction dominates, and a result measured only on
        degree-scale polygons says nothing about it. =#
        rng = MersenneTwister(4242)
        worst = Dict(2.3e-4 => 0.0, 2.3e-3 => 0.0, 1.0 => 0.0)
        for s in (2.3e-4, 2.3e-3, 1.0), lat0 in (0.0, 45.0, 84.9, -84.9), dup in (false, true), _ in 1:3
            lon0 = rand(rng) * 360 - 180
            P = cell(rng, lon0, lat0, s; dup)
            Q = cell(rng, lon0 + s * 0.9, lat0 + s * 0.35, s; dup)
            for op in (:intersection, :union, :difference)
                worst[s] = max(worst[s], s2_area_disagreement(op, P, Q))
            end
        end
        #= Degree scale runs at machine precision. Cell scale is limited by `PolyNode`
        storing lon/lat in Float64: a crossing is quantized to ~1e-14 of a degree, which is
        ~1e-10 of a level-18 cell, and that floor is what these bounds leave room for. =#
        @test worst[1.0] < 1e-11
        @test worst[2.3e-3] < 1e-7
        @test worst[2.3e-4] < 1e-7
    end

    @testset "spherical FH agrees with ConvexConvexSutherlandHodgman on convex input" begin
        #= Where both algorithms are valid they must agree; divergence means one is broken.
        The convex clipper takes unit-sphere xyz, the Foster-Hormann one lon/lat, so the
        comparison goes through the chart conversion. =#
        to_usp(g) = GO.apply(GI.PointTrait(), g) do p
            GO.UnitSpherical.UnitSphereFromGeographic()((GI.x(p), GI.y(p)))
        end
        to_ll(g) = GO.apply(GI.PointTrait(), g) do p
            GO._usp_to_lonlat(p)
        end
        sh = GO.ConvexConvexSutherlandHodgman(Spherical())
        rng = MersenneTwister(20260823)
        mk(n, lon0, lat0, r) = ring([(lon0 + r * cospi(2t), lat0 + r * sinpi(2t)) for t in sort(rand(rng, n))])
        worst, n = 0.0, 0
        for _ in 1:60
            lon0, lat0 = rand(rng) * 300 - 150, rand(rng) * 100 - 50
            r = 0.3 + rand(rng) * 3
            P = mk(8, lon0, lat0, r)
            Q = mk(8, lon0 + (rand(rng) - 0.5) * r, lat0 + (rand(rng) - 0.5) * r, r)
            ours = as_mpoly(fh(:intersection, ALG_S, P, Q))
            theirs = GO.intersection(sh, to_usp(P), to_usp(Q); target = GI.PolygonTrait())
            (ours === nothing || GI.npoint(theirs) < 4) && continue
            a1 = s2_area(write_wkb(ours))
            a2 = s2_area(write_wkb(to_ll(theirs)))
            max(a1, a2) == 0 && continue
            n += 1
            worst = max(worst, abs(a1 - a2) / max(a1, a2))
        end
        @test n > 40
        @test worst < 1e-9
    end
else
    @testset "s2geography oracle" begin
        @test_skip "S2Geography_jll unavailable on this platform"
    end
end
