# Tests for the OverlayNG phase-2b engine core (design §3): the labeller, the
# result builders (polygons with the real minimal-ring split + hole nesting,
# lines, points), and the internal `_overlay_ng` driver — end to end, over the
# phase-1 arrangement and the phase-2a graph.
#
# Equality strategy: planar results are checked against LibGEOS's own overlay
# (high-level API only) with GEOS topological `equals` (order/orientation/merge-
# granularity independent) plus `isValid`, which is independent of it — a
# self-intersecting ring can be `equals` to a valid one. Spherical is gated on
# area conservation. Ported JTS cases keep JTS's expected WKT answers, compared
# via GEOS `equals`.

using Test
include(joinpath(@__DIR__, "common.jl"))
import GeometryOps: Planar, Spherical, True

const OPS = OP_CODES

# Planar: check `_overlay_ng` against LibGEOS for one op.
#
# No area leg: GEOS `equals` is point-set equality, and equal point sets have
# equal areas under any correct `GO.area`, so an `isapprox` on the two areas
# could only fail where `equals` already had. The `isValid` leg is NOT
# dominated — a self-intersecting ring can be `equals` to a valid one.
function check_planar(op, A, B)
    r = GO._overlay_ng(Planar(), op, A, B; exact = EX)
    @test LG.isValid(lgc(r))
    @test LG.equals(lgc(r), geos_op(op, A, B))
    return r
end

check_all_ops(A, B) =
    for op in OPS
        @testset "$(opname(op))" begin check_planar(op, A, B) end
    end

# The rounded coordinate multiset of a geometry (for vertex-set equality).
function coordset(g; digits = 9)
    s = Set{Tuple{Float64, Float64}}()
    for p in GI.getpoint(g)
        push!(s, (round(GI.x(p); digits), round(GI.y(p); digits)))
    end
    return s
end

# ---------------------------------------------------------------------------
# 1. S2/S3 case suite — all four ops vs LibGEOS (planar)
# ---------------------------------------------------------------------------

@testset "overlapping squares (all ops)" begin
    A, B = SQ_A, SQ_B
    #-- the analytic areas (1 / 7 / 3 / 6) live in `api.jl`, where they run
    #-- through `@testset_implementations` across all four geometry backends
    check_all_ops(A, B)
    #-- vertex-set equality against GEOS (intersection = the overlap square)
    ri = GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, A, B; exact = EX)
    @test coordset(ri) == coordset(GO.tuples(geos_op(GO.OVERLAY_INTERSECTION, A, B)))
end

@testset "polygon-with-hole, B overlaps into the hole (§2.7 regression)" begin
    #-- the wrong-area-hole case: B reaches into A's hole. The material-interior
    #-- authority (§2.7) must give the hole the right side, or the areas invert.
    A = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
                    [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)]])
    B = GI.Polygon([[(5.0, 5.0), (12.0, 5.0), (12.0, 12.0), (5.0, 12.0), (5.0, 5.0)]])
    #-- A area = 100 - 16 = 84; overlap of B with A-material = B∩A = 21, pinned
    #-- across all four backends in `api.jl`
    check_all_ops(A, B)
end

@testset "collinear shared boundary (degenerate intersection, merged union)" begin
    A = SQ_A
    B = GI.Polygon([[(2.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 0.0)]])
    #-- intersection is the shared boundary line (1-D), not an area
    ri = GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, A, B; exact = EX)
    @test GI.trait(ri) isa GI.LineStringTrait
    @test LG.equals(lgc(ri), geos_op(GO.OVERLAY_INTERSECTION, A, B))
    #-- union merges into one 2×4 box
    ru = GO._overlay_ng(Planar(), GO.OVERLAY_UNION, A, B; exact = EX)
    @test isapprox(GO.area(ru), 8.0; rtol = 1e-12)
    @test LG.isValid(lgc(ru)) && LG.equals(lgc(ru), geos_op(GO.OVERLAY_UNION, A, B))
end

@testset "degree-6 coincident crossing (all ops)" begin
    #-- two A squares touching at (1,1) + a B triangle edge through (1,1)
    A = GI.MultiPolygon([[[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]],
                         [[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)]]])
    B = GI.Polygon([[(0.0, 0.0), (2.0, 2.0), (2.0, 0.0), (0.0, 0.0)]])
    #-- that the shared vertex really becomes one degree-6 graph node is asserted
    #-- in `overlay_graph.jl`, on the graph layer that owns it; here the fixture's
    #-- job is to drive all four ops end to end through that node
    check_all_ops(A, B)
end

#-- shape breadth: each row is one input shape run through all four ops against
#-- GEOS, and nothing distinguishes them beyond the fixture
@testset "shape breadth (all ops)" begin
    for (name, A, B) in (
        ("concave L-shapes",
         GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (1.0, 1.0), (1.0, 3.0), (0.0, 3.0), (0.0, 0.0)]]),
         GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (2.0, 3.0), (2.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])),
        ("MultiPolygon input",
         GI.MultiPolygon([[[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]],
                          [[(5.0, 5.0), (7.0, 5.0), (7.0, 7.0), (5.0, 7.0), (5.0, 5.0)]]]),
         GI.Polygon([[(1.0, 1.0), (6.0, 1.0), (6.0, 6.0), (1.0, 6.0), (1.0, 1.0)]])),
    )
        @testset "$name" begin check_all_ops(A, B) end
    end
end

# ---------------------------------------------------------------------------
# 2. Ported JTS OverlayNGTest subset (floating-safe area/line cases)
# ---------------------------------------------------------------------------
#
# Equality: GEOS `equals` against JTS's expected WKT (topological — handles
# JTS's coordinate ordering and OverlayNG's line-merging, which this engine
# emits as raw noded segments).
#
# SKIPPED, with reasons (all fixed-precision / topology-collapse tests that do
# not apply to a floating-only, non-snapping engine — design §0):
#   testTriangleFillingHoleUnion(Prec10), testBoxTri{Intersection,Union},
#   test2spikes{Intersection,Union}, testTriBoxIntersection,
#   testCollapse* / testSnapBoxGore* (topology collapse from precision rounding),
#   testVerySmallBIntersection (scale 1e8), testEdgeDisappears (scale 1e6),
#   testBcollapse* / testBNearVertexSnappingCausesInversion /
#   testBCollapsedHoleEdgeLabelledExterior (snap-rounding collapse),
#   testDisjointLinesRoundedIntersection (coordinate rounding to a point).
# `testTouchingPolyDifference` used to be skipped here as a substrate limitation
# ("the substrate does not self-node one input"). That stopped being true when
# `collect.jl` gained per-input vertex self-noding — the same fix `fuzz.jl`
# records as closing its SELF-TOUCHING INPUT defect class — and the case has been
# passing, un-run, ever since. It is in the table below.

const JTS_CASES = [
    # (name, op, A_wkt, B_wkt, expected_wkt)
    ("NestedShellsIntersection", GO.OVERLAY_INTERSECTION,
     "POLYGON ((100 200, 200 200, 200 100, 100 100, 100 200))",
     "POLYGON ((120 180, 180 180, 180 120, 120 120, 120 180))",
     "POLYGON ((120 180, 180 180, 180 120, 120 120, 120 180))"),
    ("NestedShellsUnion", GO.OVERLAY_UNION,
     "POLYGON ((100 200, 200 200, 200 100, 100 100, 100 200))",
     "POLYGON ((120 180, 180 180, 180 120, 120 120, 120 180))",
     "POLYGON ((100 200, 200 200, 200 100, 100 100, 100 200))"),
    ("AdjacentBoxesIntersection", GO.OVERLAY_INTERSECTION,
     "POLYGON ((100 200, 200 200, 200 100, 100 100, 100 200))",
     "POLYGON ((300 200, 300 100, 200 100, 200 200, 300 200))",
     "LINESTRING (200 100, 200 200)"),
    ("AdjacentBoxesUnion", GO.OVERLAY_UNION,
     "POLYGON ((100 200, 200 200, 200 100, 100 100, 100 200))",
     "POLYGON ((300 200, 300 100, 200 100, 200 200, 300 200))",
     "POLYGON ((100 100, 100 200, 200 200, 300 200, 300 100, 200 100, 100 100))"),
    ("TouchingHoleUnion", GO.OVERLAY_UNION,
     "POLYGON ((100 300, 300 300, 300 100, 100 100, 100 300), (200 200, 150 200, 200 300, 200 200))",
     "POLYGON ((130 160, 260 160, 260 120, 130 120, 130 160))",
     "POLYGON ((100 100, 100 300, 200 300, 300 300, 300 100, 100 100), (150 200, 200 200, 200 300, 150 200))"),
    ("TouchingMultiHoleUnion", GO.OVERLAY_UNION,
     "POLYGON ((100 300, 300 300, 300 100, 100 100, 100 300), (200 200, 150 200, 200 300, 200 200), (250 230, 216 236, 250 300, 250 230), (235 198, 300 200, 237 175, 235 198))",
     "POLYGON ((130 160, 260 160, 260 120, 130 120, 130 160))",
     "POLYGON ((100 300, 200 300, 250 300, 300 300, 300 200, 300 100, 100 100, 100 300), (200 300, 150 200, 200 200, 200 300), (250 300, 216 236, 250 230, 250 300), (300 200, 235 198, 237 175, 300 200))"),
    ("ATouchingNestedPolyUnion", GO.OVERLAY_UNION,
     "MULTIPOLYGON (((0 200, 200 200, 200 0, 0 0, 0 200), (50 50, 190 50, 50 200, 50 50)), ((60 100, 100 60, 50 50, 60 100)))",
     "POLYGON ((135 176, 180 176, 180 130, 135 130, 135 176))",
     "MULTIPOLYGON (((0 0, 0 200, 50 200, 200 200, 200 0, 0 0), (50 50, 190 50, 50 200, 50 50)), ((50 50, 60 100, 100 60, 50 50)))"),
    ("BoxLineIntersection", GO.OVERLAY_INTERSECTION,
     "POLYGON ((100 200, 200 200, 200 100, 100 100, 100 200))",
     "LINESTRING (50 150, 150 150)",
     "LINESTRING (100 150, 150 150)"),
    ("BoxLineUnion", GO.OVERLAY_UNION,
     "POLYGON ((100 200, 200 200, 200 100, 100 100, 100 200))",
     "LINESTRING (50 150, 150 150)",
     "GEOMETRYCOLLECTION (POLYGON ((200 200, 200 100, 100 100, 100 150, 100 200, 200 200)), LINESTRING (50 150, 100 150))"),
    ("LinePolygonUnion", GO.OVERLAY_UNION,
     "LINESTRING (50 150, 150 150)",
     "POLYGON ((100 200, 200 200, 200 100, 100 100, 100 200))",
     "GEOMETRYCOLLECTION (LINESTRING (50 150, 100 150), POLYGON ((100 200, 200 200, 200 100, 100 100, 100 150, 100 200)))"),
    ("LinePolygonUnionAlongPolyBoundary", GO.OVERLAY_UNION,
     "LINESTRING (150 300, 250 300)",
     "POLYGON ((100 400, 200 400, 200 300, 100 300, 100 400))",
     "GEOMETRYCOLLECTION (LINESTRING (200 300, 250 300), POLYGON ((200 300, 150 300, 100 300, 100 400, 200 400, 200 300)))"),
    ("LinePolygonIntersectionAlongPolyBoundary", GO.OVERLAY_INTERSECTION,
     "LINESTRING (150 300, 250 300)",
     "POLYGON ((100 400, 200 400, 200 300, 100 300, 100 400))",
     "LINESTRING (200 300, 150 300)"),
    ("PolygonLineVerticalIntersection", GO.OVERLAY_INTERSECTION,
     "POLYGON ((-200 -200, 200 -200, 200 200, -200 200, -200 -200))",
     "LINESTRING (-100 100, -100 -100)",
     "LINESTRING (-100 100, -100 -100)"),
    ("PolygonLineHorizontalIntersection", GO.OVERLAY_INTERSECTION,
     "POLYGON ((10 90, 90 90, 90 10, 10 10, 10 90))",
     "LINESTRING (20 50, 80 50)",
     "LINESTRING (20 50, 80 50)"),
    ("PolygonMultiLineUnion", GO.OVERLAY_UNION,
     "POLYGON ((100 200, 200 200, 200 100, 100 100, 100 200))",
     "MULTILINESTRING ((150 250, 150 50), (250 250, 250 50))",
     "GEOMETRYCOLLECTION (LINESTRING (150 50, 150 100), LINESTRING (150 200, 150 250), LINESTRING (250 50, 250 250), POLYGON ((100 100, 100 200, 150 200, 200 200, 200 100, 150 100, 100 100)))"),
    ("PolygonLineIntersectionOrder", GO.OVERLAY_INTERSECTION,
     "POLYGON ((1 1, 1 9, 9 9, 9 7, 3 7, 3 3, 9 3, 9 1, 1 1))",
     "MULTILINESTRING ((2 10, 2 0), (4 10, 4 0))",
     "MULTILINESTRING ((2 9, 2 1), (4 9, 4 7), (4 3, 4 1))"),
    ("AreaLineIntersection", GO.OVERLAY_INTERSECTION,
     "POLYGON ((360 200, 220 200, 220 180, 300 180, 300 160, 300 140, 360 200))",
     "MULTIPOLYGON (((280 180, 280 160, 300 160, 300 180, 280 180)), ((220 230, 240 230, 240 180, 220 180, 220 230)))",
     "GEOMETRYCOLLECTION (LINESTRING (280 180, 300 180), LINESTRING (300 160, 300 180), POLYGON ((220 180, 220 200, 240 200, 240 180, 220 180)))"),
    ("LineUnion", GO.OVERLAY_UNION,
     "LINESTRING (0 0, 1 1)", "LINESTRING (1 1, 2 2)",
     "MULTILINESTRING ((0 0, 1 1), (1 1, 2 2))"),
    ("Line2Union", GO.OVERLAY_UNION,
     "LINESTRING (0 0, 1 1, 0 1)", "LINESTRING (1 1, 2 2, 3 3)",
     "MULTILINESTRING ((0 0, 1 1), (0 1, 1 1), (1 1, 2 2, 3 3))"),
    ("Line3Union", GO.OVERLAY_UNION,
     "MULTILINESTRING ((0 1, 1 1), (2 2, 2 0))", "LINESTRING (0 0, 1 1, 2 2, 3 3)",
     "MULTILINESTRING ((0 0, 1 1), (0 1, 1 1), (1 1, 2 2), (2 0, 2 2), (2 2, 3 3))"),
    ("Line4Union", GO.OVERLAY_UNION,
     "LINESTRING (100 300, 200 300, 200 100, 100 100)",
     "LINESTRING (300 300, 200 300, 200 300, 200 100, 300 100)",
     "MULTILINESTRING ((200 100, 100 100), (300 300, 200 300), (200 300, 200 100), (200 100, 300 100), (100 300, 200 300))"),
    ("LineFigure8Union", GO.OVERLAY_UNION,
     "LINESTRING (5 1, 2 2, 5 3, 2 4, 5 5)", "LINESTRING (5 1, 8 2, 5 3, 8 4, 5 5)",
     "MULTILINESTRING ((5 1, 2 2, 5 3), (5 1, 8 2, 5 3), (5 3, 2 4, 5 5), (5 3, 8 4, 5 5))"),
    ("LineRingUnion", GO.OVERLAY_UNION,
     "LINESTRING (1 1, 5 5, 9 1)", "LINESTRING (1 1, 9 1)",
     "MULTILINESTRING ((1 1, 5 5, 9 1), (1 1, 9 1))"),
    ("PolygonFlatCollapseIntersection", GO.OVERLAY_INTERSECTION,
     "POLYGON ((200 100, 150 200, 250 200, 150 200, 100 100, 200 100))",
     "POLYGON ((50 150, 250 150, 250 50, 50 50, 50 150))",
     "POLYGON ((175 150, 200 100, 100 100, 125 150, 175 150))"),
    #-- both of these were absent from the ported set and both pass. The first
    #-- was skipped for a substrate limitation that no longer exists (above); the
    #-- second is the only JTS case whose result is a genuinely mixed collection,
    #-- exercising all three result builders in one op.
    ("TouchingPolyDifference", GO.OVERLAY_DIFFERENCE,
     "POLYGON ((200 200, 200 0, 0 0, 0 200, 200 200), (100 100, 50 100, 50 200, 100 100))",
     "POLYGON ((150 100, 100 100, 150 200, 150 100))",
     "MULTIPOLYGON (((0 0, 0 200, 50 200, 50 100, 100 100, 150 100, 150 200, 200 200, 200 0, 0 0)), ((50 200, 150 200, 100 100, 50 200)))"),
    ("AreaLinePointIntersection", GO.OVERLAY_INTERSECTION,
     "POLYGON ((100 100, 200 100, 200 150, 250 100, 300 100, 300 150, 350 100, 350 200, 100 200, 100 100))",
     "POLYGON ((100 140, 170 140, 200 100, 400 100, 400 30, 100 30, 100 140))",
     "GEOMETRYCOLLECTION (POINT (350 100), LINESTRING (250 100, 300 100), POLYGON ((100 100, 100 140, 170 140, 200 100, 100 100)))"),
]

@testset "ported JTS OverlayNGTest subset ($(length(JTS_CASES)) cases)" begin
    for (name, op, awkt, bwkt, ewkt) in JTS_CASES
        @testset "$name" begin
            r = GO._overlay_ng(Planar(), op, giwkt(awkt), giwkt(bwkt); exact = EX)
            @test LG.isValid(lgc(r))
            @test LG.equals(lgc(r), LG.readgeom(ewkt))
        end
    end
end

# ---------------------------------------------------------------------------
# 3. Ring-builder specifics
# ---------------------------------------------------------------------------

@testset "self-touching result ring (minimal-ring split)" begin
    #-- A minus a triangle B that touches A's boundary at a single point (0,3):
    #-- the result ring self-touches there and must split into shell + hole.
    A = GI.Polygon([[(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (0.0, 6.0), (0.0, 0.0)]])
    B = GI.Polygon([[(0.0, 3.0), (4.0, 1.0), (4.0, 5.0), (0.0, 3.0)]])
    r = GO._overlay_ng(Planar(), GO.OVERLAY_DIFFERENCE, A, B; exact = EX)
    @test GI.trait(r) isa GI.PolygonTrait
    @test GI.nring(r) == 2                       # shell + one split-out hole
    @test LG.isValid(lgc(r))
    @test LG.equals(lgc(r), geos_op(GO.OVERLAY_DIFFERENCE, A, B))
    @test isapprox(GO.area(r), 28.0; rtol = 1e-12)  # 36 - 8
end

@testset "free-hole assignment (strictly interior hole)" begin
    #-- difference of a strictly-interior square from a shell → one free hole
    Big = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
    Inner = GI.Polygon([[(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)]])
    r = GO._overlay_ng(Planar(), GO.OVERLAY_DIFFERENCE, Big, Inner; exact = EX)
    @test GI.trait(r) isa GI.PolygonTrait
    @test GI.nring(r) == 2
    @test LG.isValid(lgc(r))
    @test isapprox(GO.area(r), 84.0; rtol = 1e-12)
    #-- ring count and area together still permit the hole being in the wrong
    #-- place; the oracle is what pins where it went
    @test LG.equals(lgc(r), geos_op(GO.OVERLAY_DIFFERENCE, Big, Inner))
end

@testset "union of a multi-island geometry (France-class nesting)" begin
    #-- 20 disjoint island squares, unioned with a copy shifted onto their
    #-- diagonal neighbours. The result is ONE connected shell with 24 free
    #-- holes — the diagonal overlap chains every island together and the gaps
    #-- between them become the holes — so this is the many-shell *free-hole
    #-- assignment* case, not the many-shell case the comment used to claim.
    isl = Vector{Vector{Vector{Tuple{Float64, Float64}}}}()
    for i in 0:19
        x = (i % 5) * 3.0; y = (i ÷ 5) * 3.0
        push!(isl, [[(x, y), (x + 2, y), (x + 2, y + 2), (x, y + 2), (x, y)]])
    end
    MI = GI.MultiPolygon(isl)
    #-- shift by (1,1) so islands overlap their diagonal neighbours
    MI2 = GO.apply(GI.PointTrait(), MI) do p
        (GI.x(p) + 1.0, GI.y(p) + 1.0)
    end
    r = GO._overlay_ng(Planar(), GO.OVERLAY_UNION, MI, MI2; exact = EX)
    geos = geos_op(GO.OVERLAY_UNION, MI, MI2)
    @test LG.isValid(lgc(r))
    @test LG.equals(lgc(r), geos)
    @test isapprox(GO.area(r), GO.area(GO.tuples(geos)); rtol = 1e-12)
end

# ---------------------------------------------------------------------------
# 4. Spherical end-to-end
# ---------------------------------------------------------------------------

# The three spherical area-conservation identities on the (0,0)-(20,20) /
# (10,10)-(30,30) quad pair live in `api.jl`, on the same fixture, run through
# the public surface and with two extra legs this copy did not have (all four
# results strictly positive, and the spherical answer differing from the planar
# one). Kept there rather than here so there is exactly one copy — and kept in
# FULL there: trimming it to a single identity while deleting this block would
# have taken two of the three out of the repository altogether.

#=
Emission used to pick between the two antipodal candidates `±(na×nb)` for a
crossing node by evaluating `_strictly_in_arc3` in Float64. Those determinants
vanish as the crossing approaches an endpoint of either arc, so on inputs whose
vertices are ulps apart — here, one region decomposed into two components two
different ways, as an overlay result fed back in — a crossing node was emitted at
its ANTIPODE. The result kept the right combinatorial structure and the right
shell/hole roles; one of its vertices was simply half a sphere away, which made
`intersection` a clean quarter of the sphere and `symdifference` a clean half.

X and Y below denote the same region, so intersection and union are that region
and both differences are empty. Planar gets this right on the same input and is
checked alongside, since the sign choice is a spherical-only construction.
=#
@testset "spherical crossing emitted at its antipode (near-endpoint crossing)" begin
    X = GO.tuples(LG.readgeom("MULTIPOLYGON (((0 0, 0 1, 0.5 1.0000380706528735, 0.5 0.5, 0.9999999999999998 0.5000190382262165, 1 0, 0 0)), ((0.9999999999999998 1, 0.5 1.0000380706528735, 0.5 1.5000000000000004, 1.5000000000000002 1.5000000000000004, 1.5000000000000002 0.5000000000000001, 0.9999999999999998 0.5000190382262165, 0.9999999999999998 1)))"))
    Y = GO.tuples(LG.readgeom("MULTIPOLYGON (((0.5 1.0000380706528735, 0 1, 0 0, 1 0, 1 0.5000190382262163, 0.5 0.5, 0.5 1.0000380706528735)), ((0.5 1.0000380706528735, 0.9999999999999998 1, 1 0.5000190382262163, 1.5000000000000002 0.5, 1.5000000000000002 1.5000000000000002, 0.5 1.5000000000000002, 0.5 1.0000380706528735)))"))

    for m in (Planar(), Spherical())
        aX, aY = GO.area(m, X), GO.area(m, Y)
        #-- fixture premise: X and Y denote the same region (no engine involved)
        @test isapprox(aX, aY; rtol = 1e-15)
        ops = Dict(op => GO.area(m, GO._overlay_ng(m, op, X, Y; exact = EX)) for op in
                   (GO.OVERLAY_INTERSECTION, GO.OVERLAY_UNION,
                    GO.OVERLAY_DIFFERENCE, GO.OVERLAY_SYMDIFFERENCE))
        @test isapprox(ops[GO.OVERLAY_INTERSECTION], aX; rtol = 1e-12)
        @test isapprox(ops[GO.OVERLAY_UNION], aX; rtol = 1e-12)
        @test ops[GO.OVERLAY_DIFFERENCE] <= 1e-12 * aX
        @test ops[GO.OVERLAY_SYMDIFFERENCE] <= 1e-12 * aX
    end

    #-- and directly on the defect: no emitted vertex may be a hemisphere away
    #-- from the crossing's own arcs. Both emitted coordinates sit in [0, 2].
    r = GO._overlay_ng(Spherical(), GO.OVERLAY_INTERSECTION, X, Y; exact = EX)
    @test all(p -> 0 <= GI.x(p) <= 2 && 0 <= GI.y(p) <= 2, GI.getpoint(r))

    #-- and at the emitter itself, on the offending crossing: X's
    #-- (1.5, 0.5)->(1, 0.50001904) against Y's (1, 1)->(1, 0.50001904), whose
    #-- endpoints differ in the last ulp of latitude so the crossing sits ~1e-16
    #-- rad from both. Float `_strictly_in_arc3` rejected BOTH candidates here
    #-- and the old code fell through to `-d`. (The exact-arithmetic half of this
    #-- lives with the predicate, in the spherical kernel conformance suite.)
    m = Spherical()
    a0, a1, b0, b1 = (GO._to_kernel_point(m, p) for p in
                      ((1.5000000000000002, 0.5000000000000001), (0.9999999999999998, 0.5000190382262165),
                       (0.9999999999999998, 1.0), (1.0, 0.5000190382262163)))
    lon, lat = GO._emit_node_coord(GO.crossing_node(a0, a1, b0, b1))
    @test 0.99 <= lon <= 1.01 && 0.4 <= lat <= 0.6
end

@testset "spherical empty-vs-full disambiguation (§3 amendment 6)" begin
    A = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
    Bdisjoint = GI.Polygon([[(40.0, 0.0), (50.0, 0.0), (50.0, 10.0), (40.0, 10.0), (40.0, 0.0)]])
    #-- disjoint intersection → empty (NOT full-sphere)
    ri = GO._overlay_ng(Spherical(), GO.OVERLAY_INTERSECTION, A, Bdisjoint; exact = EX)
    @test GI.npoint(ri) == 0
    @test isapprox(GO.area(Spherical(), ri), 0.0; atol = 1e-9)

    #-- the disambiguation function directly: a boundaryless union that covers
    #-- everything is the full sphere and must throw; an empty intersection must not.
    #-- the throw is conjunction-gated on `_covers_everything`, so asserting the
    #-- throw asserts the predicate too; only the throws are checked here
    inp_union = GO._OverlayInput(Spherical(), A, A, 2, 2, EX, false, false, nothing, nothing)
    @test_throws ArgumentError GO._resolve_empty_result(Spherical(), GO.OVERLAY_UNION, inp_union)

    inp_int = GO._OverlayInput(Spherical(), A, Bdisjoint, 2, 2, EX, false, false, nothing, nothing)
    @test GI.npoint(GO._resolve_empty_result(Spherical(), GO.OVERLAY_INTERSECTION, inp_int)) == 0
end

# ---------------------------------------------------------------------------
# 5. Input validation + empty inputs
# ---------------------------------------------------------------------------

@testset "input validation and empty short-circuits" begin
    A = SQ_A
    #-- point inputs are routed to the phase-3 point builders (they used to be
    #-- rejected here); see overlay_points.jl for their full coverage
    @test GI.trait(GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, GI.Point((1.0, 1.0)), A;
                                  exact = EX)) isa GI.PointTrait
    #-- geometry collections are still rejected
    @test_throws ArgumentError GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION,
        GI.GeometryCollection([GI.Point((1.0, 1.0))]), A; exact = EX)
    #-- a disjoint intersection is empty. This does NOT show that the planar
    #-- envelope short circuit ran: with `_env_disjoint` forced false the full
    #-- pipeline reaches the same empty answer, so the assertion is about the
    #-- ANSWER, and the short circuit itself is a performance property.
    Far = GI.Polygon([[(100.0, 100.0), (102.0, 100.0), (102.0, 102.0), (100.0, 102.0), (100.0, 100.0)]])
    r = GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, A, Far; exact = EX)
    @test GI.npoint(r) == 0
end

# ---------------------------------------------------------------------------
# 5b. Clip-envelope pruning (the construct-free `RingClipper`)
# ---------------------------------------------------------------------------
#
# `_overlay_ng` prunes input SEGMENTS whose bbox misses a per-op clip box before
# they ever become `NodedEdge`s (`OVERLAY_INTERSECTION`: both sides to
# `env(A) ∩ env(B)`; `OVERLAY_DIFFERENCE`: B only, to `env(A)`), leaving the
# surviving chains open in the graph. The correctness claim is that a pruned
# edge's location is recovered by point-in-area against the ORIGINAL input, so
# the cases below are chosen to prune a WHOLE side away and then demand the right
# answer — which only the original-geometry fallback can supply.
#
# Every one of them asserts that pruning actually fired, by counting the noded
# edges each side contributed; without that an "A ∩ hugeB == A" test is satisfied
# by the unpruned pipeline too and proves nothing about this path.

#-- the driver's own envelope rule, so these prune exactly as `_overlay_ng` does
clip_for(op, A, B) = GO._overlay_clip_envelopes(op, GI.extent(A), GI.extent(B))

function clipped_arr(op, A, B)
    ca, cb = clip_for(op, A, B)
    return GO.NodedArrangement(Planar(), A, B; exact = EX, clip_a = ca, clip_b = cb)
end

#-- noded edges contributed by one side (`is_a = true` → A)
side_edges(arr, is_a) = count(e -> arr.segstrings[e.string_idx].is_a == is_a, arr.edges)

#-- the same pipeline with pruning switched OFF, for a like-for-like comparison
#-- (area × area only, which is all this section builds)
function overlay_unpruned(op, A, B)
    arr = GO.NodedArrangement(Planar(), A, B; exact = EX)
    g = GO.OverlayGraph(Planar(), arr; exact = EX)
    input = GO._OverlayInput(Planar(), A, B, 2, 2, EX, false, false, nothing, nothing)
    GO._compute_labelling!(g, input)
    GO._mark_result_area_edges!(g, op)
    GO._unmark_duplicate_edges_from_result_area!(g)
    return GO._extract_result(Planar(), op, g, input, nothing; exact = EX)
end

const CLIP_A_UNIT = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
const CLIP_B_HUGE = GI.Polygon([[(-100.0, -100.0), (100.0, -100.0), (100.0, 100.0),
                                 (-100.0, 100.0), (-100.0, -100.0)]])
#-- the same huge square with a hole that swallows the unit square whole
const CLIP_B_HOLE = GI.Polygon([[(-100.0, -100.0), (100.0, -100.0), (100.0, 100.0),
                                 (-100.0, 100.0), (-100.0, -100.0)],
                                [(-50.0, -50.0), (-50.0, 50.0), (50.0, 50.0),
                                 (50.0, -50.0), (-50.0, -50.0)]])

@testset "A strictly inside a huge B — every B edge prunes away" begin
    A, B = CLIP_A_UNIT, CLIP_B_HUGE
    arr = clipped_arr(GO.OVERLAY_INTERSECTION, A, B)
    #-- B's shell never enters env(A), so the arrangement holds NO B linework:
    #-- the answer below can only come from locating A's edges against the
    #-- original B (labeller pass 5), which is the whole point of the spike
    @test side_edges(arr, false) == 0
    @test side_edges(arr, true) == side_edges(GO.NodedArrangement(Planar(), A, B; exact = EX), true)

    ri = GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, A, B; exact = EX)
    @test LG.isValid(lgc(ri))
    @test LG.equals(lgc(ri), lgc(A))
    @test isapprox(GO.area(ri), 1.0; rtol = 1e-12)

    #-- union is NOT pruned (there is no box the result stays inside), so it is
    #-- unaffected — and it would be wrong if the intersection box leaked into it
    ru = GO._overlay_ng(Planar(), GO.OVERLAY_UNION, A, B; exact = EX)
    @test LG.equals(lgc(ru), lgc(B))
    @test isapprox(GO.area(ru), 40000.0; rtol = 1e-12)
    @test GO.num_edges(clipped_arr(GO.OVERLAY_UNION, A, B)) ==
          GO.num_edges(GO.NodedArrangement(Planar(), A, B; exact = EX))
end

@testset "A inside B's hole — intersection empty, difference untouched" begin
    A, B = CLIP_A_UNIT, CLIP_B_HOLE
    arr = clipped_arr(GO.OVERLAY_INTERSECTION, A, B)
    #-- neither B ring comes near env(A): the whole of B prunes, and the hole's
    #-- EXTERIOR verdict likewise comes from the original geometry
    @test side_edges(arr, false) == 0
    ri = GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, A, B; exact = EX)
    @test GI.npoint(ri) == 0

    #-- A ∖ B is all of A, and here too B prunes to nothing (to `env(A)`)
    arrd = clipped_arr(GO.OVERLAY_DIFFERENCE, A, B)
    @test side_edges(arrd, false) == 0
    rd = GO._overlay_ng(Planar(), GO.OVERLAY_DIFFERENCE, A, B; exact = EX)
    @test LG.isValid(lgc(rd))
    @test LG.equals(lgc(rd), lgc(A))
    @test isapprox(GO.area(rd), 1.0; rtol = 1e-12)
end

@testset "difference against a far-away B — B prunes to nothing, A survives" begin
    A = SQ_A
    B = GI.Polygon([[(100.0, 100.0), (102.0, 100.0), (102.0, 102.0), (100.0, 102.0), (100.0, 100.0)]])
    #-- DIFFERENCE never prunes A, and prunes B to `env(A)` — which B misses
    ca, cb = clip_for(GO.OVERLAY_DIFFERENCE, A, B)
    @test ca === nothing
    arr = clipped_arr(GO.OVERLAY_DIFFERENCE, A, B)
    @test side_edges(arr, false) == 0
    @test side_edges(arr, true) == 4
    r = GO._overlay_ng(Planar(), GO.OVERLAY_DIFFERENCE, A, B; exact = EX)
    @test LG.equals(lgc(r), lgc(A))
    @test isapprox(GO.area(r), 4.0; rtol = 1e-12)
end

@testset "MultiPolygon with a far component — only that component prunes" begin
    A = SQ_A
    #-- one component straddles A, one sits far away; env(B) spans both, so the
    #-- intersection box is env(A) ∩ env(B) and only the far component leaves it
    B = GI.MultiPolygon([[[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]],
                         [[(50.0, 50.0), (52.0, 50.0), (52.0, 52.0), (50.0, 52.0), (50.0, 50.0)]]])
    plain = GO.NodedArrangement(Planar(), A, B; exact = EX)
    arr = clipped_arr(GO.OVERLAY_INTERSECTION, A, B)
    @test side_edges(arr, false) < side_edges(plain, false)   # pruning fired ...
    @test side_edges(arr, false) > 0                          # ... but not to nothing
    for op in OPS
        r = GO._overlay_ng(Planar(), op, A, B; exact = EX)
        @test LG.isValid(lgc(r))
        @test LG.equals(lgc(r), geos_op(op, A, B))
        @test LG.equals(lgc(r), lgc(overlay_unpruned(op, A, B)))
    end
end

@testset "a truncated star outside the box does not raise a topology error" begin
    #-- FOUR valid B components meeting at one apex far outside env(A). Three of
    #-- them reach back to A with both of their apex edges (which therefore
    #-- survive); the third component's SHORT apex edge (100,5)->(99,4) has a bbox
    #-- that misses the box, so exactly one ray of its wedge is pruned and the
    #-- apex star loses an edge from the MIDDLE of its CCW order.
    #--
    #-- That is the shape `_propagate_area_locations!` cannot read: it walks the
    #-- star as a closed cycle, so the half-open wedge reads as a side-location
    #-- conflict. Verified: running pass 1 over this graph WITHOUT the
    #-- truncated-node skip in `_compute_labelling!` raises
    #-- `_OverlayTopologyError("side location conflict: arg 1")`. (The middle of
    #-- the order matters — pass 1 never re-checks its own start edge, so a wedge
    #-- broken at the star's wrap-around point is silently tolerated.)
    A = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
    apex = (100.0, 5.0)
    B = GI.MultiPolygon([
        [[(0.0, 1.0), apex, (0.0, 3.0), (0.0, 1.0)]],
        [[(0.0, 7.0), apex, (0.0, 9.0), (0.0, 7.0)]],
        [[apex, (99.0, 4.0), (0.0, 0.5), apex]],
        [[apex, (10.0, -100.0), (9.0, -100.0), apex]]])
    @test LG.isValid(lgc(B))
    arr = clipped_arr(GO.OVERLAY_INTERSECTION, A, B)
    @test side_edges(arr, false) < side_edges(GO.NodedArrangement(Planar(), A, B; exact = EX), false)
    #-- the apex is exactly the node whose star lost an edge
    @test count(arr.truncated) > 0
    for op in OPS
        r = GO._overlay_ng(Planar(), op, A, B; exact = EX)
        @test LG.isValid(lgc(r))
        @test LG.equals(lgc(r), geos_op(op, A, B))
    end
end

# ## Self-noding under the clip box
#
# The clip box also restricts BOTH self-noding passes (collect.jl), which is
# where the time actually was. The cases below are the adversarial ones: a
# same-side incidence between a KEPT segment and a PRUNED one, outside the box,
# so the kept segment is deliberately left unsplit at a point where the input
# touches itself. Each asserts the drop really happened (`self_nodes`) and then
# demands the same answer as both the unpruned pipeline and GEOS.

#-- interior nodes the self-noding passes record for one side, with and without
#-- the clip box. This is the direct measurement of "the pruned path fired".
function self_nodes(A, B, op)
    P = Planar()
    ssa = GO._overlay_segstrings(P, A, true; exact = EX)
    ssb = GO._overlay_segstrings(P, B, false; exact = EX)
    ca, cb = clip_for(op, A, B)
    ta = GO._relate_edge_index(P, ssa)
    tb = GO._relate_edge_index(P, ssb)
    function census(ss, off, tree, clip)
        t = GO.NodeTable{Tuple{Float64, Float64}}()
        sn = Dict{Tuple{Int32, Int32}, Vector{Int32}}()
        GO._collect_self_crossings!(P, t, sn, ss, Int32(off); exact = EX, clip)
        GO._collect_self_vertex_nodes!(P, t, sn, ss, Int32(off), tree; exact = EX, clip)
        return sum(length, values(sn); init = 0)
    end
    return (a_clipped = census(ssa, 0, ta, ca), a_plain = census(ssa, 0, ta, nothing),
            b_clipped = census(ssb, length(ssa), tb, cb),
            b_plain = census(ssb, length(ssa), tb, nothing))
end

function check_prune_agrees(A, B; ops = OPS)
    for op in ops
        r = GO._overlay_ng(Planar(), op, A, B; exact = EX)
        @test LG.isValid(lgc(r))
        @test LG.equals(lgc(r), geos_op(op, A, B))
    end
end

@testset "same-side vertex touch outside the box (two A components)" begin
    #-- A2's apex lies STRICTLY INSIDE A1's bottom edge at (15,0), well outside
    #-- the box. A1's bottom edge survives (its bbox spans into the box), A2's two
    #-- apex segments do not — so the self-vertex pass drops a node that the
    #-- unpruned pass finds, and A1's bottom edge is left unsplit there.
    A = GI.MultiPolygon([[[(0.0, 0.0), (20.0, 0.0), (20.0, 1.0), (0.0, 1.0), (0.0, 0.0)]],
                         [[(15.0, 0.0), (20.0, -5.0), (10.0, -5.0), (15.0, 0.0)]]])
    B = GI.Polygon([[(0.0, -1.0), (4.0, -1.0), (4.0, 2.0), (0.0, 2.0), (0.0, -1.0)]])
    @test LG.isValid(lgc(A))
    c = self_nodes(A, B, GO.OVERLAY_INTERSECTION)
    @test c.a_plain > 0                 # the touch really is a self-node ...
    @test c.a_clipped < c.a_plain       # ... and the clip box really drops it
    ri = GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, A, B; exact = EX)
    @test LG.equals(lgc(ri), lgc(overlay_unpruned(GO.OVERLAY_INTERSECTION, A, B)))
    @test isapprox(GO.area(ri), 4.0; rtol = 1e-12)   # the bar clipped to x ∈ [0,4]
    check_prune_agrees(A, B)
end

@testset "shell–hole vertex touch outside the box" begin
    #-- the hole touches the shell's bottom edge at (15,0), outside the box. The
    #-- shell edge is kept, all three hole segments prune.
    A = GI.Polygon([[(0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
                    [(15.0, 0.0), (12.0, 4.0), (18.0, 4.0), (15.0, 0.0)]])
    B = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
    @test LG.isValid(lgc(A))
    c = self_nodes(A, B, GO.OVERLAY_INTERSECTION)
    @test c.a_plain > 0
    @test c.a_clipped < c.a_plain
    ri = GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, A, B; exact = EX)
    @test LG.equals(lgc(ri), lgc(overlay_unpruned(GO.OVERLAY_INTERSECTION, A, B)))
    @test isapprox(GO.area(ri), 40.0; rtol = 1e-12)  # B is entirely inside A
    check_prune_agrees(A, B)
end

@testset "self-crossing LineString, crossing outside the box" begin
    #-- the linear all-pairs pass, not the vertex pass: A crosses ITSELF at
    #-- (15,5), outside the box, and crosses B inside it. The long horizontal
    #-- segment survives; the vertical one that produces the self-crossing does
    #-- not, so the crossing node is dropped and the survivor stays unsplit.
    A = GI.LineString([(0.0, 5.0), (20.0, 5.0), (20.0, 0.0), (15.0, 0.0), (15.0, 10.0)])
    B = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
    c = self_nodes(A, B, GO.OVERLAY_INTERSECTION)
    @test c.a_plain > 0
    @test c.a_clipped < c.a_plain
    ri = GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, A, B; exact = EX)
    @test LG.equals(lgc(ri), geos_op(GO.OVERLAY_INTERSECTION, A, B))
    @test isapprox(GO.perimeter(ri), 4.0; rtol = 1e-12)  # the (0,5)–(4,5) stub
    #-- every op, against GEOS (`overlay_unpruned` is area×area only)
    check_prune_agrees(A, B)
end

@testset "pruning is answer-preserving on the ordinary pairs" begin
    #-- the direct pruned-vs-unpruned comparison, on shapes whose linework
    #-- straddles the box rather than sitting wholly inside or outside it
    Ring = GI.Polygon([[(0.0, 0.0), (12.0, 0.0), (12.0, 12.0), (0.0, 12.0), (0.0, 0.0)],
                       [(3.0, 3.0), (3.0, 9.0), (9.0, 9.0), (9.0, 3.0), (3.0, 3.0)]])
    Bar = GI.Polygon([[(6.0, -5.0), (20.0, -5.0), (20.0, 6.0), (6.0, 6.0), (6.0, -5.0)]])
    pairs = ((SQ_A, SQ_B), (Ring, Bar), (Ring, CLIP_A_UNIT), (CLIP_A_UNIT, CLIP_B_HOLE))
    for (A, B) in pairs, op in OPS
        r = GO._overlay_ng(Planar(), op, A, B; exact = EX)
        @test LG.equals(lgc(r), lgc(overlay_unpruned(op, A, B)))
    end
end

# ---------------------------------------------------------------------------
# 6. Natural Earth shifted-self smoke, both manifolds (env-gated)
# ---------------------------------------------------------------------------

ne_ok = try
    import NaturalEarth, GeoJSON
    include(joinpath(@__DIR__, "..", "..", "..", "data", "natural_earth_pairs.jl"))
    global ne_names, ne_geoms = load_ne(110)
    length(ne_geoms) > 0
catch err
    @info "Natural Earth subset skipped (data unavailable)" err
    false
end

# ONE country, not four. `realdata_identities.jl` runs this identity over ten
# shifted-self cases at a tighter rtol, so breadth here buys nothing — except for
# Egypt, which is not subsumed: its shift is horizontal (`dy = 0`) and Egypt is
# the one pick with long exactly-horizontal border segments, so the shifted copy
# lands COLLINEAR with the original rather than crossing it. That is a different
# degeneracy from the crossing case the other three exercise.
@testset "Natural Earth shifted-self area conservation (spherical + planar)" begin
    if !ne_ok
        @test_skip "Natural Earth data unavailable"
    else
        idx = findfirst(==("Egypt"), ne_names)
        @test idx !== nothing
        if idx !== nothing
            A = ne_geoms[idx]
            B = shift_geom(A, 0.5, 0.0)
            for m in (Planar(), Spherical())
                ri = GO._overlay_ng(m, GO.OVERLAY_INTERSECTION, A, B; exact = EX)
                ru = GO._overlay_ng(m, GO.OVERLAY_UNION, A, B; exact = EX)
                @test isapprox(GO.area(m, ru) + GO.area(m, ri),
                               GO.area(m, A) + GO.area(m, B); rtol = 1e-9)
            end
        end
    end
end

# ---------------------------------------------------------------------------
# 7. `target`: narrowing the result to one dimension
# ---------------------------------------------------------------------------

# The dimension of a geometry, by trait.
gdim(g) = (t = GI.trait(g);
    (t isa GI.PolygonTrait || t isa GI.MultiPolygonTrait) ? 2 :
    (t isa GI.LineStringTrait || t isa GI.LinearRingTrait ||
     t isa GI.MultiLineStringTrait) ? 1 : 0)

# The dimension-`dim` atomic components of an UNtargeted result, which is a
# single geometry, a Multi, or a GeometryCollection.
function untargeted_parts(r, dim)
    t = GI.trait(r)
    t isa GI.GeometryCollectionTrait &&
        return collect(Iterators.flatten(untargeted_parts(g, dim) for g in GI.getgeom(r)))
    gdim(r) == dim || return Any[]
    (t isa GI.MultiPolygonTrait || t isa GI.MultiLineStringTrait ||
     t isa GI.MultiPointTrait) && return collect(GI.getgeom(r))
    return Any[r]
end

# The atomic components of a TARGETED result, whichever container it came in.
targeted_parts(r) = r isa AbstractVector ? r : collect(GI.getgeom(r))

coordlist(v) = [GI.coordinates(g) for g in v]

const TARGETS = ((GI.PolygonTrait(),         GI.MultiPolygonTrait(),    2),
                 (GI.LineStringTrait(),      GI.MultiLineStringTrait(), 1),
                 (GI.PointTrait(),           GI.MultiPointTrait(),      0))

# A pair whose intersection is genuinely mixed-dimension, i.e. the one shape that
# reaches `_create_result_geometry`'s GeometryCollection branch. Getting there
# takes a multipolygon: one component must overlap A in an area while another
# meets it in a line only, and the two must stay disjoint for the input to be
# valid. (A single polygon that overlaps and shares boundary yields just the
# polygon — the shared segment is interior to the result area and is dropped.)
const MIXED_A = SQ_A
const MIXED_B = GI.MultiPolygon([
    GI.Polygon([[(1.0, 0.0), (3.0, 0.0), (3.0, 1.0), (1.0, 1.0), (1.0, 0.0)]]),   # overlaps
    GI.Polygon([[(-2.0, 0.0), (0.0, 0.0), (0.0, 2.0), (-2.0, 2.0), (-2.0, 0.0)]]) # touches
])

@testset "target is exactly the untargeted result, filtered by dimension" begin
    A  = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B  = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    #-- shares one edge only: the intersection is a pure line
    Cedge = GI.Polygon([[(2.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 0.0)]])
    Etouch = GI.Polygon([[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0), (2.0, 2.0)]])
    L  = GI.LineString([(-1.0, 1.0), (3.0, 1.0)])
    L2 = GI.LineString([(0.5, -1.0), (0.5, 3.0)])
    P  = GI.MultiPoint([(0.5, 0.5), (5.0, 5.0), (2.0, 2.0)])
    P2 = GI.MultiPoint([(0.5, 0.5), (9.0, 9.0)])

    #-- one pair per distinct branch signature (driver path × which dimensions are
    #-- non-empty, per op) — that tuple is everything the target machinery
    #-- branches on, so pairs sharing a signature share their coverage exactly.
    pairs = [("area × area", A, B), ("shared edge", A, Cedge),
             ("corner touch", A, Etouch), ("line × area", L, A), ("area × line", A, L),
             ("line × line", L, L2), ("point × area", P, A), ("area × point", A, P),
             ("point × line", P, L), ("point × point", P, P2),
             #-- the only pair here that produces a GeometryCollection, so the
             #-- only one exercising the mixed-dimension branch on both sides
             ("mixed dims", MIXED_A, MIXED_B), ("mixed dims swapped", MIXED_B, MIXED_A)]

    #-- guard the guard: if this stops being a collection the sweep silently
    #-- stops covering `untargeted_parts`' flattening branch
    @test GI.trait(GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, MIXED_A, MIXED_B;
                                  exact = EX)) isa GI.GeometryCollectionTrait

    for (name, X, Y) in pairs, op in OPS
        untargeted = GO._overlay_ng(Planar(), op, X, Y; exact = EX)
        for (single, multi, dim) in TARGETS
            want = coordlist(untargeted_parts(untargeted, dim))
            @test coordlist(targeted_parts(GO._overlay_ng(Planar(), op, X, Y;
                                           exact = EX, target = single))) == want
            @test coordlist(targeted_parts(GO._overlay_ng(Planar(), op, X, Y;
                                           exact = EX, target = multi))) == want
        end
    end

    #=
    The container SHAPE is a function of the target alone, so it is pinned once
    per target rather than once per (pair, op) — the assertion above cannot do it,
    because `coordlist` compares two empty lists as equal whatever container they
    arrived in. It is also the only ABSOLUTE pin of the shape: everything else in
    this file compares one result's type against another's. Three fixtures, so
    that each dimension is non-empty in at least one of them.
    =#
    for (single, multi, _) in TARGETS, (X, Y) in ((A, B), (A, Cedge), (A, Etouch))
        @test GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, X, Y;
                             exact = EX, target = single) isa AbstractVector
        @test GI.trait(GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, X, Y;
                                      exact = EX, target = multi)) === multi
    end
end

@testset "target fixes the return type, empty or not" begin
    A  = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B  = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    Cedge = GI.Polygon([[(2.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 0.0)]])
    Far = GI.Polygon([[(9.0, 9.0), (10.0, 9.0), (10.0, 10.0), (9.0, 10.0), (9.0, 9.0)]])
    alg = GO.OverlayNG()

    for (single, multi, dim) in TARGETS
        #-- a nonempty-at-this-dimension case and two empty ones (one from the
        #-- pipeline, one from the disjoint-envelope short circuit)
        nonempty = dim == 2 ? GO.union(alg, A, B; target = single) :
                   dim == 1 ? GO.intersection(alg, A, Cedge; target = single) :
                              GO.intersection(alg, GI.MultiPoint([(1.0, 1.0)]),
                                              GI.MultiPoint([(1.0, 1.0)]); target = single)
        empties = (GO.intersection(alg, A, Far; target = single),
                   GO.difference(alg, A, A; target = single))
        @test !isempty(nonempty)
        for e in empties
            @test isempty(e)
            @test typeof(e) === typeof(nonempty)
        end

        nonempty_m = dim == 2 ? GO.union(alg, A, B; target = multi) :
                     dim == 1 ? GO.intersection(alg, A, Cedge; target = multi) :
                                GO.intersection(alg, GI.MultiPoint([(1.0, 1.0)]),
                                                GI.MultiPoint([(1.0, 1.0)]); target = multi)
        for e in (GO.intersection(alg, A, Far; target = multi),
                  GO.difference(alg, A, A; target = multi))
            @test GI.ngeom(e) == 0
            @test typeof(e) === typeof(nonempty_m)
        end
    end
end

@testset "target above the result dimension short-circuits before noding" begin
    A = SQ_A
    L = GI.LineString([(-1.0, 1.0), (3.0, 1.0)])
    P = GI.MultiPoint([(0.5, 0.5)])
    alg = GO.OverlayNG()

    #-- OGC result dimensions: min for intersection, max for union/symdiff, lhs
    #-- for difference. A target above that is unsatisfiable for ANY input. One
    #-- row per op — the op is the only axis `_result_dimension` reads.
    unsat = [(GO.intersection,  L, A, GI.PolygonTrait()),
             (GO.difference,    P, A, GI.MultiLineStringTrait()),
             (GO.union,         L, L, GI.PolygonTrait()),
             (GO.symdifference, P, P, GI.MultiLineStringTrait())]
    for (f, X, Y, t) in unsat
        r = f(alg, X, Y; target = t)
        @test r isa AbstractVector ? isempty(r) : GI.ngeom(r) == 0
    end

    #=
    That the answer is empty does NOT show the short circuit ran: without it the
    pipeline reaches the same empty answer. Two probes that do discriminate.

    First, an input the noder cannot process at all. A spherical edge between
    antipodal vertices has no defined great circle, so noding throws — and an
    above-dimension target must still return empty, which is only possible if
    noding never happened.
    =#
    Ant = GI.LineString([(0.0, 0.0), (180.0, 0.0)])
    SA  = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
    sph = GO.OverlayNG(Spherical())
    @test_throws ArgumentError GO.intersection(sph, Ant, SA)
    @test isempty(GO.intersection(sph, Ant, SA; target = GI.PolygonTrait()))

    #-- second, the cost: the short circuit allocates essentially nothing, where
    #-- the same pair through the pipeline allocates ~21 kB
    GO.intersection(alg, L, A; target = GI.PolygonTrait())   # warm up
    @test @allocated(GO.intersection(alg, L, A; target = GI.PolygonTrait())) < 1024

    #-- and it does NOT fire when the target is at or below the result dimension
    @test length(GO.intersection(alg, L, A; target = GI.LineStringTrait())) >= 1
end

@testset "target elides work it cannot want" begin
    A = SQ_A
    alg = GO.OverlayNG()
    mpoly, mpt = GI.MultiPolygonTrait(), GI.MultiPointTrait()
    #-- 5k points against one polygon: locating them is the whole cost of this
    #-- union, and an areal target must not pay it
    Pts = GI.MultiPoint([(4 * (i / 5000) - 1, 4 * ((i * 7919) % 5000) / 5000 - 1)
                         for i in 1:5000])

    #-- the elision is structural, not a timing accident: an areal target builds
    #-- neither the locator nor the point map, and gets the shared empty list back
    #-- by identity; a point target skips the `tuples` copy of the non-point input
    @test GO._mixed_points(Planar(), Pts, A, 2, true, mpoly; exact = EX) === GO._NO_POINTS
    let (pc, lc) = GO._mixed_components(A, 2, mpt)
        @test pc === GO._NO_COMPONENTS && lc === GO._NO_COMPONENTS
    end

    #-- end to end, in allocations rather than wall clock (CI runs under
    #-- `--code-coverage`, which taxes timing unevenly — see the perf notes).
    #-- Measured margin is ~3800×, so the ÷10 gate is nowhere near the noise.
    GO.union(alg, Pts, A); GO.union(alg, Pts, A; target = mpoly)   # warm up
    @test @allocated(GO.union(alg, Pts, A; target = mpoly)) <
          @allocated(GO.union(alg, Pts, A)) ÷ 10

    #-- and the answer is still right: A survives the union untouched. (That the
    #-- areal answer matches the untargeted one for every op is the sweep's job.)
    @test GI.ngeom(GO.union(alg, Pts, A; target = mpoly)) == 1
    @test GO.area(GO.union(alg, Pts, A; target = mpoly)) == GO.area(A)
end

@testset "target on the sphere, and the full-sphere gate" begin
    A = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
    B = GI.Polygon([[(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 15.0), (5.0, 5.0)]])
    sph = GO.OverlayNG(Spherical())
    #-- targeting is manifold-agnostic: same area as the untargeted result
    for op in OPS
        plain = GO._overlay_ng(Spherical(), op, A, B; exact = EX)
        tgt = GO._overlay_ng(Spherical(), op, A, B; exact = EX,
                             target = GI.MultiPolygonTrait())
        @test isapprox(GO.area(Spherical(), tgt), GO.area(Spherical(), plain); rtol = 1e-12)
    end
    @test GI.ngeom(GO.intersection(sph, A, B; target = GI.MultiPolygonTrait())) == 1

    #-- the full-sphere rejection is areal, so it fires for an areal target (and
    #-- for no target) and is skipped for a target that excludes areas
    #-- built by hand because no VALID pair of inputs can reach the gate: under
    #-- enclosed-region semantics a boundaryless full-sphere union forces both
    #-- operands to be exact hemispheres, whose great-circle ring is degenerate
    #-- here. (The untargeted throw is also asserted at the §4 testset above.)
    inp = GO._OverlayInput(Spherical(), A, A, 2, 2, EX, false, false, nothing, nothing)
    for t in (GI.PolygonTrait(), GI.MultiPolygonTrait())
        @test_throws ArgumentError GO._resolve_empty_result(Spherical(), GO.OVERLAY_UNION,
                                                            inp, t)
    end
    for t in (GI.LineStringTrait(), GI.MultiLineStringTrait(),
              GI.PointTrait(), GI.MultiPointTrait())
        r = GO._resolve_empty_result(Spherical(), GO.OVERLAY_UNION, inp, t)
        @test r isa AbstractVector ? isempty(r) : GI.ngeom(r) == 0
    end
end

@testset "target accepts the Foster–Hormann spellings and rejects the rest" begin
    A, B = SQ_A, SQ_B
    alg = GO.OverlayNG()
    want = GO.intersection(alg, A, B; target = GI.PolygonTrait())
    #-- trait instance, trait type, and TraitTarget all mean the same thing
    for t in (GI.PolygonTrait(), GI.PolygonTrait, GO.TraitTarget(GI.PolygonTrait()),
              GO.TraitTarget(GI.PolygonTrait))
        r = GO.intersection(alg, A, B; target = t)
        @test typeof(r) === typeof(want)
        @test coordlist(r) == coordlist(want)
    end
    #-- no target is the default, and stays the untargeted geometry
    @test GI.trait(GO.intersection(alg, A, B; target = nothing)) isa GI.PolygonTrait

    #-- everything else is rejected, including a Union of traits (no single
    #-- result dimension) and the traits this engine cannot emit
    for bad in (GI.LinearRingTrait(), GI.GeometryCollectionTrait(), GI.FeatureTrait(),
                GO.TraitTarget(GI.PolygonTrait(), GI.PointTrait()),
                GO.TraitTarget{Union{GI.PolygonTrait, GI.LineStringTrait}}(), 42, "polygon")
        @test_throws ArgumentError GO.intersection(alg, A, B; target = bad)
    end
end

#=
Type stability is the point of `target`, so it is asserted rather than assumed.
For inputs with edges the untargeted return type is a 7-member `Union` — `Point`,
`LineString`, `MultiLineString`, `MultiPoint`, `Polygon`, `MultiPolygon`,
`GeometryCollection` — because the engine returns the most specific geometry over
whatever survived, which inference cannot know. (For a `MultiPoint` operand it is
`Any`.) Targeted, every combination below infers to one concrete type.

`Base.return_types` over a closure is the instrument, not `@inferred`: it asks
what the compiler can prove at a call site where the target is a compile-time
constant (which is what `target = GI.PolygonTrait()` is), independent of what
this particular pair of inputs happens to produce.

The axes are chosen, not crossed exhaustively. The INPUT SHAPE matters, because
inference const-folds the driver's dimension branches, so the three shapes below
reach three different sets of return sites (arrangement / mixed-point /
point×point). The OP does not — every return site funnels through
`_target_result(target, …)` and `op` is only ever compared to a constant — so it
is probed once, at the level of the four public wrappers, which are four separate
methods each forwarding `target` and so the one place an op-specific inference
regression could hide.
=#
@testset "target makes the return type concrete and input-independent" begin
    A  = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B  = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    P  = GI.MultiPoint([(0.5, 0.5), (5.0, 5.0)])
    argtypes = [Tuple{typeof(A), typeof(B)},   # the arrangement
                Tuple{typeof(P), typeof(A)},   # the mixed-point path
                Tuple{typeof(P), typeof(P)}]   # the point × point path

    for (single, multi, _) in TARGETS, t in (single, multi)
        inferred = Type[]
        for m in (Planar(), Spherical()), Ts in argtypes
            probe = (a, b) -> GO.intersection(GO.OverlayNG(m), a, b; target = t)
            R = Base.return_types(probe, Ts)[1]
            @test isconcretetype(R)
            push!(inferred, R)
        end
        #-- and it is ONE type across every manifold and input shape
        @test length(unique(inferred)) == 1
        #-- which is exactly the type the call actually returns
        @test only(unique(inferred)) === typeof(GO.intersection(GO.OverlayNG(), A, B; target = t))
    end

    #-- all four public wrappers forward `target` without widening the inference
    for f in (GO.intersection, GO.union, GO.difference, GO.symdifference)
        probe = (a, b) -> f(GO.OverlayNG(), a, b; target = GI.PolygonTrait())
        @test isconcretetype(Base.return_types(probe, Tuple{typeof(A), typeof(B)})[1])
    end

    #-- the contrast: untargeted, inference can only give the 7-member Union of
    #-- everything the engine is allowed to emit
    untargeted = Base.return_types((a, b) -> GO.intersection(GO.OverlayNG(), a, b),
                                   Tuple{typeof(A), typeof(B)})[1]
    @test untargeted isa Union
    @test !isconcretetype(untargeted)
    #-- but every MEMBER of it is concrete, which is only true because the
    #-- GeometryCollection branch builds a `Vector{_ResultComponent}` rather than
    #-- a `Vector{Any}`: with `Any` elements the wrapper's hasz/hasm detection
    #-- cannot fold either, and that member degrades to a `where {_A,_B}`
    members(U) = U isa Union ? vcat(members(U.a), members(U.b)) : Any[U]
    ms = members(untargeted)
    @test length(ms) == 7
    @test all(isconcretetype, ms)
end

@testset "a mixed-dimension result is a typed collection, not Vector{Any}" begin
    gc = GO.intersection(GO.OverlayNG(), MIXED_A, MIXED_B)
    @test GI.trait(gc) isa GI.GeometryCollectionTrait
    #-- the components iterate as the engine's three-way Union, so a caller
    #-- walking them union-splits instead of dispatching dynamically per element
    #-- this is the ONLY assertion here that catches a regression to `Vector{Any}`;
    #-- that the wrapper TYPE infers concretely is asserted by `all(isconcretetype,
    #-- ms)` in the testset above (`isconcretetype(typeof(x))` would be a tautology
    #-- — `typeof` is always concrete — so it is deliberately not written here)
    @test eltype(GI.getgeom(gc)) === GO._ResultComponent
    #-- and the collection really is mixed-dimension, with its parts intact
    traits = [GI.trait(g) for g in GI.getgeom(gc)]
    @test any(t -> t isa GI.PolygonTrait, traits)
    @test any(t -> t isa GI.LineStringTrait, traits)
end

@testset "target on all four public operations, including symdifference" begin
    A, B = SQ_A, SQ_B
    alg = GO.OverlayNG()
    #-- the SINGULAR target, deliberately: it returns a `Vector`, which no
    #-- untargeted call can produce, so these fail if a wrapper drops the keyword.
    #-- (`MultiPolygonTrait` would not: the untargeted symdifference of this pair
    #-- is already a MultiPolygon, and comparing areas cannot tell the two apart
    #-- for any op, since an ignored target returns the same geometry.)
    for f in (GO.intersection, GO.union, GO.difference, GO.symdifference)
        v = f(alg, A, B; target = GI.PolygonTrait())
        @test v isa AbstractVector
        @test all(g -> GI.trait(g) isa GI.PolygonTrait, v)
    end
    #-- symdifference's algorithm-free and manifold forms take it too
    @test GO.symdifference(A, B; target = GI.PolygonTrait()) isa AbstractVector
    @test GO.symdifference(Planar(), A, B; target = GI.PolygonTrait()) isa AbstractVector
    @test GO.symdifference(Spherical(), A, B; target = GI.PolygonTrait()) isa AbstractVector
end

# ---------------------------------------------------------------------------
# ## Sub-grid ring collapse (`_ring_is_subgrid`, the second half of
# ## `_ring_is_collapsed`)
#
# Two inputs that share a boundary can disagree about it by an ULP or two — a
# cascade of unions manufactures exactly that, because every intermediate result
# is re-ingested through Float64 coordinates. The exact arrangement then contains
# genuine needle-shaped faces between the two versions of the edge, and they are
# real: the operands really do differ there. What they do not have is a faithful
# Float64 image. Emission rounds each vertex by up to half an ULP, which is
# enough to push a needle's two sides through each other, and the result is a
# self-touching or self-crossing "ring" — invalid output produced from a
# perfectly correct exact face.
#
# `_ring_is_subgrid` drops such a ring, deciding on the exact node kernel points
# and comparing against the resolution of the OUTPUT FORMAT where the ring sits
# (`_ring_grid_step`). See `maximal_edge_ring.jl` for why that is a format
# property rather than a tolerance.
#
# The fixture: two rectangles sharing a subdivided vertical edge at lon = 1,
# where B's copy of the shared vertices alternates `k` ULPs either side of A's.
# Each right-hand excursion opens a needle-shaped gap that the union must either
# represent faithfully or drop.
# ---------------------------------------------------------------------------

_ulp(x, k) = k >= 0 ? nextfloat(x, k) : prevfloat(x, -k)

function shared_edge_pair(nsub::Int, k::Int; x0 = 1.0, y0 = 49.0, y1 = 50.0)
    ys = collect(range(y0, y1; length = nsub + 2))
    A = GI.Polygon([GI.LinearRing(vcat([(x0 - 1, y0)], [(x0, y) for y in ys],
                                       [(x0 - 1, y1)], [(x0 - 1, y0)]))])
    #-- endpoints stay shared exactly; interior vertices alternate -k/+k ULPs
    pert = [(i == 1 || i == length(ys)) ? (x0, y) :
            (_ulp(x0, isodd(i) ? k : -k), y) for (i, y) in enumerate(ys)]
    B = GI.Polygon([GI.LinearRing(vcat([(x0, y0), (x0 + 1, y0), (x0 + 1, y1)],
                                       reverse(pert)))])
    return (A, B)
end

@testset "a sub-grid sliver is dropped, and the union stays valid" begin
    for nsub in (3, 5, 9)
        A, B = shared_edge_pair(nsub, 1)
        for m in (Planar(), Spherical())
            r = GO.union(GO.OverlayNG(m), A, B)
            #-- every needle is finer than the output grid, so none survives and
            #-- the union is the plain 2 x 1 rectangle
            @test GI.nring(r) == 1
            @test LG.isValid(lgc(r))
        end
        #-- the planar area is exactly the rectangle's: the dropped needles are
        #-- below Float64 resolution, so they cannot show up in it
        @test GO.area(GO.union(GO.OverlayNG(), A, B)) == 2.0
    end
end

#=
The failure mode of a collapse test is dropping geometry that IS representable,
so this is the load-bearing half of the pair: widen the same needles until the
output format can hold them, and every one must come back.

`nsub` needles at `k` ULPs, and `k = 16` is 4x the threshold — see the boundary
test below for where it actually sits.
=#
@testset "a representable sliver is kept" begin
    for (nsub, nholes) in ((3, 1), (9, 4))
        A, B = shared_edge_pair(nsub, 16)
        for m in (Planar(), Spherical())
            r = GO.union(GO.OverlayNG(m), A, B)
            @test GI.nring(r) == 1 + nholes
            @test LG.isValid(lgc(r))
        end
        #-- and each one is a legal ring of representable width, not a sliver of
        #-- nothing. (Its AREA is deliberately not asserted: at 16 ULPs of width
        #-- on coordinates of magnitude 1 it is ~1e-15, which a Float64 shoelace
        #-- over terms of magnitude 50 cannot resolve — `GO.area` returns exactly
        #-- 2.0 for the whole polygon either way. That the area is invisible to
        #-- Float64 while the ring is perfectly representable is the distinction
        #-- this whole testset turns on, so measuring the wrong one here would be
        #-- worse than measuring nothing.)
        r = GO.union(GO.OverlayNG(), A, B)
        for h in GI.gethole(r)
            pts = collect(GI.getpoint(h))
            @test length(unique(pts)) >= 3
            @test maximum(GI.x, pts) > minimum(GI.x, pts)
        end
    end
end

#=
The threshold is `_RING_GRID_MARGIN` steps of the output grid and nothing else,
so it must land in the same place on both manifolds and move only with the
format's resolution — not with the geometry, the manifold, or the operand size.

At `x0 = 1.0` the grid step is `eps(1.0) = 2^-52` (latitude's `eps(50.0)` is
coarser, and `_ring_grid_step` takes the finer axis). A `k`-ULP excursion makes a
needle of maximum width `k` steps and hence MEAN width `k/2` steps, so the
`< 4` test flips between `k = 6` (mean width 3) and `k = 8` (mean width 4).
Both manifolds, same k, because both are measuring the same output grid.
=#
@testset "the collapse threshold is the output grid, and nothing else" begin
    @test GO._RING_GRID_MARGIN == 4.0
    for m in (Planar(), Spherical())
        @test GI.nring(GO.union(GO.OverlayNG(m), shared_edge_pair(3, 6)...)) == 1
        @test GI.nring(GO.union(GO.OverlayNG(m), shared_edge_pair(3, 8)...)) == 2
    end

    #-- the grid step itself: the finer of the two axes, and on the sphere
    #-- converted from degrees to the radians the exact width is measured in
    pts = [(1.0, 49.0), (1.0, 50.0), (2.0, 50.0), (1.0, 49.0)]
    @test GO._ring_grid_step(Planar(), pts) == min(eps(2.0), eps(50.0))
    @test GO._ring_grid_step(Spherical(), pts) ≈
          min(eps(2.0) * cosd(49.5), eps(50.0)) * (π / 180)
    #-- it is genuinely LOCAL: the same shape a thousand times further out sits
    #-- on a grid a thousand times coarser, so the test is scale-free
    far = [(1024.0 * p[1], 1024.0 * p[2]) for p in pts]
    @test GO._ring_grid_step(Planar(), far) == 1024 * GO._ring_grid_step(Planar(), pts)
end

#=
Half 1 (`_ring_image_is_degenerate`, fewer than three distinct emitted vertices)
and half 2 are independent: neither implies the other, and the spike population
that motivated half 2 is invisible to half 1. Measured on the 1384-polygon
Vancouver watershed cascade, half 1 fires on ZERO of the 2974 minimal rings built
while half 2 fires on 10 — so a regression that quietly reduced the pair to half
1 would not show up as an error anywhere, only as invalid output.
=#
@testset "half 1 and half 2 of the collapse test are independent" begin
    A, B = shared_edge_pair(3, 1)
    m = Spherical()
    arr = GO.NodedArrangement(m, A, B; exact = EX)
    g = GO.OverlayGraph(m, arr; exact = EX)
    GO._compute_labelling!(g, GO._OverlayInput(m, A, B, 2, 2, EX,
                                               false, false, nothing, nothing))
    ctx = GO._build_faces(m, g; exact = EX)
    subgrid = [r for r in ctx.edge_rings if GO._ring_is_subgrid(ctx, r)]
    @test !isempty(subgrid)
    #-- the needles have a perfectly ordinary-looking emitted image: three or
    #-- more distinct vertices, so half 1 sees nothing wrong with any of them
    @test all(r -> !GO._ring_image_is_degenerate(r), subgrid)
    #-- but the full test drops them, and the one-argument method does not
    @test all(r -> GO._ring_is_collapsed(ctx, r), subgrid)
    @test all(r -> !GO._ring_is_collapsed(r), subgrid)
end

#=
Converting (lon, lat) to xyz rounds `cos(φ)cos(λ)` and `cos(φ)sin(λ)`
independently, so it does NOT preserve coplanarity along a meridian: two lat/lon
cells that share the meridian `λ = 8` do not share a great circle, they get two
that differ by ~1e-17 rad. The arrangement is exact, so it sees the gap and rings
it — legitimately. Every consumer that assumed a shared grid edge is *exactly*
shared broke on that, and this is the last of them.

The sliver's four nodes are two vertex nodes and two crossings a fraction of an
ULP from them. Rounding the crossings' exact directions to `UnitSphericalPoint`
(`_node_kernel_point`) lands one of them bit-for-bit on a vertex node and the
other one ULP away, so `_prune_loop_degeneracies` deletes the coincident pair as
a repeated vertex and the Float64 curvature answers about a degenerate triangle
instead. It reported the ring as a hole; exactly, the ring is clockwise — a
shell. Nothing then contained the "hole" and the overlay threw `unable to assign
free hole to a shell`.

Emitting a ~1e-17-rad-wide sliver is the right answer for the arrangement as
given, and is not what this tests. What it tests is that the ring's ROLE is read
off the exact node directions when the rounded ones cannot carry it.
=#
@testset "a sub-ULP result ring is oriented from exact node directions" begin
    m = Spherical()
    #-- a 1°×1° cell against a polar triangle; they meet along λ = 8. The
    #-- triangle repeats its pole vertex, as polar cells of real grids do.
    A = GI.Polygon([[(8.0, -89.0), (9.0, -89.0), (9.0, -88.0), (8.0, -88.0), (8.0, -89.0)]])
    B = GI.Polygon([[(0.0, -90.0), (0.0, -90.0), (8.0, -86.0), (4.0, -86.0), (0.0, -90.0)]])

    #-- end to end: this used to throw. The regions meet along an edge and
    #-- nowhere else, so the intersection is empty up to the meridian gap.
    r = GO._overlay_ng(m, GO.OVERLAY_INTERSECTION, A, B; exact = EX)
    @test GO.area(m, r) <= 1e-12 * GO.area(m, A)

    #-- and on the ring itself. Rebuild the one minimal ring the result has.
    arr = GO.NodedArrangement(m, A, B; exact = EX)
    g = GO.OverlayGraph(m, arr; exact = EX)
    GO._compute_labelling!(g, GO._OverlayInput(m, A, B, 2, 2, EX,
                                               false, false, nothing, nothing))
    GO._mark_result_area_edges!(g, GO.OVERLAY_INTERSECTION)
    GO._unmark_duplicate_edges_from_result_area!(g)
    rae = GO.graph_result_area_edges(g)
    P = typeof(GO._to_kernel_point(m, (0.0, 0.0)))
    ctx = GO._PolyBuilderCtx(m, g.edges, g.arr, EX, GO._MaxEdgeRing[],
                             GO._OverlayEdgeRing{P}[], Int32[], Int32[])
    for e in rae
        GO._link_result_area_max_ring_at_node!(ctx.edges, e)
    end
    mr = only(GO._build_maximal_rings!(ctx, rae))
    ids = ctx.edge_rings[only(GO._build_minimal_rings!(ctx, ctx.max_rings[mr]))].node_ids

    ks = [GO._node_kernel_point(ctx, i) for i in ids]
    dirs = [GO._node_exact_dir(ctx, i) for i in ids]
    #-- the ring is not readable at Float64: two of its four nodes round together
    @test !GO._rounded_ring_is_faithful(ks)
    @test length(unique(ks)) < length(unique(ids))
    #-- so the two tests disagree, and the exact one — every consecutive triple
    #-- turning the same way, hence convex and clockwise — is the shell
    @test GO._ring_is_ccw(m, ks; exact = EX) != GO._ring_is_ccw_dirs(dirs)
    @test all(i -> GO._dot3(GO._cross3(dirs[mod1(i - 1, 4)], dirs[i]),
                            dirs[mod1(i + 1, 4)]) < 0, 1:4)
    @test !GO._ring_is_ccw_dirs(dirs)
end

#=
The escalation above must not change any answer the Float64 path was entitled to
give, so the two agree wherever the rounded ring is readable at all — which is
every ring of an ordinary overlay.
=#
@testset "exact and rounded ring orientation agree on readable rings" begin
    m = Spherical()
    A = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
                    [(2.0, 2.0), (2.0, 5.0), (5.0, 5.0), (5.0, 2.0), (2.0, 2.0)]])
    B = GI.Polygon([[(3.0, 3.0), (13.0, 3.0), (13.0, 13.0), (3.0, 13.0), (3.0, 3.0)]])
    checked = 0
    for op in OPS
        arr = GO.NodedArrangement(m, A, B; exact = EX)
        g = GO.OverlayGraph(m, arr; exact = EX)
        GO._compute_labelling!(g, GO._OverlayInput(m, A, B, 2, 2, EX,
                                                   false, false, nothing, nothing))
        GO._mark_result_area_edges!(g, op)
        GO._unmark_duplicate_edges_from_result_area!(g)
        rae = GO.graph_result_area_edges(g)
        isempty(rae) && continue
        P = typeof(GO._to_kernel_point(m, (0.0, 0.0)))
        ctx = GO._PolyBuilderCtx(m, g.edges, g.arr, EX, GO._MaxEdgeRing[],
                                 GO._OverlayEdgeRing{P}[], Int32[], Int32[])
        for e in rae
            GO._link_result_area_max_ring_at_node!(ctx.edges, e)
        end
        for mr in GO._build_maximal_rings!(ctx, rae)
            for er in GO._build_minimal_rings!(ctx, ctx.max_rings[mr])
                ids = ctx.edge_rings[er].node_ids
                ks = [GO._node_kernel_point(ctx, i) for i in ids]
                #-- fixture premise: these rings are nowhere near the threshold
                @test GO._rounded_ring_is_faithful(ks)
                @test GO._ring_is_ccw(m, ks; exact = EX) ==
                      GO._ring_is_ccw_dirs([GO._node_exact_dir(ctx, i) for i in ids])
                checked += 1
            end
        end
    end
    @test checked >= 4
end
