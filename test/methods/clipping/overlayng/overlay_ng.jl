# Tests for the OverlayNG phase-2b engine core (design §3): the labeller, the
# result builders (polygons with the real minimal-ring split + hole nesting,
# lines, points), and the internal `_overlay_ng` driver — end to end, over the
# phase-1 arrangement and the phase-2a graph.
#
# Equality strategy: planar results are checked against LibGEOS's own overlay
# (high-level API only) with GEOS topological `equals` (order/orientation/merge-
# granularity independent), `isValid`, and area agreement at rtol 1e-12 (compared
# through the same `GO.area`, so it is a machine-precision gate). Spherical is
# gated on area conservation. Ported JTS cases keep JTS's expected WKT answers,
# compared via GEOS `equals`.

using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG
import GeometryOps: Planar, Spherical, True

const EX = True()

lgc(g) = GI.convert(LG, g)
giwkt(wkt) = GO.tuples(LG.readgeom(wkt))

const OPS = (GO.OVERLAY_INTERSECTION, GO.OVERLAY_UNION,
             GO.OVERLAY_DIFFERENCE, GO.OVERLAY_SYMDIFFERENCE)
opname(op) = op == GO.OVERLAY_INTERSECTION ? "intersection" :
             op == GO.OVERLAY_UNION ? "union" :
             op == GO.OVERLAY_DIFFERENCE ? "difference" : "symdifference"

geos_op(op, A, B) =
    op == GO.OVERLAY_INTERSECTION ? LG.intersection(lgc(A), lgc(B)) :
    op == GO.OVERLAY_UNION ? LG.union(lgc(A), lgc(B)) :
    op == GO.OVERLAY_DIFFERENCE ? LG.difference(lgc(A), lgc(B)) :
    LG.symmetricDifference(lgc(A), lgc(B))

# Planar: check `_overlay_ng` against LibGEOS for one op.
function check_planar(op, A, B; areatol = 1e-12)
    r = GO._overlay_ng(Planar(), op, A, B; exact = EX)
    geos = geos_op(op, A, B)
    @test LG.isValid(lgc(r))
    @test LG.equals(lgc(r), geos)
    @test isapprox(GO.area(Planar(), r), GO.area(Planar(), GO.tuples(geos));
                   rtol = areatol, atol = 1e-12)
    return r
end

check_all_ops(A, B; areatol = 1e-12) =
    for op in OPS
        @testset "$(opname(op))" begin check_planar(op, A, B; areatol) end
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
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    #-- analytic areas: overlap [1,2]² = 1
    @test isapprox(GO.area(GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, A, B; exact = EX)), 1.0; rtol = 1e-12)
    @test isapprox(GO.area(GO._overlay_ng(Planar(), GO.OVERLAY_UNION, A, B; exact = EX)), 7.0; rtol = 1e-12)
    @test isapprox(GO.area(GO._overlay_ng(Planar(), GO.OVERLAY_DIFFERENCE, A, B; exact = EX)), 3.0; rtol = 1e-12)
    @test isapprox(GO.area(GO._overlay_ng(Planar(), GO.OVERLAY_SYMDIFFERENCE, A, B; exact = EX)), 6.0; rtol = 1e-12)
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
    #-- A area = 100 - 16 = 84; overlap of B with A-material = B∩A = 21
    @test isapprox(GO.area(GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, A, B; exact = EX)), 21.0; rtol = 1e-12)
    check_all_ops(A, B)
end

@testset "collinear shared boundary (degenerate intersection, merged union)" begin
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(2.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 0.0)]])
    #-- intersection is the shared boundary line (1-D), not an area
    ri = GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, A, B; exact = EX)
    @test GI.trait(ri) isa GI.LineStringTrait
    @test LG.equals(lgc(ri), geos_op(GO.OVERLAY_INTERSECTION, A, B))
    #-- union merges into one 2×4 box
    ru = GO._overlay_ng(Planar(), GO.OVERLAY_UNION, A, B; exact = EX)
    @test GI.trait(ru) isa GI.PolygonTrait
    @test isapprox(GO.area(ru), 8.0; rtol = 1e-12)
    @test LG.isValid(lgc(ru)) && LG.equals(lgc(ru), geos_op(GO.OVERLAY_UNION, A, B))
end

@testset "degree-6 coincident crossing (all ops)" begin
    #-- two A squares touching at (1,1) + a B triangle edge through (1,1)
    A = GI.MultiPolygon([[[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]],
                         [[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)]]])
    B = GI.Polygon([[(0.0, 0.0), (2.0, 2.0), (2.0, 0.0), (0.0, 0.0)]])
    arr = GO.NodedArrangement(Planar(), A, B; exact = EX)
    g = GO.OverlayGraph(Planar(), arr; exact = EX)
    nid = findfirst(1:GO.num_nodes(arr)) do i
        p = GO.node_point(arr, i)
        isapprox(p[1], 1.0; atol = 1e-9) && isapprox(p[2], 1.0; atol = 1e-9)
    end
    @test nid !== nothing
    e = GO.graph_node_edge(g, nid)
    @test e != 0
    @test GO.he_degree(g.edges, e) == 6      # 2 (square 1) + 2 (square 2) + 2 (B diagonal)
    check_all_ops(A, B)
end

@testset "concave L-shapes (all ops)" begin
    A = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (1.0, 1.0), (1.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
    B = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (2.0, 3.0), (2.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
    check_all_ops(A, B)
end

@testset "MultiPolygon input (all ops)" begin
    A = GI.MultiPolygon([[[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]],
                         [[(5.0, 5.0), (7.0, 5.0), (7.0, 7.0), (5.0, 7.0), (5.0, 5.0)]]])
    B = GI.Polygon([[(1.0, 1.0), (6.0, 1.0), (6.0, 6.0), (1.0, 6.0), (1.0, 1.0)]])
    check_all_ops(A, B)
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
# SKIPPED (substrate limitation, not precision): testTouchingPolyDifference — a
#   single input whose hole touches its own shell, where the DIFFERENCE splits the
#   result into point-touching polygons; the substrate does not self-node one
#   input (design §2.2), so it yields a correct-area but non-simple result.

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
end

@testset "union of a multi-island geometry (France-class nesting)" begin
    #-- 20 disjoint island squares, unioned with a shifted copy (overlaps) — the
    #-- many-component case the spike prototypes faked, exercising minimal-ring
    #-- split + hole nesting across many shells.
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

@testset "spherical overlapping quads — area conservation" begin
    A = GI.Polygon([[(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0), (0.0, 0.0)]])
    B = GI.Polygon([[(10.0, 10.0), (30.0, 10.0), (30.0, 30.0), (10.0, 30.0), (10.0, 10.0)]])
    ai = GO.area(Spherical(), GO._overlay_ng(Spherical(), GO.OVERLAY_INTERSECTION, A, B; exact = EX))
    au = GO.area(Spherical(), GO._overlay_ng(Spherical(), GO.OVERLAY_UNION, A, B; exact = EX))
    aA = GO.area(Spherical(), A); aB = GO.area(Spherical(), B)
    @test isapprox(au + ai, aA + aB; rtol = 1e-12)
    #-- difference + intersection reconstruct A
    ad = GO.area(Spherical(), GO._overlay_ng(Spherical(), GO.OVERLAY_DIFFERENCE, A, B; exact = EX))
    @test isapprox(ad + ai, aA; rtol = 1e-12)
    #-- symdifference = union - intersection
    asd = GO.area(Spherical(), GO._overlay_ng(Spherical(), GO.OVERLAY_SYMDIFFERENCE, A, B; exact = EX))
    @test isapprox(asd, au - ai; rtol = 1e-12)
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
    inp_union = GO._OverlayInput(Spherical(), A, A, 2, 2, EX, false, false, nothing, nothing)
    @test GO._covers_everything(Spherical(), GO.OVERLAY_UNION, inp_union)
    @test_throws ArgumentError GO._resolve_empty_result(Spherical(), GO.OVERLAY_UNION, inp_union)

    inp_int = GO._OverlayInput(Spherical(), A, Bdisjoint, 2, 2, EX, false, false, nothing, nothing)
    @test !GO._covers_everything(Spherical(), GO.OVERLAY_INTERSECTION, inp_int)
    rr = GO._resolve_empty_result(Spherical(), GO.OVERLAY_INTERSECTION, inp_int)
    @test GI.npoint(rr) == 0
end

# ---------------------------------------------------------------------------
# 5. Input validation + empty inputs
# ---------------------------------------------------------------------------

@testset "input validation and empty short-circuits" begin
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    #-- point inputs are rejected (phase 3)
    @test_throws ArgumentError GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, GI.Point((1.0, 1.0)), A; exact = EX)
    #-- disjoint intersection short-circuits to empty (planar envelope)
    Far = GI.Polygon([[(100.0, 100.0), (102.0, 100.0), (102.0, 102.0), (100.0, 102.0), (100.0, 100.0)]])
    r = GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, A, Far; exact = EX)
    @test GI.npoint(r) == 0
end

# ---------------------------------------------------------------------------
# 6. Spherical NE shifted-self smoke (env-gated, phase-1 smoke pattern)
# ---------------------------------------------------------------------------

ne_ok = false
ne_names = String[]; ne_geoms = Any[]
try
    import NaturalEarth, GeoJSON
    fc = NaturalEarth.naturalearth("admin_0_countries", 110)
    for f in fc
        gg = GeoJSON.geometry(f)
        (gg === nothing || GI.npoint(gg) == 0) && continue
        nm = try; string(f.NAME); catch; "?"; end
        push!(ne_names, nm); push!(ne_geoms, GO.tuples(gg))
    end
    global ne_ok = length(ne_geoms) > 0
catch err
    @info "Natural Earth subset skipped (data unavailable)" err
end

@testset "Natural Earth shifted-self area conservation (spherical + planar)" begin
    if !ne_ok
        @test_skip "Natural Earth data unavailable"
    else
        picks = String["Brazil", "France", "Egypt", "Australia"]
        tested = 0
        for nm in picks
            idx = findfirst(==(nm), ne_names)
            idx === nothing && continue
            A = ne_geoms[idx]
            LG.isValid(lgc(A)) || continue
            B = GO.apply(GI.PointTrait(), A) do p
                (GI.x(p) + 0.5, GI.y(p))
            end
            tested += 1
            for m in (Planar(), Spherical())
                ri = GO._overlay_ng(m, GO.OVERLAY_INTERSECTION, A, B; exact = EX)
                ru = GO._overlay_ng(m, GO.OVERLAY_UNION, A, B; exact = EX)
                aA = GO.area(m, A); aB = GO.area(m, B)
                @test isapprox(GO.area(m, ru) + GO.area(m, ri), aA + aB; rtol = 1e-9)
            end
        end
        @test tested >= 2
    end
end
