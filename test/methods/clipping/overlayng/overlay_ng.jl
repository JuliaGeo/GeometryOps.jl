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
    #-- from the crossing's own arcs. The emitted longitudes all sit in [0, 2].
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
    #-- point inputs are routed to the phase-3 point builders (they used to be
    #-- rejected here); see overlay_points.jl for their full coverage
    @test GI.trait(GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION, GI.Point((1.0, 1.0)), A;
                                  exact = EX)) isa GI.PointTrait
    #-- geometry collections are still rejected
    @test_throws ArgumentError GO._overlay_ng(Planar(), GO.OVERLAY_INTERSECTION,
        GI.GeometryCollection([GI.Point((1.0, 1.0))]), A; exact = EX)
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
const MIXED_A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
const MIXED_B = GI.MultiPolygon([
    GI.Polygon([[(1.0, 0.0), (3.0, 0.0), (3.0, 1.0), (1.0, 1.0), (1.0, 0.0)]]),   # overlaps
    GI.Polygon([[(-2.0, 0.0), (0.0, 0.0), (0.0, 2.0), (-2.0, 2.0), (-2.0, 0.0)]]) # touches
])

@testset "target is exactly the untargeted result, filtered by dimension" begin
    A  = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B  = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    #-- shares one edge only: the intersection is a pure line
    Cedge = GI.Polygon([[(2.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 0.0)]])
    #-- overlaps AND shares a boundary segment: a mixed-dimension result
    Dmixed = giwkt("POLYGON ((1 0, 4 0, 4 2, 1 2, 1 0))")
    Etouch = GI.Polygon([[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0), (2.0, 2.0)]])
    L  = GI.LineString([(-1.0, 1.0), (3.0, 1.0)])
    L2 = GI.LineString([(0.5, -1.0), (0.5, 3.0)])
    P  = GI.MultiPoint([(0.5, 0.5), (5.0, 5.0), (2.0, 2.0)])
    P2 = GI.MultiPoint([(0.5, 0.5), (9.0, 9.0)])

    pairs = [("area × area", A, B), ("shared edge", A, Cedge), ("overlap+boundary", A, Dmixed),
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
            vec  = GO._overlay_ng(Planar(), op, X, Y; exact = EX, target = single)
            mul  = GO._overlay_ng(Planar(), op, X, Y; exact = EX, target = multi)
            @test coordlist(targeted_parts(vec)) == want
            @test coordlist(targeted_parts(mul)) == want
            #-- the two container shapes agree with each other, and the Multi one
            #-- really is a Multi of the right dimension
            @test GI.trait(mul) === multi
            @test vec isa AbstractVector
        end
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
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    L = GI.LineString([(-1.0, 1.0), (3.0, 1.0)])
    P = GI.MultiPoint([(0.5, 0.5)])
    alg = GO.OverlayNG()

    #-- OGC result dimensions: min for intersection, max for union/symdiff, lhs
    #-- for difference. A target above that is unsatisfiable for ANY input.
    unsat = [(GO.intersection, L, A, GI.PolygonTrait()),
             (GO.intersection, A, L, GI.MultiPolygonTrait()),
             (GO.intersection, P, A, GI.LineStringTrait()),
             (GO.intersection, P, A, GI.PolygonTrait()),
             (GO.difference,   L, A, GI.PolygonTrait()),
             (GO.difference,   P, A, GI.MultiLineStringTrait()),
             (GO.union,        L, L, GI.PolygonTrait()),
             (GO.symdifference, P, P, GI.MultiLineStringTrait())]
    for (f, X, Y, t) in unsat
        r = f(alg, X, Y; target = t)
        @test r isa AbstractVector ? isempty(r) : GI.ngeom(r) == 0
    end

    #-- the short circuit reads the input dimensions only, so it survives inputs
    #-- the pipeline would reject outright
    @test isempty(GO.intersection(alg, GI.LineString([(0.0, 0.0), (1.0, 1.0)]), A;
                                  target = GI.PolygonTrait()))
    #-- and it does NOT fire when the target is at or below the result dimension
    @test length(GO.intersection(alg, L, A; target = GI.LineStringTrait())) >= 1
end

@testset "target elides work it cannot want" begin
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
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
    @test GO._mixed_points(Planar(), Pts, A, 2, true, mpt; exact = EX) !== GO._NO_POINTS
    let (pc, lc) = GO._mixed_components(A, 2, mpt)
        @test pc === GO._NO_COMPONENTS && lc === GO._NO_COMPONENTS
    end
    @test first(GO._mixed_components(A, 2, mpoly)) !== GO._NO_COMPONENTS

    #-- end to end, in allocations rather than wall clock (CI runs under
    #-- `--code-coverage`, which taxes timing unevenly — see the perf notes)
    GO.union(alg, Pts, A); GO.union(alg, Pts, A; target = mpoly)   # warm up
    @test @allocated(GO.union(alg, Pts, A; target = mpoly)) <
          @allocated(GO.union(alg, Pts, A)) ÷ 10

    #-- and the answer is still right: A survives the union untouched
    @test GI.ngeom(GO.union(alg, Pts, A; target = mpoly)) == 1
    @test GO.area(GO.union(alg, Pts, A; target = mpoly)) == GO.area(A)

    #-- the builder elision leaves the areal answer identical to the untargeted one
    B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    for op in OPS
        plain = GO._overlay_ng(Planar(), op, A, B; exact = EX)
        tgt = GO._overlay_ng(Planar(), op, A, B; exact = EX, target = GI.MultiPolygonTrait())
        @test isapprox(GO.area(tgt), GO.area(plain); rtol = 1e-12)
    end
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
    inp = GO._OverlayInput(Spherical(), A, A, 2, 2, EX, false, false, nothing, nothing)
    @test_throws ArgumentError GO._resolve_empty_result(Spherical(), GO.OVERLAY_UNION, inp)
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
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
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
Untargeted, the return type is a 7-member `Union` — `Point`, `LineString`,
`MultiLineString`, `MultiPoint`, `Polygon`, `MultiPolygon`, `GeometryCollection`
— because the engine returns the most specific geometry over whatever survived,
which inference cannot know. Targeted, every combination below infers to one
concrete type, and it is the same one for every input shape, op and manifold.

`Base.return_types` over a closure is the instrument, not `@inferred`: it asks
what the compiler can prove at a call site where the target is a compile-time
constant (which is what `target = GI.PolygonTrait()` is), independent of what
this particular pair of inputs happens to produce.
=#
@testset "target makes the return type concrete and input-independent" begin
    A  = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B  = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    MB = GI.MultiPolygon([B])
    L  = GI.LineString([(-1.0, 1.0), (3.0, 1.0)])
    P  = GI.MultiPoint([(0.5, 0.5), (5.0, 5.0)])
    argtypes = [Tuple{typeof(A), typeof(B)}, Tuple{typeof(A), typeof(MB)},
                Tuple{typeof(L), typeof(A)}, Tuple{typeof(P), typeof(A)}]

    for (single, multi, _) in TARGETS, t in (single, multi)
        inferred = Type[]
        for m in (Planar(), Spherical()),
            f in (GO.intersection, GO.union, GO.difference, GO.symdifference),
            Ts in argtypes
            probe = (a, b) -> f(GO.OverlayNG(m), a, b; target = t)
            R = Base.return_types(probe, Ts)[1]
            @test isconcretetype(R)
            push!(inferred, R)
        end
        #-- and it is ONE type across every op, manifold and input shape
        @test length(unique(inferred)) == 1
        #-- which is exactly the type the call actually returns
        @test only(unique(inferred)) === typeof(GO.intersection(GO.OverlayNG(), A, B; target = t))
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
    @test eltype(GI.getgeom(gc)) === GO._ResultComponent
    @test isconcretetype(typeof(gc))
    #-- and the collection really is mixed-dimension, with its parts intact
    traits = [GI.trait(g) for g in GI.getgeom(gc)]
    @test any(t -> t isa GI.PolygonTrait, traits)
    @test any(t -> t isa GI.LineStringTrait, traits)
    @test all(g -> g isa GO._ResultComponent, GI.getgeom(gc))
end

@testset "target on all four public operations, including symdifference" begin
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
    alg = GO.OverlayNG()
    for f in (GO.intersection, GO.union, GO.difference, GO.symdifference)
        r = f(alg, A, B; target = GI.MultiPolygonTrait())
        @test GI.trait(r) isa GI.MultiPolygonTrait
        @test isapprox(GO.area(r), GO.area(f(alg, A, B)); rtol = 1e-12)
    end
    #-- symdifference's algorithm-free and manifold forms take it too
    @test GI.trait(GO.symdifference(A, B; target = GI.MultiPolygonTrait())) isa GI.MultiPolygonTrait
    @test GI.trait(GO.symdifference(Planar(), A, B; target = GI.MultiPolygonTrait())) isa GI.MultiPolygonTrait
    @test isapprox(GO.area(Spherical(),
                           GO.symdifference(Spherical(), A, B; target = GI.MultiPolygonTrait())),
                   GO.area(Spherical(), GO.symdifference(Spherical(), A, B)); rtol = 1e-12)
end
